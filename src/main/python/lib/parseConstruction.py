"""Parses inputted construction file into standardized Working Form .json representation
Input: String containing contents of construction file to be parsed
Output: json representation of input"""
import json
from collections import OrderedDict
import re
import sys

#Add in SOEing and PCA and gelPurify
#Note: Add in all supported reactions
#Note: What's in firstWords are currently supported reactions
firstWords = ("pcr", "eipcr", "digest", "ligate", "gelpurify")

def parse(fileString):
    steps = []
    dictionary = OrderedDict()
    contents = fileString.split("\n")
    for lineNum,line in list(enumerate(contents)):
        words = line.split()
        if len(words) == 0:
            continue
        #Identify function based on firstWord in each line
        #   - PCR = PCR reaction step
        #   - EIPCR = EIPCR reaction step
        #   - Digest = Digestion step
        #   - Ligate = Ligation step
        #   - > = Following characters are a name mapping to sequence in next line
        #       - Rest of words in current line are description for aforementioned name
        firstWord = words[0].lower()
        #If invalid format, skip line
        #Note: MUST HAVE firstWords BE ROBUST
        #   - Throw error or print something out to notify that this has occured
        if firstWord[0] != ">" and firstWord not in firstWords:
            continue

        #Need to worry about lengths, catching invalid, etc. For now assume valid
        #TODO: Distinguish between pcr and eipcr for pcr_predict input
        if firstWord == "pcr" or firstWord == "eipcr":
            oligos = words[1]

            if "/" not in oligos:
                #Throw error stating format is incorrect
                continue
            #Need to deal with following structures of oligo:
            #   PCR ca569/ca567R on pBAC562 (154 bp)            Simple case
            #   PCR ca845F/R on pHS4108     (814 bp, BstEII)    Oligos not explicitly listed
            #Assume not explicitly listed case first
            oligo1name = oligos[:-3] + 'F'
            oligo2name = oligos[:-3] + 'R'
            #Oligos explicitly listed case, so reassign
            if "f/r" not in oligos.lower():
                index = oligos.index("/")
                oligo1name = oligos[:index]
                oligo2name = oligos[index+1:]

            templateStartIndex = 0
            for index, word in list(enumerate(words)):
                if "on" == word.lower():
                    templateStartIndex = index + 1
                    break
            if templateStartIndex == 0:
                #Throw error stating format is incorrect
                continue
            templateEndIndex = 0
            for index, word in list(enumerate(words)):
                if "(" in word.lower():
                    templateEndIndex = index
                    break
            if templateEndIndex == 0:
                #Throw error stating format is incorrect
                continue
            templatename = ' '.join(words[i] for i in range(templateStartIndex,templateEndIndex))

            args = [w.strip('(),') for w in words[templateEndIndex:]]
            size = "-1"
            pcrProduct = ""

            for arg in args:
                if arg.isdigit():
                    productSize = arg
                else:
                    if "product" in arg.lower() or "pdt" in arg.lower():
                        pcrProduct = arg

            info = OrderedDict([("reaction", "pcr"),
                                ("input", [templatename,oligo1name,oligo2name]),
                                ("productSize", int(size)),
                                ("output", pcrProduct)])
            if size == "-1":
                del info["productSize"]
            steps.append(info)

        elif firstWord == "digest" or firstWord == "dig":
            templateEndIndex = 0
            for index, word in list(enumerate(words)):
                if "(" in word.lower():
                    templateEndIndex = index
                    break
            if templateEndIndex == 0:
                #Throw error stating format is incorrect
                continue
            templatename = ' '.join(words[i] for i in range(1,templateEndIndex))

            inputs = words[templateEndIndex:]

            sizes = ""
            choice = ""
            enzymes = ""
            digestProduct = ""
            inputs = [item.strip('(),') for item in inputs]
            for item in inputs:
                if "+" in item: #Figure out better way, i.e. regex re.findall(r'\d+', 'hello 42 I\'m a 32 string 30')
                    #Remove unncessary characters in string
                    sizes = item.replace("+"," ").split()
                elif len(item) == 1:
                    choice = item
                elif "/" in item or item[-1] == "I":
                    enzymes = item.split('/')
                else:
                    digestProduct = item
            if enzymes == "":
                inputs = [templatename]
            else:
                inputs = [templatename, enzymes]

            info = OrderedDict([("reaction", "digest"),
                ("input", inputs),
                ("sizes", sizes),
                ("choice", choice),
                ("output", digestProduct)])
            #Remove any fields not specified (note: some fields are mandatory)
            if sizes == "":
                del info["sizes"]
            if choice == "":
                del info["choice"]
            steps.append(info)

        elif firstWord == "ligate" or firstWord == "lig":
            inputsEndIndex = 0
            for index, word in list(enumerate(words)):
                if "(" in word.lower():
                    inputsEndIndex = index
                    break
            if inputsEndIndex == 0:
                #Throw error stating format is incorrect
                continue
            inputs = []
            for i in range(1,inputsEndIndex):
              inputs = inputs + words[i].split("/")

            #ligateProduct name is last word in line and surrounded by parantheses
            ligateProduct = words[-1][1:-1]

            info = OrderedDict([("reaction", "ligate"),
                                ("input", [inputs]),
                                ("output", ligateProduct)])
            steps.append(info)

        elif firstWord[0] == ">":
            key = words[0][1:]
            description = ' '.join(words[i] for i in range(1,len(words)))
            value = contents[lineNum+1].split()[0]
            dictionary[key] = (value, description)

        elif firstWord == "gelpurify":
            templatename = words[1]

            argsStartIndex = 0
            for index, word in list(enumerate(words)):
                if "(" in word.lower():
                    argsStartIndex = index
                    break
            if argsStartIndex == 0:
                #Throw error stating format is incorrect
                continue
            args = [w.strip('(),') for w in words[argsStartIndex:]]
            productSize = ""
            gelPurifyProduct = ""

            for arg in args:
                if arg.isdigit():
                    productSize = arg
                #elif "bp" in arg.lower(): Figure out better way, i.e. regex
                #REFINE SIZE PLACING
                #Note: Current Max supports only "L" (i.e. change "Largest" to "L")
                elif arg.lower() == "largest" or arg.lower() == "smallest" or arg.lower() == "l" or arg.lower() == "s":
                    productSize = arg
                elif "dig" in arg.lower():
                    gelPurifyProduct = arg

            info = OrderedDict([("reaction", "gelpurify"),
                                ("input", [templatename]),
                                ("length", productSize),
                                ("output", gelPurifyProduct)])
            #Remove any fields not specified (note: some fields are mandatory)
            if productSize == "":
                del info["length"]

            steps.append(info)

    dictContents = OrderedDict([(k, OrderedDict([("name",k),("description",v[1]),("sequence",v[0])])) for k,v in dictionary.iteritems()])
    #Remove Description field for those items with no description (empty string)
    for k,v in dictionary.iteritems():
    	if dictContents[k]["description"] == "":
    		del dictContents[k]["description"]

    data = OrderedDict([("schema", "org.clothocad.model.ConstructionFile"),
                        ("steps", steps),
                        ("dictionary", dictContents)])
    return data

if __name__ == '__main__':
    fileString = sys.argv[1]
    parse(fileString)
