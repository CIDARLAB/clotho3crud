"""PCR predictor algorithm"""
import Bio.pairwise2
import re

#NOTE: Sanitize inputs as all uppercase for Biopython
Debug = False
EIPCR = True


def revComp(sequence):
    """Returns the reverse complement of a DNA sequence.

    sequence -- A string containing only 'A', 'T', 'C', 'G'
    """
    #Accepts degenerate codes
    #https://en.wikipedia.org/wiki/Nucleic_acid_notation
    baseComplement = {'a':'t','c':'g','t':'a','g':'c',
                    'u':'a',
                    'w':'w','s':'s','m':'k','k':'m','r':'y','y':'r',
                    'b':'v','d':'h','h':'d','v':'b',
                    'n':'n',' ':' ',
                    '^':'_', '_':'^', '|':'|'} #To not crash with sticky ends
    return ''.join([baseComplement[b] for b in sequence.lower()])[::-1]


def matchForward(seq1, seq2):
    """Returns a list of indices where end of HEXAMER of 3' end of Forward primer seq1 matches seq2.

    >>> seq1 = "aaaaaa"
    >>> seq2 = "gaaaaaagcggggcct"
    >>> matchForward(seq1, seq2)
    [7]
    """
    lenCheck = 6
    matchSeq = seq1[-lenCheck:]
    if EIPCR:
        temp = seq2+seq2[:lenCheck-1]
        return [(m.start()+lenCheck)%len(seq2) for m in re.finditer('(?='+matchSeq+')',temp)]
    else:
        return [(m.start()+lenCheck) for m in re.finditer('(?='+matchSeq+')',seq2)]

def matchReverse(seq1, seq2):
    """Returns a list of indices where HEXAMER of 3' end of Reverse primer matches seq2.

    >>> seq1 = "tttttt"
    >>> seq2 = "gggcccaaaaaa"
    >>> matchReverse(seq1, seq2)
    [5]
    """
    lenCheck = 6
    seq1 = revComp(seq1)
    matchSeq = seq1[0:lenCheck]
    if EIPCR:
        temp = seq2+seq2[:lenCheck-1]
        return [(m.start()-1)%len(seq2) for m in re.finditer('(?='+matchSeq+')',temp)]
    else:
        return [(m.start()-1) for m in re.finditer('(?='+matchSeq+')',seq2)]

#NOTE: Edge cases: oligo falls off end, pairwise align entire thing?
def calcBestForward(oligo, indexes, template):
    #FORWARD MATCHING CASE
    if Debug:
        print "CalcBestForward Running:"
    offset = len(oligo)
    oligoCheck = oligo[:-6]
    bestMatch = (0, 0) #Tuple of (index, alignmentScore)
    bestScore = -float("inf")
    for i in indexes:
        #Bounds check
        if len(template) <= i:
            break
        #Calculate pairwise alignment score of rest of oligo with rest of portion of template (up to length of oligo)
        leftOffset = i - offset
        if leftOffset < 0:
            if EIPCR:
                #Circular DNA so leftOffset is in "other" end of DNA
                leftOffset = len(template)+leftOffset
            else:
                leftOffset = 0
        if EIPCR and leftOffset > i: #Note: Add in bound errors if template too long
            templateCheck = (template[leftOffset:]+template[:i-6])[:len(oligoCheck)]
        else:
            templateCheck = template[leftOffset:i-6]
        score = Bio.pairwise2.align.globalms(oligoCheck, templateCheck, 2, -2, -1, -1, score_only = True)
        if Debug:
            print (i, score)
        if score > bestScore:
            bestMatch = (i, score)
            bestScore = score
    return bestMatch

def calcBestReverse(oligo, indexes, template):
    #REVERSE MATCHING CASE
    if Debug:
        print "CalcBestReverse Running:"
    oligo = revComp(oligo)
    offset = len(oligo)
    oligoCheck = oligo[6:]
    bestMatch = (0, 0) #Tuple of (index, alignmentScore)
    bestScore = -float("inf")
    for i in indexes:
        #Bounds check
        if len(template) <= i:
            break
        #Calculate pairwise alignment score of rest of oligo with rest of portion of template (up to length of oligo)
        rightOffset = i+1+offset
        if EIPCR and rightOffset > len(template):
            start = 0
            #Oligo hexamer extends over origin
            if (i+1+6) > len(template):
                start = (i+1+6)%len(template)
            templateCheck = template[i+1+6:rightOffset]+template[start:rightOffset%len(template)]
        else:
            templateCheck = template[i+1+6:rightOffset]
        score = Bio.pairwise2.align.globalms(oligoCheck, templateCheck, 2, -2, -1, -1, score_only = True)
        if Debug:
            print (i, score)
        if score > bestScore:
            bestMatch = (i, score)
            bestScore = score
    return bestMatch

#All functions have implicit assumption that inputs are 5' to 3'
def predict(oligo1, oligo2, template):
    """Takes two primers and a template and returns list of possible products.
    - This is done by simply running the algorithm twice
        - Case 1: Forward primer = oligo1, Reverse primer = oligo2
        - Case 2: Forward primer = oligo1, Reverse primer = oligo2, but use template = revComp(template)
    """
    oligo1 = oligo1.lower()
    oligo2 = oligo2.lower()
    template = template.lower()

    #CASE 1
    m1c1 = matchForward(oligo1, template)
    m2c1 = matchReverse(oligo2, template)
    best1F = calcBestForward(oligo1, m1c1, template)
    best1R = calcBestReverse(oligo2, m2c1, template)

    #CASE 2
    rc = revComp(template)
    m1c2 = matchForward(oligo1, rc)
    m2c2 = matchReverse(oligo2, rc)
    best2F = calcBestForward(oligo1, m1c2, rc)
    best2R = calcBestReverse(oligo2, m2c2, rc)

    #Currently best cut is highest sum of scores
    C1Score = best1F[1] + best1R[1]
    C2Score = best2F[1] + best2R[1]

    if Debug:
        print "m1c1: " + str(m1c1)
        print "m2c1: " + str(m2c1)
        print ("Best 1F: " + str(best1F) + " Best 1R: " + str(best1R))
        print "m1c2: " + str(m1c1)
        print "m2c2: " + str(m2c1)
        print ("Best 2F: " + str(best2F) + " Best 2R: " + str(best2R))

    #ADD IN NOTE BELOW
    #NOTE: If index of BestForward > index of BestReverse, WILL BE JUNK (SIGNAL ERROR, etc.)

    #NOTE: score = pairwise2.align.globalms(oligoCheck, templateCheck, 2, -2, -1, -1, score_only = True)
    #NOTE: Threshold for viable PCR product is score >= 30
    #threshold = 24
    #if C1Score < threshold and C2Score < threshold:
        #return "PCR Product score too low, less than threshold %s" % threshold
    #    return ""
    #Find best overall, check if product otherwise FAIL
    #Case 1 better
    if C1Score >= C2Score:
        forwardIndex1 = best1F[0]
        reverseIndex1 = best1R[0]
        retC1 = oligo1 + template[forwardIndex1: reverseIndex1+1] + revComp(oligo2)
        #EIPCR UNREFINED (Still need to test oligos over origin)
        if EIPCR:
            if forwardIndex1 > reverseIndex1:
                retC1 = oligo1 + template[forwardIndex1:] + template[:reverseIndex1+1] + revComp(oligo2)

        if Debug:
            print "C1Score (%s) >= C2Score (%s), so we choose Case 1 PCR Product" %(C1Score, C2Score)
            print retC1

        return retC1.replace("_^|","") #In case of stick end special characters (remove them)

    #Case 2 must be better
    forwardIndex2 = best2F[0]
    reverseIndex2 = best2R[0]
    retC2 = oligo1 + rc[forwardIndex2: reverseIndex2+1] + revComp(oligo2)
    #EIPCR UNREFINED (Still need to test oligos over origin)
    if EIPCR:
        if forwardIndex2 > reverseIndex2:
            retC2 = oligo1 + rc[forwardIndex2:] + rc[:reverseIndex2+1] + revComp(oligo2)
    if Debug:
        print "C1Score (%s) < C2Score (%s), so we choose Case 2 PCR Product" %(C1Score, C2Score)
        print retC2

    return retC2.replace("_^|","") #In case of stick end special characters (remove them)

    #NOTE: Add on warnings if scores are equal, etc.
    #   - EIPCR where input templates are short -> add error checking (i.e. calcBestForward loops around too much?)

def SOEing(oligo1, oligo2, *templates): #Takes in 2 oligos and arbitary number of templates
    pass
