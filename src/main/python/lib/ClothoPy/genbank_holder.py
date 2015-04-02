# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014
#
# For more information the Genbank file format, see:
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3245039/
# https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
# ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt
#
# This is a more easily queriable object type to hold SeqRecord objects
# from Entrez.
#
# The JSON-ified version of Genbank is called Polynucleotide, which is
# the schema name used in Clotho.

from Bio import Alphabet, SeqIO
import ClothoSeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord, _RestrictedDict
from Bio.SeqFeature import SeqFeature, Reference, FeatureLocation
import re

class Genbank:
    def __init__(self, seq_record=None):
        if seq_record is None:
            self.record = None
            self.description = "<unknown description>"
            self.sequence = MutableSeq("")
            self.name = "<unknown name>"
            self.id = "<unknown id>"
            self.annotations = {'references': []}
            self.features = []
            self.dbxrefs = []
            self.letter_annotations = _RestrictedDict(length=0)
            self.firstLine = ""
        elif isinstance(seq_record, (str, unicode)):
            self.readFile(seq_record)
        elif isinstance(seq_record, SeqRecord):
            self.record = seq_record
            self.description = seq_record.description
            self.sequence = seq_record.seq
            self.name = seq_record.name
            if seq_record.id is None or seq_record.id == "":
                self.id = seq_record.name
            else:
                self.id = seq_record.id
            self.annotations = seq_record.annotations
            self.features = seq_record.features
            self.dbxrefs = seq_record.dbxrefs
            self.letter_annotations = seq_record.letter_annotations
            self.firstLine = ""
        else:
            raise TypeError

# Most of the methods below this point aren't entirely necessary, and
# in truth, I've taken to just using the dot notation to get and set
# the attributes.
#
# The two methods that are particular useful are:
# readFile() - see line 235
# and
# writeRecord() - see line 219

    def setDescription(self, desc = ""):
        self.description = desc

    def setSeq(self, seq = "", alphabet = Alphabet.generic_alphabet):
        s = Seq(seq, alphabet)
        self.sequence = s

    def setName(self, name = ""):
        self.name = name
        self.annotations['accessions'] = [name]

    def setID(self, id = ""):
        self.id = id

    "***These are all annotation related from here.****"

    def setComment(self, comm = ""):
        if isinstance(comm, str):
            self.annotations['comment'] = comm 
        else:
            print "Please enter in a string."

    def setVersion(self,ver = 0):
        if isinstance(ver, int):
            self.annotations['sequence_version'] = ver
        else:
            print "Please enter in an int."

    def setSource(self, sour = ""):
        if isinstance(sour, str):
            self.annotations['source'] = sour
        else:
            print "Pleae enter in a string."

    def setTaxonomy(self, tax = []):
        if hasattr(tax, '__iter__') and not hasattr(tax, 'strip'):
            self.annotations['taxonomy'] = tax
        else:
            print "Please enter a list."

    def setKeywords(self, key=[]):
        if hasattr(tax, '__iter__') and not hasattr(tax, 'strip'):
            self.annotations['keywords'] = key
        else:
            print "Please enter a list."

    def setDataFileDivision(self, div=""):
        if isinstance(div, str):
            self.annotations['data_file_division'] = div
        else:
            print "Please enter in a string."

    def setDate(self, date=""):
        if isinstance(date, str):
            p = re.compile("\d+-\w+-\d+")
            m = p.match(date)
            if m is None:
                pass
                # print "Please enter a string in the format of: DAY-MONTH-YEAR."
            else:
                self.annotations['date'] = date
        else:
            print "Please enter a string in the format of: DAY-MONTH-YEAR."
            print "Example: '24-FEB-2014'"

    def setOrganism(self, org=""):
        if isinstance(org, str):
            self.annotations['organism'] = org
        else:
            print "Please enter in a string."

    def setGI(self, gi):
        if isinstance(gi, str):
            self.annotations['gi'] = gi
        else:
            print "Please enter in a string."


    "***Reference related methods.***"

    def addReference(self, title="", authors="", journal="", start=0, end=0, strand=None, comment="", consrtm="", medline_id="", pubmed_id=""):
        if title=="" and author=="":
            print "Please enter in at least an author and/or a title."
        else:
            ref = Reference()
            ref.title = title
            ref.authors = authors
            ref.journal = journal
            ref.comment = comment
            ref.consrtm = consrtm
            ref.medline_id = medline_id
            ref.pubmed_id = pubmed_id
            if end != 0:
                ref.location.append(FeatureLocation(start, end, strand))
            self.annotations['references'].append(ref)

    def removeReference(self, title, authors):
        for ref in self.annotations['references']:
            if ref.title == title and ref.authors == authors:
                self.annotations['references'].remove(ref)

    def addLocationToReference(self, title, authors, start, end, strand=None):
        for ref in self.annotations['references']:
           if ref.title == title and ref.authors == authors:
                ref.location.append(FeatureLocation(start, end, strand)) 

    def removeLocationFromReference(self, title, authors, start, end, strand=None):
        for ref in self.annotations['references']:
            if ref.title == title and ref.authors == authors:
                for loc in ref.location:
                    if loc.start == start and loc.end == end and loc.strand == strand:
                        ref.location.remove(loc)

    "***Features.***"

#this still needs work because qualifiers is a dictionary, sub_features is an array - want to do better parsing to avoid errors
    def addFeatures(self, start=0, end=0, strand=None, type="", location_operator="", id="<unknown id>", qualifiers="", sub_features=None, ref=None, ref_db=None):
        if start==0 and end==0:
            print "Please enter a valid argument."
        if not isinstance(qualifiers, dict):
            print "Please enter a dictionary of qualifiers."
        else:
            seq_feature = SeqFeature(FeatureLocation(start, end, strand), type, location_operator, strand, id, qualifiers, sub_features, ref, ref_db)
            self.features.append(seq_feature)

    def removeFeature(self, start, end, strand=None, type=""):
        for feat in self.features:
            loc = feat.location
            print loc
            if loc.start == start and loc.end == end and loc.strand == strand \
            and feat.type == type:
                print feat
                self.features.remove(feat)

    """Set dbxrefs."""
    #This is a list of database cross references (a list of strings)
    #But I haven't seen a single example of this occuring
    def setDBxrefs(self, item = None):
        if item is not None:
            self.dbxrefs.append(item)

    """Set letter_annotations."""
    #This is a restricted dictionary of ptyhon sequences (like lists and strings and tuples) whose length matches the sequence
    #Typically used to hold a list of integers representing sequencing quality scores
    #Or a string representing the secondary structure
    #NOTE: Seems to be obsolete, as far as testing shows.
    def setLetterAnnotations(self, ann = None):
        if ann is not None:
            self.letter_annotations = ann

    """Used for if you simply want to return the record for
    other uses."""
    def getRecord(self):
        self.record = SeqRecord(self.sequence, id=self.id, name=self.name, description=self.description, dbxrefs=self.dbxrefs, features=self.features, annotations=self.annotations, letter_annotations=self.letter_annotations)
        return self.record

    """Used for writing the GenBank record stored here into a
    usable file."""
    def writeRecord(self, file_name):
        self.record = SeqRecord(self.sequence, id=self.id, name=self.name, description=self.description, dbxrefs=self.dbxrefs, features=self.features, annotations=self.annotations, letter_annotations=self.letter_annotations)
        SeqIO.write(self.record, file_name, "gb")
        if self.firstLine == "":
            lines = open(file_name, 'r').readlines()
            self.firstLine = lines[0]
            #file = open(file_name, 'w')
            #for line in lines:
            #    file.write(line)
            #file.close()

    BP_AA_RC = 41
    LIN_CIR = 54

    """Good for either when you already have a GenBank object
    and want to read off from a file, or if you're initializing
    a new GenBank object and pass in the file name as an
    argument."""
    def readFile(self, file_name):
        with open(file_name, 'r') as f:
            lines = f.readlines()
            self.firstLine = lines[0]
            seqStart = False
            featStart = False
            seqString = ""
            INDEX = 1
            sub = 0
            for i in range(len(lines)):
                build = lines[i]
                if (not seqStart):
                    if lines[i].startswith("ORIGIN"):
                        seqStart = True
                        featStart = False
                else:
                    import re
                    if len(build.strip()) > 0:
                        if re.search("\d+", build.split()[0]) is None and not build.startswith("//"):
                            indexStr = str(INDEX)
                            temp = lines[i]
                            build = ""
                            for i in range(9 - len(indexStr)):
                                build += " "
                            build += indexStr + " "
                            for group in range(6):
                                for ten in range(10):
                                    if len(temp) == 0:
                                        break
                                    build += temp[0]
                                    temp = temp[1:]
                                if group != 5:
                                    build += " "
                            INDEX += 60
                            build += "\n"
                            build = build.lower()
                        elif len(build.rstrip()) > 75 and build[0] == " ":
                            if sub == 0:
                                sub = len(build.rstrip())-75
                        build = build[sub:]
                if (not featStart):
                    if lines[i].startswith("FEATURES"):
                        featStart = True
                else:
                    if len(build) > 0 and build[0] == " ":
                        spl = build.split(' /')
                        if len(spl) == 3:
                            build = spl[0] + ' /' + spl[1] + "\n"  + \
                            "                     " + spl[2] + "\n" 
                            #potentially want to add spl[2] back in as + spl[2] + "\n"
                    else:
                        featStart = False
                seqString += build
        with open(file_name, 'w') as f:
            f.write(seqString)
        handle = open(file_name, "rU")
        hold = list(ClothoSeqIO.parse(handle, "gb"))
        hold = hold[0]
        self.record = hold
        #everything into variables to make it queryable
        self.description = hold.description
        self.sequence = hold.seq
        self.name = hold.name
        self.id = hold.id
        self.annotations = hold.annotations
        self.features = hold.features
        self.dbxrefs = hold.dbxrefs
        self.letter_annotations = hold.letter_annotations


    """Allows you to slice a section of the sequence."""
    def seqSlice(self, start=0, length=""):
        if length == "":
            return self.sequence[start:]
        return self.sequence[start:start+length]

    
