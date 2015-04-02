# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014
#
# This is a more easily queriable object type to hold SeqRecord objects
# from Entrez. Specifically, querying from the Protein database.

from Bio import Alphabet, SeqIO
import ClothoSeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord, _RestrictedDict
from Bio.SeqFeature import SeqFeature, Reference, FeatureLocation
import re
import json
from StringIO import StringIO

class Polypeptide:
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
            out_handle = StringIO()
            ClothoSeqIO.write(seq_record, out_handle, "gb")
            gb_data = out_handle.getvalue().split('\n')
            self.firstLine = gb_data[0]
        else:
            raise TypeError

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

    d = {}
    json = ""

    highlighted = {'ApEinfo_fwdcolor':'forColor', 'ApEinfo_revcolor':'revColor', \
    'inference':'inference', 'label': 'description', 'note': 'notes'}

    def toJSON(self):
        high = self.feature()
        self.d =  {'schema': 'org.clothocad.model.Polypeptide',
        'description': self.description, \
        'sequence': str(self.sequence), \
        'name': self.name, \
        'id': self.name if self.id is None else self.id, \
        'accession': None if (self.annotations['accessions'] == [] or self.annotations['accessions'] is None) else self.annotations['accessions'], \
        'date': self.annotations['date'], \
        'highlight': high, \
        'isLinear': False if 'circular' in self.firstLine.split() else True, #not sure about this? \
        'isSingleStranded': 'ss-' in self.firstLine }
        self.json = json.dumps(self.d, indent=4)
        return self.json

    
    def feature(self):
        import random
        r = lambda: random.randint(0,255)
        #color1 = '#%02X%02X%02X' % (r(),r(),r())
        #color2 = '#%02X%02X%02X' % (r(),r(),r())
        high = self.features
        highs = []
        for f in high:
            color1 = '#%02X%02X%02X' % (r(),r(),r())
            color2 = '#%02X%02X%02X' % (r(),r(),r())
            loc = f.location
            qual = f.qualifiers
            highs.append( 
            {'start':loc.start.position, 'end': loc.end.position, \
            'forColor': color1, 'revColor': color2, 'inference': None, \
            'notes': None, 'description': None }) #notes is supposed to be the qualifier notes if it exists
            #'strand': loc.strand, 'type':f.type, \
            #'id:': f.id, \
            #'location operator':f.location_operator, 'ref': f.ref,\
            #'ref db': f.ref_db} )
            for q in qual.keys():
                if q in self.highlighted.keys():
                    highs[len(highs) - 1][self.highlighted[q]] = qual[q][0]
        return highs