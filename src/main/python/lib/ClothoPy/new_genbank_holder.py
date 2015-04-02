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
# This class holds a new version of the Genbank file format to be used
# in NucSeq.

from Bio.Alphabet import ProteinAlphabet, DNAAlphabet, RNAAlphabet
from genbank_holder import Genbank

class New_Genbank:
    def __init__ (self, record):
        self.GB = Genbank(record) # A lot of the work is done already in GenBank

        self.highs = []
        
        self.name = self.GB.name
        self.description = self.GB.description

        self.sequence = str(self.GB.sequence)
    
        self.date = None
        if 'date' in self.GB.annotations.keys():
            self.date = self.GB.annotations['date']
        if self.GB.id is None or self.GB.id == "":
            self.id == self.GB.name
        else:
            self.id = self.GB.id
        if 'gi' in self.GB.annotations.keys():
            self.id = self.GB.annotations['gi']
        self.accn = None
        if 'accessions' in self.GB.annotations.keys():
            self.accn = self.GB.annotations['accessions']

        self.highlight = self.GB.features
        
        self.firstLine = self.GB.firstLine
        self.isLinear = True
        if 'circular' in self.firstLine.split():
            self.isLinear = False

        self.isSingleStranded = 'ss-' in self.firstLine

    """Fairly self explanatory. Just extracts what the Pubmed IDs are,
    if there are any."""
    def extract_refs(self):
        for ref in self.references:
            pm = ref.pubmed_id
            if pm != "":
                self.pubmed.append(pm)

    def get_record(self):
        return self.GB.recordrecord

    highlighted = {'ApEinfo_fwdcolor':'forColor', 'ApEinfo_revcolor':'revColor', \
    'inference':'inference', 'label': 'description', 'note': 'notes'}

    def highlighter(self):
        import random
        r = lambda: random.randint(0,255)
        color1 = '#%02X%02X%02X' % (r(),r(),r())
        color2 = '#%02X%02X%02X' % (r(),r(),r())
        high = self.highlight
        count = 1;
        for f in high:
            loc = f.location
            qual = f.qualifiers
            self.highs.append( 
            {'start':loc.start.position, 'end': loc.end.position, \
            'forColor': color1, 'revColor': color2, 'inference': None, \
            'notes': None, 'description': None }) #notes is supposed to be the qualifier notes if it exists
            #'strand': loc.strand, 'type':f.type, \
            #'id:': f.id, \
            #'location operator':f.location_operator, 'ref': f.ref,\
            #'ref db': f.ref_db} )
            for q in qual.keys():
                if q in self.highlighted.keys():
                    self.highs[len(self.highs) - 1][self.highlighted[q]] = qual[q][0]
            count += 1

    def _json(self):
        self.highlighter()
        return {'schema': 'org.clothocad.model.Polynucleotide',
        'description': self.description, \
        #'type': self.gb.type, \
        'sequence': self.sequence, \
        'name': self.name, \
        'id': self.id, \
        'accession': self.accn[0], \
        #'organism': self.gb.organism, \
        'date': self.date, \
        'highlight': self.highs, \
        #'data_file_division': self.gb.data_file_division, \
        #'pubmed': self.gb.pubmed, \
        'isLinear': self.isLinear, \
        'isSingleStranded': self.isSingleStranded }