# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014
#
# This is a more easily queriable object type to hold Blast Record objects
# from Entrez.

from Bio.Blast.Record import Blast
import json

class Alignment:
    def __init__(self):
        self.accession = None
        self.evalue = None
        self.identity = None
        self.positives = None
        self.bits = None
        self.score = None
        self.match_start = None
        self.match_end = None
        self.query_start = None
        self.query_end = None
    def as_dict(self):
        return {'accession': self.accession, \
            'evalue': self.evalue, \
            'identity': self.identity, \
            'positives': self.positives, \
            'bits': self.bits, \
            'score': self.score, \
            'match_start': self.match_start, \
            'match_end': self.match_end, \
            'query_start': self.query_start, \
            'query_end': self.query_end}

class Blast_Record:
    def __init__(self, blast=None, seq=""):
        if blast is None:
            self.id = ""
            self.query = ""
            self.alignments = []
        elif isinstance(blast, Blast):
            spl = blast.query_id.split('|')
            if len(spl) >= 4:
                self.id = "ncbi_blast_" + blast.query_id.split('|')[3]
            else:
                self.id = "ncbi_blast_" + spl[0]
            self.query = seq
            self.alignments = []
            for align in blast.alignments:
                temp = Alignment()
                temp.accession = align.accession
                temp.evalue = align.hsps[0].expect
                temp.identity = align.hsps[0].identities
                temp.positives = align.hsps[0].positives
                temp.bits = align.hsps[0].bits
                temp.score = align.hsps[0].score
                temp.match_start = align.hsps[0].sbjct_start
                temp.match_end = align.hsps[0].sbjct_end
                temp.query_start = align.hsps[0].query_start
                temp.query_end = align.hsps[0].query_end
                self.alignments.append(temp)
    def align_as_dict_list(self):
        dictList = []
        for a in self.alignments:
            temp = a.as_dict()
            dictList.append(temp)
        return dictList
    def _json(self):
        align = self.align_as_dict_list()
        return {'schema': 'org.clothocad.model.BlastRecord', \
            'id': self.id, \
            'query': self.query, \
            'alignments': align}