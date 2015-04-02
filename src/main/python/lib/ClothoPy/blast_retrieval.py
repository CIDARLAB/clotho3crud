from TBlastN import TBlastN
from blast_holder import Blast_Record
from new_genbank_holder import New_Genbank
from Bio import Entrez, SeqIO

class call_blast:
    def __init__(self, sequence, email):
        Entrez.email = email
        self.database = 'nucleotide'
        self.return_type = 'gb'
        self.email = email
        self.records = []
        self.sequence = sequence

    def retrieve_gb(self):
        uni = TBlastN([self.sequence, None])
        #print dir(uni)
        # return uni
        #print uni.alignments
        seq_start = uni.alignments[0].query_start
        seq_end = uni.alignments[0].query_end
        if uni.alignments[0].identity / float(uni.alignments[0].positives) > 0.9:
            # print "start:\t" + str(seq_start)
            # print "end:\t" + str(seq_end)
            try:
                ncbi = uni.alignments[0].accession # get ncbi number
                result_handle = Entrez.efetch(db=self.database, rettype=self.return_type, id=ncbi, seq_start=seq_start, seq_stop=seq_end)
                for seq_record in SeqIO.parse(result_handle, self.return_type):
                    self.records.append(New_Genbank(seq_record))
                record = None
                if len(self.records) != 0:
                    record = self.records[0]
                self.records = []
                return record # now we have the ncbi record

            except Exception:
                raise Exception("No corresponding orf for %s" % self.sequence)
        else:
            raise Exception("No corresponding orf for %s" % self.sequence)