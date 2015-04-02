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
# This is used to query for particular lists of accession IDs from Entrez.
# Requires a user email address to pass.

from Bio import Entrez, SeqIO
from protein_holder import Polypeptide

class call_accn:
    def __init__(self, search_database, return_type, email):
        Entrez.email = email
        self.database = search_database
        self.return_type = return_type
        self.email = email
        self.records = []

    """Grab records associated with the accession numbers and add them to
    self.records as GenBank objects.
    Prints out which numbers are successful and which ones are unsuccessful."""
    def retrieve_gb(self, accnList):
        success = []
        failed = []
        for acc in accnList:
            try:
                result_handle = Entrez.efetch(db=self.database, rettype=self.return_type, id=acc)
                for seq_record in SeqIO.parse(result_handle, self.return_type):
                    self.records.append(Polypeptide(seq_record))
                    #print seq_record
                result_handle.close()
                success.append(acc)
            except:
                #print("Failed.")
                failed.append(acc)
        #print "Finished. %d successful retrievals and %d unsuccessful retrievals." % (len(success), len(failed))
        #print "Successful: ", success
        #print "Failed: ", failed

    """Empty the records."""
    def clearRecords(self):
        self.records = []

    """For the string representation of this object."""
    def __str__(self):
        count = 0
        s = ""
        for r in self.records:
            s += "%d\t%s\n" % (count, r.id)
            count += 1
        return s