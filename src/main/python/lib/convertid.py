# ID to Genbank

from ClothoPy.accn_retrieval import call_accn
from Bio import SeqIO
from StringIO import StringIO


def _convertID(accession_id):
	retriever = call_accn('nucleotide', 'gb', 'nobody@example.com')
	retriever.retrieve_gb([accession_id])
	out_handle = StringIO()
	SeqIO.write(retriever.records[0].GB.record, out_handle, "gb")
	gb_data = out_handle.getvalue()
	return gb_data

def run(accession_id):
    return _convertID(accession_id) #map(_convertID, accession_ids)
