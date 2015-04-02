from ClothoPy.accn_retrieval import call_accn
from ClothoPy.new_genbank_holder import New_Genbank
from ClothoPy.new_gb_to_json import New_GB_Converter

def _convert_genbank_to_nucseq(genbank_obj):
    return genbank_obj._json()

def run(accession_id):
    retriever = call_accn('nucleotide', 'gb', 'nobody@example.com')
    retriever.retrieve_gb(accession_id)
    return _convert_genbank_to_nucseq(retriever.records[0]) #map(_convert_genbank_to_nucseq, retriever.records)
