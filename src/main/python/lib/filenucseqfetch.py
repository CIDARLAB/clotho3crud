from ClothoPy.new_genbank_holder import New_Genbank
from ClothoPy.new_gb_to_json import New_GB_Converter
from ClothoPy.genbank_holder import Genbank
from ClothoPy.accn_retrieval import call_accn

def _convert_file_to_nucseq(file_name):
	genbank = GenBank(file_name)
	con = New_GB_Converter(NewGenBank(gen.record))
    con.convert()
	return con.d

def run(accession_id):
	retriever = call_accn('nucleotide', 'gb', 'nobody@example.com')
    retriever.retrieve_gb(accession_ids)
    return _convert_file_to_nucseq(retriver.records[0]) #map(_convert_genbank_to_nucseq, retriever.records)
