# File to Genbank

from ClothoPy.genbank_holder import Genbank
from StringIO import StringIO
import ClothoPy.ClothoSeqIO

def _convert_file(file_name):
	genbank = Genbank(file_name)
	out_handle = StringIO()
	ClothoPy.ClothoSeqIO.write(genbank.record, out_handle, "gb")
	gb_data = out_handle.getvalue()
	return gb_data 

def run(accession_id):
    return _convert_file(accession_id) #map(_convertFile, accession_ids)