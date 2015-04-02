# Genbank to Polynucleotide (New Genbank)

from ClothoPy.genbank_holder import Genbank
from ClothoPy.new_genbank_holder import New_Genbank
from ClothoPy.new_gb_to_json import New_GB_Converter
from Bio.SeqRecord import SeqRecord
import ClothoPy.ClothoSeqIO

def _convertGB(gb):
    from StringIO import StringIO
    gb_handle = StringIO(gb)
    record = ClothoPy.ClothoSeqIO.read(gb_handle, 'gb')
    con = New_GB_Converter(New_Genbank(record))
    con.convert()
    return con.d

def run(accession_id):
    return _convertGB(accession_id) #map(_convertGB, accession_ids)