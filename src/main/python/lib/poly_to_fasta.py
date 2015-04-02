from StringIO import StringIO
from convertpoly import _convert_polynucleotide_to_genbank
import ClothoPy.ClothoSeqIO

def _poly_to_fasta(poly):
	gen_string = _convert_polynucleotide_to_genbank(poly)
	gen_handle = StringIO(gen_string)
	gen_record = ClothoPy.ClothoSeqIO.read(gen_handle, 'gb')
	out_handle = StringIO()
	ClothoPy.ClothoSeqIO.write(gen_record, out_handle, "fasta")
	fasta_data = out_handle.getvalue()
	return fasta_data


def run(poly):
    return _poly_to_fasta(poly) #map(_poly_to_fasta, polys)
