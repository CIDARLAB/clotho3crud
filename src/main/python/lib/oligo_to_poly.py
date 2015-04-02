from StringIO import StringIO
import ClothoPy.ClothoSeqIO

def _oligo_to_poly(oligo):
	fin = {}
	fin['description'] = oligo.description
	fin['sequence'] = oligo.sequence
	fin['name'] = oligo.name
	if hasattr(oligo, "id"):
		fin['id'] = oligo.id
	else:
		fin['id'] = oligo.name
	fin['isLinear'] = True
	fin['isSingleStranded'] = True
	fin['submissionDate'] = None
	fin['accession'] = None
	fin['highlights'] = []
	fin['schema'] = 'org.clothocad.model.Polynucleotide'
	return fin


def run(poly):
    return _oligo_to_poly(poly) #map(_oligo_to_poly, oligo)
