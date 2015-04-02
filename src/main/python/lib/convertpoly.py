# Polynucleotide to Genbank

from ClothoPy.genbank_holder import Genbank
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_rna, generic_protein
from Bio import SeqIO
from StringIO import StringIO

def _convert_polynucleotide_to_genbank(poly):
	gen = Genbank()
	
	if hasattr(poly, "description"):
		gen.description = poly.description
	if hasattr(poly, "name"):
		gen.setName(poly.name)
	if hasattr(poly, "id"):
		gen.setID(poly.id)
	date = ""
	if hasattr(poly, "submissionDate"):
		MONTHS = {'01':'JAN', '02':'FEB', '03':'MAR', '04':'APR', '05':'MAY', '06':'JUN', \
		'07':'JUL', '08':'AUG', '09':'SEP', '10':'OCT', '11':'NOV', '12':'DEC'}
		dateSplit = poly.submissionDate.split("-")
		date = dateSplit[2][:2] + '-' + MONTHS[dateSplit[1]] + '-' + dateSplit[0]
		gen.setDate(date)
	if hasattr(poly, "annotations"):
		if "accessions" in gen.annotations.keys():
			gen.annotations['accessions'] = poly.accession
			gen.setVersion(int(poly.accession[len(poly.accession) - 1:len(poly.accession)]))
	if hasattr(poly, "highlights"):
		for high in poly.highlights:
			h_start = None
			h_end = None
			h_strand = None
			h_type = "misc_feature"
			h_loc = ""
			h_id = "<unknown id>"
			h_for = ""
			h_rev = ""
			h_inf = ""
			h_desc = ""
			h_note = []
			h_ref = ""
			if hasattr(high, "start"):
				h_start = high.start
			if hasattr(high, "end"):
				h_end = high.end
			if hasattr(high, "plusStrand"):
				h_strand = high.plusStrand
			if hasattr(high, "id"):
				h_id = high.id
			if hasattr(high, "forColor"):
				h_for = high.forColor
			if hasattr(high, "revColor"):
				h_rev = high.revColor
			if hasattr(high, "inference"):
				h_inf = high.inference
			if hasattr(high, "description"):
				h_desc = high.description
			if hasattr(high, "notes"):
				h_note = high.notes
			if hasattr(high, "refSeq"):
				h_ref = high.refSeq
			gen.addFeatures(h_start, h_end, h_strand, "misc_feature", "", h_id, \
				{'ApEinfo_fwdcolor':h_for, 'ApEinfo_revColor':h_rev, \
				'inference':h_inf, 'label':h_desc, 'note':h_note}, \
				None, h_refe)
	
	typ = 'DNA'
	alpha = generic_dna
	if 'u' in poly.sequence or 'U' in poly.sequence:
		typ = 'RNA'
		alpha = generic_rna
	if 'm' in poly.sequence or 'M' in poly.sequence or 'r' in poly.sequence or 'R' in poly.sequence:
		typ = ''
		alpha = generic_protein
	
	gen.sequence = MutableSeq(poly.sequence, alpha)
	
	h_single = False
	if hasattr(poly, "isSingleStranded"):
		h_single = poly.isSingleStranded
	h_linear = True
	if hasattr(poly, "isLinear"):
		h_linear = poly.isLinear
	gen.firstLine = "LOCUS       " + gen.name + "          " + str(len(poly.sequence)) \
		+ " bp " + ("ss-" if h_single else "  ") + typ + "    " + \
		("linear" if h_linear else "circular") + "      " + date

	for feat in gen.features:
		for key in feat.qualifiers.keys():
			if feat.qualifiers[key] == '' or feat.qualifiers[key] == None or feat.qualifiers[key] == []:
				del feat.qualifiers[key]

	gen.record = SeqRecord(gen.sequence, gen.id, gen.name, \
        gen.description, None, \
        gen.features, gen.annotations)

	out_handle = StringIO()
	SeqIO.write(gen.record, out_handle, "gb")
	gb_data = out_handle.getvalue()

	return gb_data #con.d

def run(accession_id):
	return _convert_polynucleotide_to_genbank(accession_id) #map(_convert_polynucleotide_to_genbank, accession_ids)