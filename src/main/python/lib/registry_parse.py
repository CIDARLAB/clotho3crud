# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

import xml.etree.ElementTree as ET
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_rna, generic_protein
from ClothoPy.genbank_holder import Genbank

"""This function assumes that the XML is formatted as is found on iGem."""
def _reg_to_gen(filename):

	gen = Genbank()

	tree = ET.parse(filename)
	root = tree.getroot()

	part = root[0][0]

	for child in part:
		if child.tag == 'part_id':
			gen.setID(child.text)
		elif child.tag == 'part_name':
			gen.setName(child.text)
		elif child.tag == 'part_short_desc':
			gen.setDescription(child.text)
		elif child.tag == 'part_entered':
			MONTHS = {'01':'JAN', '02':'FEB', '03':'MAR', '04':'APR', '05':'MAY', '06':'JUN', \
			'07':'JUL', '08':'AUG', '09':'SEP', '10':'OCT', '11':'NOV', '12':'DEC'}
			dateSplit = child.text.split("-")
			date = dateSplit[2] + '-' + MONTHS[dateSplit[1]] + '-' + dateSplit[0]
			gen.setDate(date)
		elif child.tag == "features":
			for feat in child:
				fID = feat.findall('id')[0].text
				fTitle = feat.findall('title')[0].text
				fType = feat.findall('type')[0].text
				fDir = feat.findall('direction')[0].text
				if fDir == "forward":
					fDirection = 1
				elif fDir == "reverse":
					fDirection = -1
				else:
					fDirection = None
				fStart = feat.findall('startpos')[0].text
				fEnd = feat.findall('endpos')[0].text
				gen.addFeatures(int(fStart), int(fEnd), fDirection, "misc_feature", "", fID, {}, None, fTitle)
		elif child.tag == "sequences":
			sequence = "".join(child[0].text.split())
			typ = 'DNA'
			alpha = generic_dna
			if 'u' in sequence or 'U' in sequence:
				typ = 'RNA'
				alpha = generic_rna
			if 'm' in sequence or 'M' in sequence or 'r' in sequence or 'R' in sequence:
				typ = ''
				alpha = generic_protein
		
			gen.sequence = MutableSeq(sequence, alpha)
			
	"""
	for feat in gen.features:
			for key in feat.qualifiers.keys():
				if feat.qualifiers[key] == '' or feat.qualifiers[key] == None or feat.qualifiers[key] == []:
					del feat.qualifiers[key]
	"""

	gen.record = SeqRecord(gen.sequence, gen.id, gen.name, gen.description, None, \
	        gen.features, gen.annotations)

	gen.writeRecord('temp.gb')
	g = open('temp.gb', 'rU').read()
	os.remove('temp.gb')
	return g

def run(file_name):
	return _reg_to_gen(file_name) #map(_reg_to_gen, files)