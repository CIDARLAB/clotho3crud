# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

import urllib
import io
import itertools as IT
import xml.etree.ElementTree as ET

def myfromstring(content):
    try:
        tree = ET.fromstring(content)
    except ET.ParseError as err:
        lineno, column = err.position
        line = next(IT.islice(io.BytesIO(content), lineno))
        caret = '{:=>{}}'.format('^', column)
        err.msg = '{}\n{}\n{}'.format(err, line, caret)
        raise 
    return tree

"""
This function fetches a UniProt record using its ID. Only UniProt and not UniRef, etc.
"""
def fetch_uniprot(ID):
	en = 'http://www.uniprot.org/uniprot/' #base url
	url = en + ID + ".xml"
	out = urllib.urlopen(url)
	s = out.read() #xml file
	out.close()

	#with open("OUT.txt", "w") as f:
	#	f.write(s)

	myfromstring(s)
	root = ET.fromstring(s)
	entry = root.getchildren()[0]

	# these things you just have to go through uniprot xml files to realize
	sequence = "".join(entry.find('{http://uniprot.org/uniprot}sequence').text.split('\n'))
	organism = entry.find('{http://uniprot.org/uniprot}organism').getchildren()[0].text #scientific name

	protein = entry.find('{http://uniprot.org/uniprot}protein')
	recommendedName = None
	submittedName = None
	name = '<no name found>'


	if protein is not None:
		recommendedName = protein.find('{http://uniprot.org/uniprot}recommendedName')

	if recommendedName is not None:
		name = recommendedName.find('{http://uniprot.org/uniprot}fullName').text
	else:
		submittedName = protein.find('{http://uniprot.org/uniprot}submittedName')

	if submittedName is not None and name is None:
		if submittedName.text != "\n":
			name = submittedName.text
			print name

	return {'uniprot': ID, 'name': name, 'source': organism, 'sequence': sequence}

def run(ids):
    return fetch_uniprot(ids) #map(fetch_uniprot, ids)