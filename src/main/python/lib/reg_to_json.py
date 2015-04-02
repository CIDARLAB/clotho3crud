# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

import xml.etree.ElementTree as ET

"""
Converts from registry (using the filename) to JSON.
"""
def reg_to_json(filename):
	tree = ET.parse(filename)
	root = tree.getroot()
	return root

def rec_traverse(root):
	if len(root.getchildren()) == 0:
		return root.text
	else:
		tree = {}
		for child in root:
			if len(child.getchildren()) == 0:
				tree[child.tag] = rec_traverse(child)
			elif child.tag not in tree.keys():
				tree[child.tag] = [rec_traverse(child)]
			else:
				tree[child.tag].append(rec_traverse(child))
		return tree

def reg_parse(inName, outName):
	root = reg_to_json(inName)
	fin = {root.tag: [rec_traverse(root)]}
	import json
	j = json.dumps(fin, indent=4)
	with open(outName, 'w') as f:
		f.write(j)

# reg_parse('registry/BBa_A340620.xml', 'BBa_A340620.json')
# reg_parse('registry/BBa_B0000.xml', 'BBa_B0000.json')
# reg_parse('registry/BBa_B0017.xml', 'BBa_B0017.json')
# reg_parse('registry/BBa_B0104.xml', 'BBa_B0104.json')