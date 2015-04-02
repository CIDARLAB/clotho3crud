# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

import urllib
import xml.etree.ElementTree as ET

"""
Obtain the root of the JSON tree.
"""
def reg_root(data):
	root = ET.fromstring(data)
	return root

"""
Traverse from the root of the JSON tree.
"""
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

"""
Parse the string representation of JSON to a Python dictionary.
"""
def reg_parse(inString):
	root = reg_root(inString)
	fin = rec_traverse(root)['part_list'][0]['part'][0]
	fin["schema"] = "org.registry.model.Part"
	fin['id'] = 'org.registry.part.' + fin['part_name']
	import json
	j = json.dumps(fin, indent=4)
	return j

"""
Query the iGem server to obtain the registry entry.
"""
def _grab_registry(theID):
	base = "http://parts.igem.org/xml/part."
	url = base + theID
	filename = theID + ".xml"
	socket = urllib.urlopen(url)
	hold = ""
	for line in socket:
		hold += line
	jReg = reg_parse(hold)
	return jReg
	
def run(ids):
	return _grab_registry(ids) #map(_grab_registry, ids)
