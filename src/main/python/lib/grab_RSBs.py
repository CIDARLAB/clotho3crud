# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

import urllib
import json

def grab_RSBs():
	# do i need to specify if using localhost or if using clotho website?
	u = urllib.urlopen('https://localhost:8443/data/query?schema=org.clothocad.model.RBS')
	s = u.read()
	u.close()
	j = json.loads(s)
	data = j['data']
	rbs_list = []
	for d in data:
		rbs_list.append({'sequence': d['utr5'], 'id': d['id']})
		#saving both the sequence and the name (or should it be in the ID?)
	return rbs_list


"""
This is meant as a tester, since RSBs are not into the Clotho system yet.
"""
def grab_polynucleotides():
	# do i need to specify if using localhost or if using clotho website?
	u = urllib.urlopen('https://localhost:8443/data/query?schema=org.clothocad.model.Polynucleotide')
	s = u.read()
	u.close()
	j = json.loads(s)
	data = j['data']
	rbs_list = []
	for d in data:
		rbs_list.append({'sequence': d['sequence'], 'id': d['id'], 'name': d['name']})
		#saving both the sequence and the name (or should it be in the ID?)
	return rbs_list

def run():
    return fetch_uniprot