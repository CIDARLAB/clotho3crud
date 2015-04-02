# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

import urllib,urllib2

def uniprot_to_ncbi(query):
	url = 'http://www.uniprot.org/mapping/'
	# use REFSEQ_NT_ID for nucleotide, P_REFSEQ_AC is for protein
	params = {'from':'ACC', 'to':'P_REFSEQ_AC','format':'tab','query':query}

	data = urllib.urlencode(params)
	request = urllib2.Request(url, data)
	contact = "me@example.com" # Please set your email address here to help us debug in case of problems.
	request.add_header('User-Agent', 'Python %s' % contact)
	try:
		response = urllib2.urlopen(request)
		page = response.read(200000)
		out = parser(page)
		#print "success"
		return out
	except Exception:
		#print "The query %s failed." % query.strip()
		return {}

def parser(page):
	split_on_break = page.split("\n")
	fin = {}
	for line in split_on_break:
		if line.startswith("From"):
			pass
		elif line == "":
			pass
		else:
			split_on_tab = line.split("\t")
			fin[split_on_tab[0]] = split_on_tab[1]
	return fin