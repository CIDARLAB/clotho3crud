# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

import urllib
import json

"""
Query from ACT to obtain JSON.
"""
def act_query(chemical):
	chemical = chemical.encode('ascii','ignore')
	en = urllib.quote("".join(str({"chemical":chemical}).strip()))
	#print ('http://act.20n.com:27080/api/sampleFAB/_find?criteria=' + en).replace("%27", "%22")
	u = urllib.urlopen(('http://act.20n.com:27080/api/sampleFAB/_find?criteria=' + en).replace("%27", "%22"))
	s = u.read()
	u.close()
	#print s
	j = json.loads(s)
	return j