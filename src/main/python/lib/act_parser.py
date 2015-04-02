# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

import json
from fetch_uniprot import fetch_uniprot

chem_dict = {}

"""
This is the main function, which takes in a dictionary of the result from act.20.
"""
def act_parser(j):
	global chem_dict
	final_pathways = []
	pathways = j['results'][0]['pathways']['pathways']
	chemicals = j['results'][0]['pathways']['chemicals']

	description = j['results'][0]['pathways']['description']
	target = j['results'][0]['pathways']['target']
	version = j['results'][0]['pathways']['version']
	name = target + " Single Pathway"
	
	for chemical in chemicals:

		chem_dict[str(int(chemical['id']))] = ('InChI', chemical['InChI'], chemical['name'])
	for path in pathways:
		single_pathway = {'description': description, 'target': target, 'version': version, 'name': name}
		interms = [] # final list of intermediates
		reacts = [] # final list of reactants
		count = 0
		reactions = path['rxns_steps']
		for reaction in reactions:
			intermediates = reaction['substrates']
			products = reaction['products']
			reaction_desc = reaction['readable']

			# add this reaction's substrates to interms
			interm_temp, count = make_intermediates(interms, intermediates, count)
			interms = interms + interm_temp

			# add this reaction's products to interms
			react_temp, count = make_intermediates(interms, products, count)
			interms = interms + react_temp

			#make reactions to add to reacts
			r = make_reaction(interm_temp, react_temp, reaction['sequences'], reaction_desc)

			for i in range(len(interm_temp)):
				for j in range(len(react_temp)):
					reacts.append( [ interm_temp[i]['index'], react_temp[j]['index'], r ] )
		for i in interms:
			i.pop("index", None)
		single_pathway['intermediates'] = interms
		single_pathway['reactions'] = reacts
		single_pathway['schema'] = "org.clothocad.model.SinglePathway"
		final_pathways.append(single_pathway)
	return final_pathways

"""
Create intermediates.
"""
def make_intermediates(interms, intermediates, count):
	global chem_dict

	interm_temp = []

	for chem in intermediates['chemicals']:
		interm = {"index": count}
		if chem in chem_dict.keys():
			if not in_intermediates(interms, chem_dict[chem][1]) and not in_intermediates(interm_temp, chem_dict[chem][1]):
				interm[chem_dict[chem][0]] = chem_dict[chem][1]		# InChI or SMILES : 
				interm['name'] = chem_dict[chem][2]					# name :
				interm_temp.append(interm)
				count += 1
		else:
			if not in_intermediates2(interms, chem) and not in_intermediates(interm_temp, chem):
				interm['name'] = chem
				interm_temp.append(interm)
				count += 1

	# there isn't always a 'cofactors' field
	if 'cofactors' in intermediates.keys():
		for chem in intermediates['cofactors']:
			interm = {"index": count}
			if chem in chem_dict.keys():
				if not in_intermediates(interms, chem_dict[chem][1]) and not in_intermediates(interm_temp, chem_dict[chem][1]):
					interm[chem_dict[chem][0]] = chem_dict[chem][1]		# InChI or SMILES : 
					interm['name'] = chem_dict[chem][2]					# name :
					interm_temp.append(interm)
					count += 1
			else:
				if not in_intermediates2(interms, chem) and not in_intermediates(interm_temp, chem):
					interm['name'] = chem 								# name :
					interm_temp.append(interm)
					count += 1
	return (interm_temp, count)

"""
Is this in intermediates?
Specifically for InChI and SMILES check.
"""
def in_intermediates(intermediates, check):
	for i in intermediates:
		if 'InChI' in i.keys():
			if i['InChI'] == check:
				return True
		# if 'SMILES' in i.keys():
		# 	if i['SMILES'] == check:
		#		return True
	return False

"""
Is this in intermediates?
Specifically for name check, aka no InChI or SMILES.
"""
def in_intermediates2(intermediates, name):
	for i in intermediates:
		if i['name'] == name:
			return True
	return False

"""
Make a reaction.
"""
def make_reaction(reactants, products, enzymes, description):
	fin = {"schema": "org.clothocad.model.Reaction", "description": description}
	reacts = [] #create the reactants
	for react in reactants:
		if 'InChI' in react.keys():
			reacts.append( "InChI=" + react['InChI'] )
		else:
			reacts.append( react['name'] )
	fin['reactants'] = reacts
	prods = [] #create the products
	for prod in products:
		if 'InChI' in prod.keys():
			prods.append( "InChI=" + prod['InChI'] )
		else:
			reacts.append( react['name'] )
	fin['products'] = prods
	fin['enzymes'] = make_enzyme(enzymes)
	return fin

"""
Make the enzyme list for reaction.
"""
def make_enzyme(enzymes):
	fin = []
	for enzyme in enzymes:
		temp = {"schema": "org.clothocad.model.Enzyme"}
		temp['polypeptides'] = make_polypeptide(enzyme['seq'])
		temp['observations'] = []
		fin.append(temp)
	return fin

"""
Make a polypeptide.
"""
def make_polypeptide(sequence): # should probably have some more inputs, like uniprot and such
	fin = []
	for seq in sequence:
		uniprot = fetch_uniprot(seq)
		fin.append( {'sequence': uniprot['sequence'], 'uniprot': seq, \
			'source': uniprot['source'], 'schema': 'org.clothocad.model.Polypeptide', \
			'name': uniprot['name']} )
	return fin