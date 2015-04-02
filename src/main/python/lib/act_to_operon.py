from act_parser import act_parser
from act_query import act_query
from ClothoPy.protein_holder import Polypeptide
from uniprot_to_ncbi import uniprot_to_ncbi
from ClothoPy.protein_retrieval import call_accn
from TBlastN import TBlastN
from ClothoPy.blast_holder import Blast_Record
from protein_to_orf import _protein_to_orf
from design_operon import design_operon
from grab_rbs import grab_rbs
from ClothoPy.blast_retrieval import call_blast
from ClothoPy.new_genbank_holder import New_Genbank
import json

"""
Converts a Pathway object into a list of ORFs.
"""
def pathway_to_orf(pathway):
	enzyme_list = pathway['reactions'][0][2]['enzymes']
	orf_list = []
	for enzyme in enzyme_list:
		orf = process_polypeptide(enzyme)
		if orf != [None]:
			orf_list += orf
	return orf_list

"""
For each Polypeptide in the Enzyme object, obtain its ORF.
"""
def process_polypeptide(enzyme):
	retriever = call_accn('protein', 'gb', 'me@example.com')
	orfs = []
	# pairs = {}
	for poly in enzyme['polypeptides']:
		# print "UniProt:\t" + str(poly['uniprot'])
		uni = uniprot_to_ncbi(str(poly['uniprot']))
		if uni != {}:
			# print poly['uniprot']
			ncbi = uni[str(poly['uniprot'])] # get ncbi number
			# print "Attempt UniProt to NCBI protein:\t" + ncbi

			retriever.retrieve_gb([ncbi]) # get from ncbi
			record = retriever.records[0] # now we have the ncbi record
			orf = _protein_to_orf(record)
			if orf is not None:
				orfs.append(orf)
				# print "NCBI:\t" + orf.name
				# orf.GB.writeRecord(orf.id + ".gb")
				# pairs[str(poly['uniprot'])] = orf.id
			else:
				blaster = call_blast(poly['sequence'], 'me@example.com')
				record = blaster.retrieve_gb()
				if record is not None:
					# record.GB.writeRecord(record.id + ".gb")
					orfs.append(record)
					# print "Blast:\t" + record.name
					# pairs[str(poly['uniprot'])] = record.id
				else:
					print 'failed'
			retriever.clearRecords()
		else:
			blaster = call_blast(poly['sequence'], 'me@example.com')
			record = blaster.retrieve_gb()
			if record is not None:
				# record.GB.writeRecord(record.id + ".gb")
				orfs.append(record)
				# print "Blast2:\t" + record.name
				# pairs[str(poly['uniprot'])] = record.id
		# print pairs
	return orfs

demo_results = {'X1JXA6': '306490827', 'W8UAZ9': '374362062', 'H2FX92': '374333820', 'B4QPF7': '195591271', 'R4VV39': '507519518', 'B4LE71': '195375794', 'K0I1J3': '407936729', 'R9VNN5': '512647881', 'Q291N5': '198456965'}
"""
Added for demo purposes.
"""
def locate_polypeptide(pathway):
	enzyme_list = pathway['reactions'][0][2]['enzymes']
	orf_list = []
	for enzyme in enzyme_list:
		orfs = []
		for poly in enzyme['polypeptides']:
			if str(poly['uniprot']) in demo_results.keys():
				orf_id = demo_results[str(poly['uniprot'])]
				orf_file = "/Users/mina/Documents/clothopy/boston_demo/" + orf_id + ".gb"
				gen = New_Genbank(orf_file)
				orfs.append(gen)
		if orfs != [None]:
			orf_list += orfs
	return orf_list

"""
Trivial selector.
"""
def select_pathway(pathways):
	return pathways[0] #this is a trivial selection method, which should be changed

RBSs = [{'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR15', 'name': u'UTR15', 'sequence': u'AATTGTGAGCGGATAACAATTTCACACAGGAAACAGCTATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR12', 'name': u'UTR12', 'sequence': u'aattgtgagcggataacaatttcacacaggaaacaCATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR16', 'name': u'UTR16', 'sequence': u'AATTGTGAGCGGATAACAATTTTAAGGAGATATACATATG'}, {'isIntergenic': True, 'id': u'org.clothocad.model.RBS.interRBS1', 'name': u'interRBS1', 'sequence': u'taaGAGCTCAAGGAAAGAAAAatg'}, {'isIntergenic': True, 'id': u'org.clothocad.model.RBS.interRBS2', 'name': u'interRBS2', 'sequence': u'tagATAACGAGGGCAAAAAatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR1', 'name': u'UTR1', 'sequence': u'ACCCGTTTTTTTGGGCTAACAGGAGGAATTAACCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR10', 'name': u'UTR10', 'sequence': u'ACCCGTTTTTTTGGGCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR11', 'name': u'UTR11', 'sequence': u'GGGCAGATCTTCGAATGCATCGCGCGCACCGTACGTCTCGAGGAATTCCTGCAGGATATCTGGATCCACGAAGCTTCCCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR13', 'name': u'UTR13', 'sequence': u'GGGAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR14', 'name': u'UTR14', 'sequence': u'GGCTAGCCTCGAGAATTCACGCGTGGTACCTCTAGAGTCGACCCGGGCGGCCGCAAAGTTCCCTTTAGTGAGGGTTAATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR17', 'name': u'UTR17', 'sequence': u'GGGCGAATTGGGTACCGGGCCCAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR18', 'name': u'UTR18', 'sequence': u'ACCCGTTTTTTTGGGCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATACCCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR19', 'name': u'UTR19', 'sequence': u'acccgtttttttgggctagcaccgcctatctcgtgtgagataggcggagatacgaactttaagaaggagatatacacatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR2', 'name': u'UTR2', 'sequence': u'GGGAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR20', 'name': u'UTR20', 'sequence': u'AATTGTGAGCGGATAACAATTTCACACAGGAAACAGCCAGTCCGTTTAGGTGTTTTCACGAGCACTTCACCAACAAGGACCATAGcatatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR21', 'name': u'UTR21', 'sequence': u'aattgtgagcggataacaatttcacaCAAGGAGGAAACAGCTatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR22', 'name': u'UTR22', 'sequence': u'aattgtgagcggataacaatttcacacaggaaacagaccatggaattcAGGAGCGACTACatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR23', 'name': u'UTR23', 'sequence': u'aattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgaattcGAGGAAGTGGTATatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR24', 'name': u'UTR24', 'sequence': u'gggcgaattgggtaccgggccccccctcgaggtcgacggtatcgataagcttgatatcgaattcctgcagcccgggGAGGAGAGAAATTatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR3', 'name': u'UTR3', 'sequence': u'GGGCGAATTCGTTAACTTTAAGAAGGAGATATACCCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR4', 'name': u'UTR4', 'sequence': u'accCGTTTTTTTGGGCTAACATCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR5', 'name': u'UTR5', 'sequence': u'accCGTTTTTTTGGGCTATTAATTAAGCGGCCGCCCTGCAGGACTCGAGTTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR6', 'name': u'UTR6', 'sequence': u'ACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR7', 'name': u'UTR7', 'sequence': u'acccgtttttttgggctagcaccgcctatctcgtgtgagataggcggagatacgaactttaagaaggagatatacccatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR8', 'name': u'UTR8', 'sequence': u'ACCCGTTTTTTGGGCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR9', 'name': u'UTR9', 'sequence': u'ATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}]

"""
Main function.
"""
def act_to_operon(query):
	# try:
	# 	#print "regular query"
	# 	j = act_query(query)  #replace with opening 1-Butanol_example.json
	# except Exception:
		#print "abnormal query"
	f = open("/Users/mina/Documents/clothopy/boston_demo/1-Butanol.json", "r")
	j = json.loads(f.read())
	f.close()
	# try:
	# 	paths = act_parser(j)
	# 	selected = select_pathway(paths)
	# except:
	f = open("/Users/mina/Documents/clothopy/boston_demo/paths.json", "r")
	selected = json.loads(f.read())
	f.close()
	# try:
	# 	print "regular query"
	# 	orf_list = pathway_to_orf(selected)
	# except Exception as e:
	print "abnormal query"
	# print e
	orf_list = locate_polypeptide(selected[0])
	# for orf in orf_list:
	# 	print "ORF:\t" + orf.name + "\t" + orf.sequence + "\n"
	try:
		print "regular query"
		rbs_list = grab_rbs()
		if rbs_list == []:
			print "was empty"
			rbs_list = RBSs
	except Exception:
		print "abnormal query"
	rbs_list = [{'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR15', 'name': u'UTR15', 'sequence': u'AATTGTGAGCGGATAACAATTTCACACAGGAAACAGCTATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR12', 'name': u'UTR12', 'sequence': u'aattgtgagcggataacaatttcacacaggaaacaCATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR16', 'name': u'UTR16', 'sequence': u'AATTGTGAGCGGATAACAATTTTAAGGAGATATACATATG'}, {'isIntergenic': True, 'id': u'org.clothocad.model.RBS.interRBS1', 'name': u'interRBS1', 'sequence': u'taaGAGCTCAAGGAAAGAAAAatg'}, {'isIntergenic': True, 'id': u'org.clothocad.model.RBS.interRBS2', 'name': u'interRBS2', 'sequence': u'tagATAACGAGGGCAAAAAatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR1', 'name': u'UTR1', 'sequence': u'ACCCGTTTTTTTGGGCTAACAGGAGGAATTAACCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR10', 'name': u'UTR10', 'sequence': u'ACCCGTTTTTTTGGGCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR11', 'name': u'UTR11', 'sequence': u'GGGCAGATCTTCGAATGCATCGCGCGCACCGTACGTCTCGAGGAATTCCTGCAGGATATCTGGATCCACGAAGCTTCCCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR13', 'name': u'UTR13', 'sequence': u'GGGAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR14', 'name': u'UTR14', 'sequence': u'GGCTAGCCTCGAGAATTCACGCGTGGTACCTCTAGAGTCGACCCGGGCGGCCGCAAAGTTCCCTTTAGTGAGGGTTAATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR17', 'name': u'UTR17', 'sequence': u'GGGCGAATTGGGTACCGGGCCCAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR18', 'name': u'UTR18', 'sequence': u'ACCCGTTTTTTTGGGCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATACCCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR19', 'name': u'UTR19', 'sequence': u'acccgtttttttgggctagcaccgcctatctcgtgtgagataggcggagatacgaactttaagaaggagatatacacatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR2', 'name': u'UTR2', 'sequence': u'GGGAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR20', 'name': u'UTR20', 'sequence': u'AATTGTGAGCGGATAACAATTTCACACAGGAAACAGCCAGTCCGTTTAGGTGTTTTCACGAGCACTTCACCAACAAGGACCATAGcatatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR21', 'name': u'UTR21', 'sequence': u'aattgtgagcggataacaatttcacaCAAGGAGGAAACAGCTatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR22', 'name': u'UTR22', 'sequence': u'aattgtgagcggataacaatttcacacaggaaacagaccatggaattcAGGAGCGACTACatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR23', 'name': u'UTR23', 'sequence': u'aattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgaattcGAGGAAGTGGTATatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR24', 'name': u'UTR24', 'sequence': u'gggcgaattgggtaccgggccccccctcgaggtcgacggtatcgataagcttgatatcgaattcctgcagcccgggGAGGAGAGAAATTatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR3', 'name': u'UTR3', 'sequence': u'GGGCGAATTCGTTAACTTTAAGAAGGAGATATACCCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR4', 'name': u'UTR4', 'sequence': u'accCGTTTTTTTGGGCTAACATCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR5', 'name': u'UTR5', 'sequence': u'accCGTTTTTTTGGGCTATTAATTAAGCGGCCGCCCTGCAGGACTCGAGTTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR6', 'name': u'UTR6', 'sequence': u'ACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR7', 'name': u'UTR7', 'sequence': u'acccgtttttttgggctagcaccgcctatctcgtgtgagataggcggagatacgaactttaagaaggagatatacccatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR8', 'name': u'UTR8', 'sequence': u'ACCCGTTTTTTGGGCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR9', 'name': u'UTR9', 'sequence': u'ATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}]
	print rbs_list
	operon = design_operon(orf_list, rbs_list, 6)
	return operon

def run(query):
	return act_to_operon(query)