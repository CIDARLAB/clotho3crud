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
from fetch_uniprot import fetch_uniprot
import json

"""
Take a straight list of polypeptides.
"""
def process_polypeptide2(poly_list):
	retriever = call_accn('protein', 'gb', 'me@example.com')
	orfs = []
	# pairs = {}
	for poly in poly_list:

		# try UniProt
		next_check = False
		try:
			attempt = fetch_uniprot(poly)
			uni = uniprot_to_ncbi(poly)
			if poly in uni.keys():
				ncbi = uni[poly]
				retriever.retrieve_gb([ncbi])
				record = retriever.records[0]
				orf = _protein_to_orf(record)
				# if orf is None:
					# print "uni none"
				orfs.append(orf)
			else:
				# print "check me"
				next_check = True
		except Exception:
			# print "check me"
			next_check = True

		# NCBI ids
		if next_check:
			retriever.retrieve_gb([poly]) # get from ncbi
			# print len(retriever.records)
			# print "try ncbi"
			if len(retriever.records) > 0:
				record = retriever.records[0]
				orf = _protein_to_orf(record)
				if orf is None:
					# print "ncbi none " + record.id
					# explicit
					next_check = True
				else:
					orfs.append(orf)
					next_check = False
			else:
				# explicit
				next_check = True

		if next_check:
			# print "try blast"
			seq = fetch_uniprot(poly)['sequence']
			blaster = call_blast(seq, 'me@example.com')
			record = blaster.retrieve_gb()
			if record is not None:
				orfs.append(record)
			else:
				raise Exception("No corresponding orf for %s" % poly)

		retriever.clearRecords()
	return orfs

"""
Main function.
"""
def ids_to_operon(query):
	orf_list = process_polypeptide2(query)
	# print orf_list
	rbs_list = [{'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR15', 'name': u'UTR15', 'sequence': u'AATTGTGAGCGGATAACAATTTCACACAGGAAACAGCTATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR12', 'name': u'UTR12', 'sequence': u'aattgtgagcggataacaatttcacacaggaaacaCATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR16', 'name': u'UTR16', 'sequence': u'AATTGTGAGCGGATAACAATTTTAAGGAGATATACATATG'}, {'isIntergenic': True, 'id': u'org.clothocad.model.RBS.interRBS1', 'name': u'interRBS1', 'sequence': u'taaGAGCTCAAGGAAAGAAAAatg'}, {'isIntergenic': True, 'id': u'org.clothocad.model.RBS.interRBS2', 'name': u'interRBS2', 'sequence': u'tagATAACGAGGGCAAAAAatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR1', 'name': u'UTR1', 'sequence': u'ACCCGTTTTTTTGGGCTAACAGGAGGAATTAACCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR10', 'name': u'UTR10', 'sequence': u'ACCCGTTTTTTTGGGCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR11', 'name': u'UTR11', 'sequence': u'GGGCAGATCTTCGAATGCATCGCGCGCACCGTACGTCTCGAGGAATTCCTGCAGGATATCTGGATCCACGAAGCTTCCCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR13', 'name': u'UTR13', 'sequence': u'GGGAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR14', 'name': u'UTR14', 'sequence': u'GGCTAGCCTCGAGAATTCACGCGTGGTACCTCTAGAGTCGACCCGGGCGGCCGCAAAGTTCCCTTTAGTGAGGGTTAATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR17', 'name': u'UTR17', 'sequence': u'GGGCGAATTGGGTACCGGGCCCAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR18', 'name': u'UTR18', 'sequence': u'ACCCGTTTTTTTGGGCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATACCCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR19', 'name': u'UTR19', 'sequence': u'acccgtttttttgggctagcaccgcctatctcgtgtgagataggcggagatacgaactttaagaaggagatatacacatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR2', 'name': u'UTR2', 'sequence': u'GGGAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR20', 'name': u'UTR20', 'sequence': u'AATTGTGAGCGGATAACAATTTCACACAGGAAACAGCCAGTCCGTTTAGGTGTTTTCACGAGCACTTCACCAACAAGGACCATAGcatatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR21', 'name': u'UTR21', 'sequence': u'aattgtgagcggataacaatttcacaCAAGGAGGAAACAGCTatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR22', 'name': u'UTR22', 'sequence': u'aattgtgagcggataacaatttcacacaggaaacagaccatggaattcAGGAGCGACTACatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR23', 'name': u'UTR23', 'sequence': u'aattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgaattcGAGGAAGTGGTATatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR24', 'name': u'UTR24', 'sequence': u'gggcgaattgggtaccgggccccccctcgaggtcgacggtatcgataagcttgatatcgaattcctgcagcccgggGAGGAGAGAAATTatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR3', 'name': u'UTR3', 'sequence': u'GGGCGAATTCGTTAACTTTAAGAAGGAGATATACCCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR4', 'name': u'UTR4', 'sequence': u'accCGTTTTTTTGGGCTAACATCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR5', 'name': u'UTR5', 'sequence': u'accCGTTTTTTTGGGCTATTAATTAAGCGGCCGCCCTGCAGGACTCGAGTTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR6', 'name': u'UTR6', 'sequence': u'ACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR7', 'name': u'UTR7', 'sequence': u'acccgtttttttgggctagcaccgcctatctcgtgtgagataggcggagatacgaactttaagaaggagatatacccatg'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR8', 'name': u'UTR8', 'sequence': u'ACCCGTTTTTTGGGCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}, {'isIntergenic': False, 'id': u'org.clothocad.model.RBS.UTR9', 'name': u'UTR9', 'sequence': u'ATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG'}]
	operon = design_operon(orf_list, rbs_list, 6)
	return operon

def run(query):
	return ids_to_operon(query)