# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014
#
# A set of functions used to convert from a SinglePathway into a list
# of ORFs as GenBank objects. Convert to "Nucleotide" by initializing
# NewGenBank objects with the .record of the GenBank objects.

from ClothoPy.protein_holder import Polypeptide
from pathways.uniprot_to_ncbi import uniprot_to_ncbi
from ClothoPy.accn_retrieval import call_accn
from ClothoPy.BlastP import BlastP
from ClothoPy.blast_holder import Blast_Record
from protein_to_orf import _protein_to_orf

"""
Converts a Pathway object into a list of ORFs.
"""
def pathway_to_orf(pathway):
	enzyme_list = pathway.reactions[0][2].enzymes
	orf_list = []
	#for reaction in pathway.reactions:
	#	enzyme_list.append(chooseEnzyme(reaction))
	for enzyme in enzyme_list:
		orf = process_polypeptide(enzyme)
		if orf != [None]:
			orf_list += orf
	return orf_list

"""
Trivial selector to choose Enzyme.
"""
def choose_enzyme(reaction):
	# for now, this is a trivial selector
	return reaction[2]['enzymes'][0]

"""
For each Polypeptide in the Enzyme object, obtain its ORF.
"""
def process_polypeptide(enzyme):
	retriever = call_accn('nucleotide', 'gb', 'me@example.com')
	orfs = []
	for polypeptide in enzyme.polypeptides:
		for poly in polypeptide:
			# need to put in a check here.
			uni = uniprot_to_ncbi(poly.uniprot)
			if uni != {}:
				ncbi = uniprot_to_ncbi(poly.uniprot)[poly.uniprot] # get ncbi number
				retriever.retrieve_gb([ncbi]) # get from ncbi
				record = retriever.records[0] # now we have the ncbi record

				record.GB.writeRecord(record.id + ".gb")

				orf = _protein_to_orf(record)
				orfs.append(orf)
				retriever.clearRecords()
			else:
				uni = Blast_Record(TBlastN([poly.sequence, None]))
				try:
					ncbi = uni.alignments[0].accession # get ncbi number
					retriever.retrieve_gb([ncbi]) # get from ncbi
					record = retriever.records[0] # now we have the ncbi record

					record.GB.writeRecord(record.id + ".gb")
					
					orfs.append(record)
				except Exception:
					print "No corresponding orf for %s" % poly.uniprot

	return orfs

def run(pathway):
    return pathway_to_orf(pathway) #map(pathway_to_orf, pathways)