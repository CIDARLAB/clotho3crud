from Bio import Entrez, SeqIO
from ClothoPy.protein_retrieval import call_accn

"""
Expects @terms in the form of [organism, gene_name, retmax].
retmax should be None if the user wants all.
"""
def _protein_by_gene(terms):
    Entrez.email = 'nobody@example.com'
    call = call_accn('protein', 'gb', 'nobody@example.com')
    term_string = terms[0] + '[ORGN] AND ' + terms[1] + '[GENE]'
    ids = []
    lastlen = 0
    start = 0

    if terms[2] is None:
        while True:
            handle = Entrez.esearch(db="protein", term=term_string, retstart=start, retmax=50)
            record = Entrez.read(handle)
            handle.close()
            ids = ids + record['IdList']
            if lastlen == len(ids):
                break
            else:
                lastlen = len(ids)
            start = start + 50
    else:
        retmax = 50
        while True:
            if retmax > terms[2]:
                retmax = terms[2]
            handle = Entrez.esearch(db="protein", term=term_string, retstart=start, retmax=retmax)
            record = Entrez.read(handle)
            handle.close()
            ids = ids + record['IdList']
            if start + retmax == terms[2]:
                break
            start = start + 50
            if start + retmax > terms[2]:
                retmax = terms[2] - (start)

    call.retrieve_gb(ids)

    prot_list = "["
    for record in call.records:
        prot_list = "\n" + prot_list + record.toJSON() + ","
    prot_list = prot_list[:len(prot_list)-1] + "\n]"

    return prot_list

def run(term):
    return _protein_by_gene(term) #map(_protein_by_gene, terms)