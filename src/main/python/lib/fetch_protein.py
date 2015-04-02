from ClothoPy.protein_retrieval import call_accn

def _fetch_protein(polypeptide):
	return polypeptide.toJSON()

def run(accession_id):
	retriever = call_accn('protein', 'gb', 'nobody@example.com')
    retriever.retrieve_gb(accession_ids)
    return _fetch_protein(retriever.records[0]) #map(_fetch_protein, retriever.records)