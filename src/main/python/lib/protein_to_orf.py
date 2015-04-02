from Bio import Entrez, SeqIO
from ClothoPy.new_genbank_holder import New_Genbank
from ClothoPy.protein_retrieval import call_accn
from Bio.Seq import Seq

"""
Expecting a Polypeptide object.
"""
def _protein_to_orf(prot):

    if hasattr(prot, "record"):
        taken = prot.record
    else:
        retriever = call_accn('protein', 'gb', 'me@example.com')
        retriever.retrieve_gb([prot.id])
        taken = retriever.records[0].record
    coded_by = []
    trigger = False

    for feat in taken.features:
        if feat.type == 'CDS' and 'coded_by' in feat.qualifiers.keys():
            coded_by = feat.qualifiers['coded_by']
            trigger = True
        else:
            pass

    if not trigger:
        return None
    
    transform = ""
    spl = coded_by[0].split('.')
    if '(' in coded_by[0]:
        transform = spl[0].split('(')[0]
        key = spl[0].split('(')[1]
        seqStart = spl[1].split(':')[1]
        seqEnd = spl[3][:len(spl[3])-1]
    else:
        key = spl[0]
        seqStart = spl[1].split(':')[1]
        seqEnd = spl[3]

    records = []

    try:
        result_handle = Entrez.efetch(db='nucleotide', rettype='gb', id=key, seq_start=seqStart, seq_stop=seqEnd)
        for seq_record in SeqIO.parse(result_handle, 'gb'):
            records.append(New_Genbank(seq_record))
        result_handle.close()
    except:
        return None

    # to deal with complements
    if transform == "complement":
        records[0].GB.sequence = records[0].GB.sequence.reverse_complement()
        records[0].sequence = str(records[0].GB.sequence)
        records[0].name = "reverse_complement(" + records[0].name + ")"
        records[0].GB.name = "reverse_complement(" + records[0].name + ")"

    return records[0]

def run(term):
    return _protein_to_orf(term) #map(_protein_to_orf, terms)