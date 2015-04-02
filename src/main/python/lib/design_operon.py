# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014

from find_local_hairpins import overall_scoring
from ClothoPy.genbank_holder import Genbank
from StringIO import StringIO
import ClothoPy.ClothoSeqIO
from Bio import Alphabet
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord

def poly_to_gen(poly):
    gen = Genbank()
    
    gen.description = ""
    gen.setName(poly["name"])
    gen.setID(poly["id"])
    date = ""
    gen.setDate(date)
    for high in poly["highlights"]:
        h_start = None
        h_end = None
        h_strand = None
        h_type = "misc_feature"
        h_loc = ""
        h_id = "<unknown id>"
        h_for = ""
        h_rev = ""
        h_inf = ""
        h_desc = ""
        h_note = []
        h_ref = ""

        h_start = high["start"]
        h_end = high["end"]
        h_strand = high["plusStrand"]
        h_for = high["forColor"]
        h_rev = high["revColor"]
        h_desc = high["description"]
        #h_ref = high["refSeq"]

        gen.addFeatures(h_start, h_end, h_strand, "misc_feature", "", h_id, \
            {'ApEinfo_fwdcolor':h_for, 'ApEinfo_revColor':h_rev, \
            'inference':h_inf, 'label':h_desc, 'note':h_note}, \
            None, h_ref)
    
    gen.sequence = poly["sequence"]
    
    typ = "DNA"

    h_single = poly["isSingleStranded"]
    h_linear = poly["isLinear"]
    gen.firstLine = "LOCUS       " + gen.name + "          " + str(len(poly["sequence"])) \
        + " bp " + ("ss-" if h_single else "  ") + typ + "    " + \
        ("linear" if h_linear else "circular") + "      " + date

    for feat in gen.features:
        for key in feat.qualifiers.keys():
            if feat.qualifiers[key] == '' or feat.qualifiers[key] == None or feat.qualifiers[key] == []:
                del feat.qualifiers[key]

    gen.record = SeqRecord(Seq(gen.sequence, Alphabet.IUPAC.IUPACAmbiguousDNA()), gen.id, gen.name, \
        gen.description, None, \
        gen.features, gen.annotations)

    out_handle = StringIO()
    ClothoPy.ClothoSeqIO.write(gen.record, out_handle, "gb")
    gb_data = out_handle.getvalue()

    # clotho.say(gb_data)

    return gb_data #con.d

def concatenate(str1, str2):
    # maxn = 0
    # for n in range(1, 1 + min(len(str1), len(str2))):
    #     suffix = str1[-n:]
    #     prefix = str2[:n]
    #     if prefix == suffix:
    #         maxn = n
    # if maxn >= 3:
    length = len(str1)
    front = str1[:length-3]
    back = str2[3:]
    overlap = str2[:3]
    concat = front + overlap + back
    return concat

def design_operon(orf_list, rbs_list, limit):
    operon = ""
    operon_name = ""
    highlights = []
    location = 0
    first = True

    for orf in orf_list:

        #print orf + " :: before" + str(len(orf))
        orf_seq = orf.sequence.lower()
        if len(orf_seq) >= 38:
            orf_seq = orf_seq[:37]
        #print orf + " :: after" + str(len(orf))
        min_score = float("infinity")
        min_rbs = ""
        is_intergenic_count = 0
        for r in rbs_list:
            #print first
            is_intergenic_count = 0
            is_not_intergenic_count = 0
            if first:
                #print r
                if not r['isIntergenic']:
                    rbs = r['sequence']
                    temp_score = overall_scoring(rbs, orf_seq, limit, False)
                    if temp_score < min_score:
                        min_score = temp_score
                        min_rbs = r
                    is_not_intergenic_count += 1
            else:
                if r['isIntergenic']:
                    rbs = r['sequence']
                    temp_score = overall_scoring(rbs, orf_seq, limit, False)
                    if temp_score < min_score:
                        min_score = temp_score
                        min_rbs = r
                    is_intergenic_count += 1
        
        min_rbs_seq = min_rbs['sequence'].lower()
        orf.sequence = orf.sequence.lower()
        if first:
            first = False
            addition = min_rbs_seq[:len(min_rbs_seq)-3] + 'tag' + orf.sequence[3:len(orf.sequence)-3] + 'taa'
            # print "rbs:\t" + min_rbs_seq[:len(min_rbs_seq)-3] + 'tag'
            # print "orf:\t" + 'tag' + orf.sequence[3:len(orf.sequence)-3] + 'taa'
            # print addition
            # addition = concatenate(min_rbs_seq[:len(min_rbs_seq)-3] + 'tag', 'tag' + orf.sequence[3:len(orf.sequence)-3] + 'taa')
            rbs_start = location
            rbs_end = location + len(min_rbs_seq)
            orf_start = rbs_end - 3
            orf_end = rbs_end - 3 + len(orf.sequence)
        else:
            addition = 'taa' + min_rbs_seq[3:len(min_rbs_seq)-3] + 'tag' + orf.sequence[3:len(orf.sequence)-3] + 'taa'
            # print "rbs:\t" + 'taa' + min_rbs_seq[3:len(min_rbs_seq)-3] + 'tag'
            # print "orf:\t" + 'tag' + orf.sequence[3:len(orf.sequence)-3] + 'taa'
            # print addition
            rbs_start = location - 3
            rbs_end = location - 3 + len(min_rbs_seq)
            orf_start = rbs_end - 3
            orf_end = rbs_end - 3 + len(orf.sequence)
        if is_intergenic_count > 1:
            rbs_list.remove(min_rbs) #remove entry from list of RBSs
        

        # rbs_extra = 0
        # orf_extra = 0
        # if min_rbs_seq.endswith("atg") and orf_seq.startswith("atg"):
        #     # min_rbs_seq = min_rbs_seq[:len(min_rbs_seq) - 3]
        #     orf.sequence = orf.sequence.lower()[3:]
        #     orf_extra = -3
        # else:
        #     rbs_extra = 3
        # addition, str1, str2, overlap1 = concatenate(min_rbs_seq, orf.sequence)
        # print "concat: " + addition + "\nfront: " + str1 + "\nback: " + str2 + "\noverlap: " + overlap1
        # operon, str3, str4, overlap2 = concatenate(operon, addition)

        operon = concatenate(operon, addition)

        operon_name += min_rbs['name'] + "." + orf.name + "."

        # print "rbs:\t" + str(rbs_start) + "\t" + str(rbs_end)
        # print "orf:\t" + str(orf_start) + "\t" + str(orf_end)
        location = len(operon)
        # print "location:\t" + str(location)

        import random
        r = lambda: random.randint(0,255)
        # add RBS
        color1 = '#%02X%02X%02X' % (r(),r(),r())
        color2 = '#%02X%02X%02X' % (r(),r(),r())
        
        highlights.append( {'start': rbs_start, \
            'end': rbs_end, \
            #'refSeq': min_rbs['id'], \
            #'sequence': addition, \
            'schema': 'org.clothocad.model.Highlight', \
            'plusStrand': True, \
            'description': min_rbs['name'], \
            'forColor': color1, \
            'revColor': color2} )
        # add ORF
        color1 = '#%02X%02X%02X' % (r(),r(),r())
        color2 = '#%02X%02X%02X' % (r(),r(),r())

        highlights.append( {'start': orf_start, \
            'end': orf_end, \
            #'refSeq': min_rbs['id'], \
            #'sequence': addition, \
            'schema': 'org.clothocad.model.Highlight', \
            'plusStrand': True, \
            'description': orf.name, \
            'forColor': color1, \
            'revColor': color2} )

        # location += len(addition)
        # print operon

        # name should be a concatenation of rbs1.cds1.rbs2.cds2.etc
        # (orf_list is going to be a list of polynucleotides and not just a straight list)
        # isSingleStranded = False
        # isLinear = True
        # no accession number or submissionDate
        # description?
    operon_name = operon_name[:len(operon_name)-1]
    # print { \
    #     'id': ".", \
    #     'sequence': operon, \
    #     'schema': 'org.clothocad.model.Polynucleotide', \
    #     'highlights': highlights, \
    #     'isLinear': True, \
    #     'isSingleStranded': False, \
    #     'name': operon_name }
    

    return poly_to_gen ({ \
        'id': ".", \
        'sequence': str(operon), \
        'schema': 'org.clothocad.model.Polynucleotide', \
        'highlights': highlights, \
        'isLinear': True, \
        'isSingleStranded': False, \
        'name': operon_name })