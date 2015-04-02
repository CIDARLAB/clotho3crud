# Attribution Information:
# Written by Mina Li, working under Professor Christopher Anderson.
# For any inqueries, please email li.mina888@berkeley.edu.
# (C) 2014
#
# A set of functions used to find all local hairpins.

# How many hydrogen bonds there are per basepair.
hydrogen_bonds = \
    {'aa': 0, \
    'at': 2, \
    'ac': 0, \
    'ag': 0, \
    'ta': 2, \
    'tt': 0, \
    'tc': 0, \
    'tg': 1, \
    'ca': 0, \
    'ct': 0, \
    'cc': 0, \
    'cg': 3, \
    'ga': 0, \
    'gt': 1, \
    'gc': 3, \
    'gg': 0, \
    'a': 0, \
    't': 0, \
    'c': 0,\
    'g': 0 }

# Holds all possible hairpins at the end.
all_hairpins = []

"""
Loops over the sequence to generate the basepairs.
"""
def obtain_basepairs(seq, loop_size):
    pairs = []
    length = len(seq)
    if length == 8 + loop_size: # full size
        stem1 = seq[:4]
        loop = seq[4: 4 + loop_size]
        stem2 = seq[4 + loop_size :]
        for i in range(4):
            pairs.append( stem1[4 - i - 1] + stem2[i] )
    elif length > 4 +loop_size: # missing complete right stem
        stem1 = seq[:4]
        loop = seq[4: 4 + loop_size]
        stem2 = seq[4 + loop_size :]
        for i in range(4):
            if i < len(stem2):
                pairs.append( stem1[4 - i - 1] + stem2[i] )
            else:
                pairs.append( stem1[4 - i - 1] )
    elif length >= 4: # missing right stem all together
                        # might be missing complete loop
        stem1 = seq[:4]
        loop = seq[4:]
        for i in range(4):
            pairs.append( stem1[4 - i - 1] )
    elif length < 4: # missing complete loop
        stem1 = seq
        for i in range(length):
            pairs.append( stem1[length - i - 1] )
    return pairs

"""
Calculates the number of hydrogen bonds using the mapping at the top.
"""
def basepair_total(bps):
    total = 0
    for bp in bps:
        if bp in hydrogen_bonds.keys():
            total += hydrogen_bonds[bp] # lookup bonds number in hydrogen_bonds
    return total

"""
Hashses out the name for each potential hairpin loop.
"""
def name(loop_size, index, total):
    hairpin_str = str(loop_size) + "loop-" + str(index) + "-" + str(total)
    all_hairpins.append( hairpin_str ) # just want to hold all the possible hairpins
    return hairpin_str

"""
Takes in the 5' UTR, CDS 37 fragment, the loop size, and the hydrogen bond
limit to qualify for further consideration.
We return only the hashes of the hairpins that have more hydrogen bonds than
the limit, although all potential hashes are stored in all_hairpins.
"""
def hairpin(utr_5, cds_37, loop_size, limit):
    utr_5 = utr_5.lower()
    cds_37 = cds_37.lower()
    utr_len = len(utr_5)
    if utr_5.endswith("atg"):
        utr_5 = utr_5[:utr_len - 3]
    sequence = utr_5 + cds_37

    max_size = loop_size + 8
    length = len(sequence)

    seq = ""
    fin_list = []

    for i in range(length):

        if i + max_size < length:
            seq = sequence[i:i + max_size]
        else:
            seq = sequence[i:]

        bps = obtain_basepairs(seq, loop_size)
        total = basepair_total(bps)

        hairpin_str = name(loop_size, utr_len - i, total)
        if total > limit:
            fin_list.append(hairpin_str)
    #print fin_list
    return fin_list

"""
Takes in the 5' UTR, CDS 37 fragment, and the hairpins generated from hairpin().
We return the total hydrogen bonds present on loops that overlap that region from
the start codon to 12 basepairs before it on the 5' UTR.
"""
def hairpin_scoring(utr_5, cds_37, hairpins):
    utr_5 = utr_5.lower()
    cds_37 = cds_37.lower()
    if hairpins != []:
        loop_size = int(hairpins[0][:1])
        loop_length = 8 + loop_size
        fin_total = 0
        for hairpin in hairpins:
            split = hairpin.split("-")
            if "" in split:
                index = int(split[2])
                bonds = int(split[3])
            else:
                index = int(split[1])
                bonds = int(split[2])
            total = index + loop_length
            if (total <= 12 and total > 0) or (index <= 12 and index > 0):
                fin_total += bonds
        return fin_total
    return 0

"""
Scoring hairpin scores from loop_size 3 through 9.
"""
def overall_scoring(utr_5, cds_37, limit, val_bool):
    total = 0
    for loop_size in range(3, 10):
        pins = hairpin(utr_5, cds_37, loop_size, limit)
        val = hairpin_scoring(utr_5, cds_37, pins)
        #comment out if you don't want to see individual results
        #print "LOOP SIZE: " + str(loop_size) + "\t" + "SCORE: " + str(val)
        total += val
    return total

"""
Scoring entryway.
"""
def _scorer(arr):
    return overall_scoring(arr[0], arr[1], arr[2], False)

def run(arr):
    return _scorer(arr) #map(_scorer, arrs)