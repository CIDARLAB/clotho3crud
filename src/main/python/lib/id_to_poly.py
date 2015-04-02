from registry_collector import _grab_registry
import json
import random

def _id_to_poly(number):
    registry = _grab_registry(number)
    return changer(registry)

def changer(registry):
    change = json.loads(registry)
    fin = {}
    fin['isLinear'] = True
    fin['isSingleStranded'] = False
    fin['schema'] = 'org.clothocad.model.Polynucleotide'
    for key in change.keys():
        if key == 'part_short_desc':
            fin['description'] = change[key]
        if key == 'part_name':
            fin['accession'] = change[key]
            fin['name'] = change[key]
        if key == 'part_entered':
            MONTHS = {'01':'JAN', '02':'FEB', '03':'MAR', '04':'APR', '05':'MAY', '06':'JUN', \
            '07':'JUL', '08':'AUG', '09':'SEP', '10':'OCT', '11':'NOV', '12':'DEC'}
            dateSplit = change[key].split("-")
            date = dateSplit[2] + '-' + MONTHS[dateSplit[1]] + '-' + dateSplit[0]
            fin['submissionDate'] = date
        if key == 'sequences':
            try:
                fin['sequence'] = "".join(change[key][0]['seq_data'].split())
            except Exception:
                pass
        if key == 'id':
            fin['id'] = change[key]
        if key == 'features':
            fin['highlights'] = []
            if change[key] is not None:
                for f in change[key][0]['feature']:
                    #print f['title'], f['type']
                    r = lambda: random.randint(0,255)
                    color1 = '#%02X%02X%02X' % (r(),r(),r())
                    color2 = '#%02X%02X%02X' % (r(),r(),r())

                    temp = {}
                    temp['forColor'] = color1
                    temp['recColor'] = color2
                    temp['inference'] = None
                    temp['description'] = None
                    temp['notes'] = []
                    if f['type'] != "BioBrick" or \
                        (f['type'] == "BioBrick" and not f['title'].startswith('BBa_') \
                        or f['title'] == change['part_name']):

                        temp['refSeq'] = None
                        for q in f.keys():
                            if q == 'startpos':
                                temp['start'] = f[q]
                            if q == 'endpos':
                                temp['end'] = f[q]
                            if q == 'title':
                                if f[q] == None or f[q] == "":
                                    temp['name'] = "unnamed"
                                else:
                                    temp['name'] = f[q]
                            if q == 'direction':
                                if f[q] == 'forward':
                                    temp['plusStrand'] = True
                                else:
                                    temp['plusStrand'] = False
                        fin['highlights'].append(temp)
                    else:
                        temp['start'] = f['startpos']
                        temp['end'] = f['endpos']
                        if f['title'] == None or f['title'] == "":
                            temp['name'] = "unnamed"
                        else:
                            temp['name'] = f['title']
                        if f['direction'] == 'forward':
                            temp['plusStrand'] = True
                        else:
                            temp['plusStrand'] = False
                        jj = _grabRegistry(f['title'])
                        res = changer(jj)
                        load = json.loads(res)
                        concat = dict({'refSeq': load}.items() + temp.items())
                        fin['highlights'].append(concat)
    return json.dumps(fin, indent=4)

def run(ids):
    return _id_to_poly(ids) #map(_id_to_poly, ids)
