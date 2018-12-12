import pdbParser as pdb
import templateFiltering
from operator import itemgetter

def score(entry):
    #print(entry)
    score = 1.0
    try:
        resolution = float(entry['Resolution'])
    except:
        #print("Resolution buggy in template:")
        #print(entry)
        resolution = 5.0
        entry['Resolution'] = resolution
    try:
        seq_id = float(entry['Seq_Id'])/100
    except:
        #print("Sequence_id buggy in template:")
        #print(entry)
        seq_id = 0.35
        entry['Seq_Id'] = seq_id
    try:
        rel_aln_length = float(entry['Coverage'])
    except:
        #print("Aln_length buggy in template:")
        #print(entry)
        rel_aln_length = 0.35
        entry['Coverage'] = rel_aln_length
    
    score = templateFiltering.qualityScore(resolution,rel_aln_length,seq_id)
    entry.append(score)
    #Expand Template by placeholders for ligand distance and anotation
    entry.append(-0.1)
    entry.append(0.0)
    return entry

"""
def sortByScore(entries):
    entries = sorted(entries, key=itemgetter(7),reverse = True)
    return entries
"""

def selectTemplates(structures,pdb_path):
    if len(structures) == 0:
        return {}
    #print entries
    filtered_structures = {}
    for (pdb_id,chain) in structures:
        resolution = pdb.getInfo(pdb_id,chain,pdb_path)
        if resolution == None:
            continue
        filtered_structures[(pdb_id,chain)] = structures[(pdb_id,chain)]
        filtered_structures[(pdb_id,chain)]['Resolution'] = resolution

    return filtered_structures
    

