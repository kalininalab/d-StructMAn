import pdbParser as pdb
import templateFiltering
from operator import itemgetter

def addPdbInfo(entry,pdb_path):
    pdb_id = entry[0]
    chain = entry[2]
    [resolution,r_value] = pdb.getInfo(pdb_id,chain,pdb_path)
    if resolution == None:
        return None
    entry.append(resolution)
    #Placeholder for Interaction_partners
    entry.append([])
    entry.append(r_value)
    return entry

def score(entry):
    #print(entry)
    score = 1.0
    try:
        resolution = float(entry[4])
    except:
        #print("Resolution buggy in template:")
        #print(entry)
        resolution = 5.0
        entry[4] = resolution
    try:
        seq_id = float(entry[1])/100
    except:
        #print("Sequence_id buggy in template:")
        #print(entry)
        seq_id = 0.35
        entry[1] = seq_id
    try:
        rel_aln_length = float(entry[3])
    except:
        #print("Aln_length buggy in template:")
        #print(entry)
        rel_aln_length = 0.35
        entry[3] = rel_aln_length
    try:
        r_value = float(entry[6])
    except:
        #print("R_value buggy in template:")
        #print(entry)
        r_value = 0.63
        entry[6] = r_value
    
    score = templateFiltering.qualityScore(resolution,rel_aln_length,seq_id,r_value)
    entry.append(score)
    #Expand Template by placeholders for ligand distance and anotation
    entry.append(-0.1)
    entry.append(0.0)
    #Expand Template by placeholders for Original_chains
    entry.append("")
    entry.append("")
    return entry

def sortByScore(entries):
    entries = sorted(entries, key=itemgetter(7),reverse = True)
    return entries

def selectTemplates(entries,pdb_path):
    if len(entries) == 0:
        return [],[]
    #print entries
    new_entries = []
    del_entries = []
    for entry in entries:
        #print entry
        entry = addPdbInfo(entry,pdb_path)
        if not entry == None:
            entry = score(entry)
            new_entries.append(entry)
        else:
            del_entries.append(entry)
    #entries = sortByScore(entries)
    return new_entries,del_entries
    

