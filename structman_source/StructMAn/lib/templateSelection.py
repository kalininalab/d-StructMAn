import pdbParser as pdb
import templateFiltering
from operator import itemgetter
import multiprocessing

"""
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
    intermediate_results = {}
    for (pdb_id,chain) in structures:
        if pdb_id in intermediate_results:
            resolution,homomer_dict = intermediate_results[pdb_id]
        else:
            resolution,homomer_dict = pdb.getInfo(pdb_id,pdb_path)
            intermediate_results[pdb_id] = resolution,homomer_dict
        if resolution == None:
            continue
        filtered_structures[(pdb_id,chain)] = structures[(pdb_id,chain)]
        filtered_structures[(pdb_id,chain)]['Resolution'] = resolution
        filtered_structures[(pdb_id,chain)]['Oligo'] = homomer_dict[chain]

    return filtered_structures

def paraGetInfo(lock,input_queue,out_queue,pdb_path):
    with lock:
        input_queue.put(None)
    while True:
        with lock:
            inp = input_queue.get()
        if inp == None:
            return
        pdb_id = inp

        resolution,homomer_dict = pdb.getInfo(pdb_id,pdb_path)

        with lock:
            out_queue.put((pdb_id,resolution,homomer_dict))

def filterRawStructureMap(raw_structure_map,pdb_ids,pdb_path,option_res_thresh,n_processes):

    manager = multiprocessing.Manager()
    lock = manager.Lock()

    input_queue = manager.Queue()
    out_queue = manager.Queue()

    info_map = {}
    
    for pdb_id in pdb_ids:
        input_queue.put(pdb_id)

    processes = {}
    for i in range(1,n_processes + 1):

        p = multiprocessing.Process(target=paraGetInfo, args=(lock,input_queue,out_queue,pdb_path))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()

    out_queue.put(None)
    while True:
        out = out_queue.get()
        if out == None:
            break
        (pdb_id,resolution,homomer_dict) = out
        info_map[pdb_id] = (resolution,homomer_dict)

    filtered_structures = {}
    for gene in raw_structure_map:
        filtered_structures[gene] = {}
        for pdb_id,chain in raw_structure_map[gene]:
            if not pdb_id in info_map:
                continue
            resolution,homomer_dict = info_map[pdb_id]
            if resolution == None:
                continue
            if resolution > option_res_thresh:
                continue
            filtered_structures[gene][(pdb_id,chain)] = raw_structure_map[gene][(pdb_id,chain)]
            filtered_structures[gene][(pdb_id,chain)]['Resolution'] = resolution
            filtered_structures[gene][(pdb_id,chain)]['Oligo'] = homomer_dict[chain]

    return filtered_structures

