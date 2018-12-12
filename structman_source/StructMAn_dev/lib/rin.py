import gzip
import os

def getIAmap(interaction_score_file):
    f = gzip.open(interaction_score_file,'r')
    lines = f.readlines()
    f.close()

    IAmap = {}

    for line in lines[1:]:
        #A:39:_:ASP (cnt:mc_sc) A:41:_:LEU	2.220600
        if line.count('\t') == 0:
            print interaction_score_file
            break
        score = float(line.split('\t')[1])
        edge = line.split('\t')[0]
        [res_a,interaction_type,res_b] = edge.split()
        [chain_a,res_nr_a,insertioncode_a,res_name_a] = res_a.split(':')
        [chain_b,res_nr_b,insertioncode_b,res_name_b] = res_b.split(':')
        [interaction_base_type,interaction_sub_type] = interaction_type[1:-1].split(':')
        res_nr_a = "%s%s" % (res_nr_a,insertioncode_a.replace('_',''))
        res_nr_b = "%s%s" % (res_nr_b,insertioncode_b.replace('_',''))

        if not chain_a in IAmap:
            IAmap[chain_a] = {}
        if not res_nr_a in IAmap[chain_a]:
            IAmap[chain_a][res_nr_a] = {}
        if not interaction_base_type in IAmap[chain_a][res_nr_a]:
            IAmap[chain_a][res_nr_a][interaction_base_type] = {}
        if not interaction_sub_type in IAmap[chain_a][res_nr_a][interaction_base_type]:
            IAmap[chain_a][res_nr_a][interaction_base_type][interaction_sub_type] = {}
        IAmap[chain_a][res_nr_a][interaction_base_type][interaction_sub_type][(chain_b,res_nr_b)] = score

        if not chain_b in IAmap:
            IAmap[chain_b] = {}
        if not res_nr_b in IAmap[chain_b]:
            IAmap[chain_b][res_nr_b] = {}
        if not interaction_base_type in IAmap[chain_b][res_nr_b]:
            IAmap[chain_b][res_nr_b][interaction_base_type] = {}
        if not interaction_sub_type in IAmap[chain_b][res_nr_b][interaction_base_type]:
            IAmap[chain_b][res_nr_b][interaction_base_type][interaction_sub_type] = {}
        IAmap[chain_b][res_nr_b][interaction_base_type][interaction_sub_type][(chain_a,res_nr_a)] = score

    return IAmap

def getProfile(interaction_map,residue,ligands,res_contig_map):
    profile = {'ligand':[0,0.0],'interchain':[0,0.0],'neighbor':[0,0.0],'short':[0,0.0],'long':[0,0.0]}
    (chain,res) = residue

    if not chain in interaction_map:
        return profile
    if not res in interaction_map[chain]:
        return profile
    if not chain in res_contig_map:
        return profile
    if not res in res_contig_map[chain]:
        return profile
    for (chain_b,res_b) in interaction_map[chain][res]['combi']['all_all']:
        score = interaction_map[chain][res]['combi']['all_all'][(chain_b,res_b)]
        if (chain_b,res_b) in ligands:
            profile['ligand'][0] += 1
            profile['ligand'][1] += score
        elif chain != chain_b:
            profile['interchain'][0] += 1
            profile['interchain'][1] += score
        else:
            if not res_b in res_contig_map[chain]:
                continue
            res_dist = abs(res_contig_map[chain][res][0]-res_contig_map[chain][res_b][0])
            if res_dist < 2:
                profile['neighbor'][0] += 1
                profile['neighbor'][1] += score
            elif res_dist < 6:
                profile['short'][0] += 1
                profile['short'][1] += score
            else:
                profile['long'][0] += 1
                profile['long'][1] += score
    return profile

#called by templateFiltering
def lookup(pdb_id,residues,chain,ligands,res_contig_map,base_path):
    pdb_id = pdb_id.replace('_AU','').lower()
    folder_path = "%s/%s/%s" % (base_path,pdb_id[1:-1],pdb_id)

    network_file = "%s/%s.sif.gz" % (folder_path,pdb_id)
    interaction_score_file = "%s/%s_intsc.ea.gz" % (folder_path,pdb_id)
    interaction_count_file = "%s/%s_nrint.ea.gz" % (folder_path,pdb_id)
    residue_file = "%s/%s_res.txt.gz" % (folder_path,pdb_id)

    if not os.path.isfile(interaction_score_file):
        error = "Did not find RIN: %s" % folder_path
        return {},error

    interaction_map = getIAmap(interaction_score_file)

    profile_map = {}
    if residues != None:
        for residue in residues:
            profile = getProfile(interaction_map,residue,ligands,res_contig_map)
            if profile == None:
                print pdb_id,residue
                continue
            profile_map[residue] = profile
        return profile_map,None
    else:
        for res_nr in res_contig_map[chain]:
            profile = getProfile(interaction_map,(chain,res_nr),ligands,res_contig_map)
            if profile == None:
                print pdb_id,(chain,res_nr)
                continue
            profile_map[(chain,res_nr)] = profile
        return profile_map,None


