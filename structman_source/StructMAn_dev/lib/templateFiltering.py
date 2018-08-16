import os
import sys
import traceback
import subprocess
import math

from Bio.PDB import *
from Bio.SubsMat import MatrixInfo
from operator import itemgetter

import pdbParser as pdb
import sphecol as siss
import rin

residue_max_acc = { 
# Miller max acc: Miller et al. 1987 http://dx.doi.org/10.1016/0022-2836(87)90038-6 
# Wilke: Tien et al. 2013 http://dx.doi.org/10.1371/journal.pone.0080635 
# Sander: Sander & Rost 1994 http://dx.doi.org/10.1002/prot.340200303 
    'Miller': { 
        'ALA': 113.0, 'ARG': 241.0, 'ASN': 158.0, 'ASP': 151.0, 
        'CYS': 140.0, 'GLN': 189.0, 'GLU': 183.0, 'GLY': 85.0, 
        'HIS': 194.0, 'ILE': 182.0, 'LEU': 180.0, 'LYS': 211.0, 
        'MET': 204.0, 'PHE': 218.0, 'PRO': 143.0, 'SER': 122.0, 
        'THR': 146.0, 'TRP': 259.0, 'TYR': 229.0, 'VAL': 160.0 
    }, 
    'Wilke': { 
        'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 
        'CYS': 167.0, 'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0, 
        'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0, 
        'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0, 
        'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0 
    }, 
    'Sander': { 
        'ALA': 106.0, 'ARG': 248.0, 'ASN': 157.0, 'ASP': 163.0, 
        'CYS': 135.0, 'GLN': 198.0, 'GLU': 194.0, 'GLY': 84.0, 
        'HIS': 184.0, 'ILE': 169.0, 'LEU': 164.0, 'LYS': 205.0, 
        'MET': 188.0, 'PHE': 197.0, 'PRO': 136.0, 'SER': 130.0, 
        'THR': 142.0, 'TRP': 227.0, 'TYR': 222.0, 'VAL': 142.0 
    } 
}

#called by templateSelection
def qualityScore(resolution,rel_aln_length,seq_id,r_value,resolution_wf=0.25,rel_aln_length_wf=0.5,seq_id_wf=1.0,r_value_wf=0.1):
    seq_id = float(seq_id)
    resolution = float(resolution)
    r_value = float(r_value)
    rel_aln_length = float(rel_aln_length)
    if seq_id > 1.0:
        seq_id = seq_id/100.0
    if rel_aln_length > 1.0:
        rel_aln_length = rel_aln_length/100.0
    
    #project the criteria to [0,1] via logistic regression
    resolution_value = (1+math.exp((1.5*resolution)-4))**(-1)
    r_value = 1.0 - r_value
    seq_value = (1+math.exp(10*(0.4-seq_id)))**(-1)

    ws = sum((resolution_wf,rel_aln_length_wf,seq_id_wf,r_value_wf))
    #print resolution_value,rel_aln_length,seq_id,r_value
    quality = sum((resolution_value*resolution_wf,rel_aln_length*rel_aln_length_wf,seq_value*seq_id_wf,r_value*r_value_wf))/ws
    return quality

#called by database
def candidateScore(lig_sub_dist,chain_sub_dist,lig_wf=1.0,chain_wf=1.0,useBlosum=False,aac="",blosum_wf=0.5):
    if useBlosum:
        if aac == "":
            raise NameError("If useBlosum is True, an aac is needed")
        #print MatrixInfo.blosum62
        try:
            blosum_value =  0.6 - float(MatrixInfo.blosum62[(aac[0],aac[-1])])/10
        except:
            blosum_value =  0.6 - float(MatrixInfo.blosum62[(aac[-1],aac[0])])/10
        if blosum_value < 0.0:
            blosum_value = 0.0
        #print blosum_value
    else:
        blosum_value = 1.0
        blosum_wf = 0.0
        
    #project the criteria to [0,1] via logistic regression
    lig_value = (1+math.exp(lig_sub_dist-10))**(-1)
    chain_value = (1+math.exp(chain_sub_dist-10))**(-1)

    if lig_sub_dist == -1:
        lig_value = 0.0
    if chain_sub_dist == -1:
        chain_value = 0.0

    ws = sum((lig_wf,chain_wf,blosum_wf))

    candidate = sum((lig_value*lig_wf,chain_value*chain_wf,blosum_value*blosum_wf))/ws
    return candidate 

def calcDSSP(path,DSSP,angles=False):

    dssp_dict = {}
    errorlist = []

    p = subprocess.Popen([DSSP,path], universal_newlines=True, 
    stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    out, err = p.communicate() 

    # Alert user for errors 
    if err.strip():
        if not out.strip():
            errorlist.append(('DSSP failed to produce an output\n%s\n%s\n' % (path,err),'','',''))
            return dssp_dict,errorlist

    lines = out.split('\n')
    
    i = 1
    for line in lines:
        words = line.split()
        if words[0] == '#':
            break
        i += 1


    for line in lines[i:]:
        #print line
        #print line[:-1]
        if len(line) < 38:
            continue
        if line[9] == " ": 
            continue 

        res = (line[5:11]).replace(" ","") #this includes the insertion code
        insertion_code = line[10]
        chain = line[11]
        ssa = line[16]
        acc = float(line[34:38])
        aa_type_one_letter = line[13]
        try:
            aa = Polypeptide.one_to_three(aa_type_one_letter)
            macc = residue_max_acc['Sander'][aa]
            racc = acc/macc
        except:
            racc = None

        if len(line) > 115:
            phi = float(line[103:109].replace(" ",""))
            psi = float(line[109:115].replace(" ",""))

        if not chain in dssp_dict:
            dssp_dict[chain] = {}
        dssp_dict[chain][res] = (racc,ssa)
        if angles:
            dssp_dict[chain][res] = (racc,ssa,phi,psi)

    return dssp_dict,errorlist

def parsePDB(input_page,chain_id_map):
    """
    Parses a PDB-file and takes all atomic coordinates.

    Input:
    input_page: String ; content of a pdb file

    Output:
    coordinate_map: {Chain:[{Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]},{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}
    """
    #siss_map: {String:[String,{String:(String,float,float,float)}]} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.

    lines = input_page.split('\n')

    coordinate_map = {}
    siss_map = {}
    modres_map = {}
    res_contig_map = {}
    contig_help_map = {}
    ligands = set()

    for line in lines:
        if len(line) > 5:
            record_name = line[0:6].replace(" ","")
            if record_name == "ENDMDL":
                #print len(coordinate_map)
                break
        #ignore short lines
        if len(line) > 20:
            atom_nr = line[6:11].replace(" ","")
            if record_name.count('ATOM') > 0 and record_name != 'ATOM': #100k atom bug fix
                atom_nr = '%s%s' % (record_name[4:],atom_nr)
                record_name = 'ATOM'
            atom_name = line[12:16].replace(" ","")
            res_name = line[17:20].replace(" ","")
        
            if len(line) > 21:
                if record_name == 'MODRES':
                    chain_id = line[16]
                    res_nr = line[18:23].replace(" ","")
                    res_name = line[24:27].replace(" ","")
                    if not (chain_id,res_nr) in modres_map:
                        modres_map[(chain_id,res_nr)] = res_name

                chain_id = line[21]
                res_nr = line[22:27].replace(" ","") #[22:27] includes the insertion_code
                #insertion_code = line[26]
                
            if record_name == "ATOM":
                if len(line) > 50 and not atom_name[0] in ('H','D'):
                    x = float(line[30:38].replace(" ",""))
                    y = float(line[38:46].replace(" ",""))
                    z = float(line[46:54].replace(" ",""))
                    if not chain_id in coordinate_map:
                        coordinate_map[chain_id] = [{},{}]
                    if res_nr not in coordinate_map[chain_id][0]:
                        coordinate_map[chain_id][0][res_nr] = [res_name,{}]
                    coordinate_map[chain_id][0][res_nr][1][atom_nr] = (atom_name,x,y,z)

                    if not chain_id in siss_map:
                        siss_map[chain_id] = {}
                    if not res_nr in siss_map[chain_id]:
                        siss_map[chain_id][res_nr] = [res_name,{}]
                    siss_map[chain_id][res_nr][1][atom_nr] = (atom_name,x,y,z)
                    
                    if not chain_id_map[chain_id] in res_contig_map:
                        res_contig_map[chain_id_map[chain_id]] = {res_nr:1}
                        contig_help_map[chain_id_map[chain_id]] = 1
                    elif not res_nr in res_contig_map[chain_id_map[chain_id]]:
                        contig_help_map[chain_id_map[chain_id]] += 1
                        res_contig_map[chain_id_map[chain_id]][res_nr] = contig_help_map[chain_id_map[chain_id]]

            if record_name == "HETATM":
                if len(line) > 50:
                    x = float(line[30:38].replace(" ",""))
                    y = float(line[38:46].replace(" ",""))
                    z = float(line[46:54].replace(" ",""))
                    if not chain_id in coordinate_map:
                        coordinate_map[chain_id] = [{},{}]
                    
                    if (chain_id,res_nr) in modres_map: #If it is modified residue, than add it to the normal residues...
                        if res_nr not in coordinate_map[chain_id][0]:
                            coordinate_map[chain_id][0][res_nr] = [res_name,{}]
                        coordinate_map[chain_id][0][res_nr][1][atom_nr] = (atom_name,x,y,z)

                        if not chain_id_map[chain_id] in res_contig_map:
                            res_contig_map[chain_id_map[chain_id]] = {res_nr:1}
                            contig_help_map[chain_id_map[chain_id]] = 1
                        elif not res_nr in res_contig_map[chain_id_map[chain_id]]:
                            contig_help_map[chain_id_map[chain_id]] += 1
                            res_contig_map[chain_id_map[chain_id]][res_nr] = contig_help_map[chain_id_map[chain_id]]
                    else:    
                        if res_nr not in coordinate_map[chain_id][1]:#If not, then add it to the ligands
                            coordinate_map[chain_id][1][res_nr] = [res_name,{}]
                        coordinate_map[chain_id][1][res_nr][1][atom_nr] = (atom_name,x,y,z)
                        ligands.add((chain_id_map[chain_id],res_nr))
                        
    return coordinate_map,siss_map,res_contig_map,ligands

def dist(coord1,coord2):
    diff = [coord1[0]-coord2[0],coord1[1]-coord2[1],coord1[2]-coord2[2]]
    d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
    return d

def getMinSubDist(c_map,target_res_id,res_chain,chain):
    atomlist = c_map[res_chain][0][target_res_id][1]
    test_coord = atomlist.values()[0][1:]

    top10 = []
    for res in c_map[chain][0]:
        test_coord_2 = c_map[chain][0][res][1].values()[0][1:]
        d = dist(test_coord,test_coord_2)
        if len(top10) < 9:
            top10.append((res,d))
        elif len(top10) == 9:
            top10.append((res,d))
            top10 = sorted(top10,key = lambda x:x[1])
        elif d < top10[9][1]:
            top10[9] = (res,d)
            top10 = sorted(top10,key = lambda x:x[1])

    min_sub_d = None
    for (res,d) in top10:
        for atomnr in c_map[chain][0][res][1]:
            #exclude Hydrogens
            if c_map[chain][0][res][1][atomnr][0][0] != 'H':
                coord1 = c_map[chain][0][res][1][atomnr][1:]
                for atomnr2 in atomlist:
                    if atomlist[atomnr2][0][0] != 'H':
                        coord2 = atomlist[atomnr2][1:]
                        d = dist(coord1,coord2)
                        if min_sub_d == None or d < min_sub_d:
                            min_sub_d = d
                            min_res = res
                            min_atom_sub = atomnr2
                            min_atom_chain = atomnr

    return min_sub_d,min_res,min_atom_sub,min_atom_chain

#coordinate_map: {Chain:[{Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]},{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}
def getMinDist(c_map,res,chain,res2,chain2,het=0,het2=0):
    min_d = None
    if not res in c_map[chain][het]:
        return None,None,None
    if not res2 in c_map[chain2][het2]:
        return None,None,None

    for atomnr in c_map[chain][het][res][1]:
        coord = c_map[chain][het][res][1][atomnr][1:]
        for atomnr2 in c_map[chain2][het2][res2][1]:
            coord2 = c_map[chain2][het2][res2][1][atomnr2][1:]
            d = dist(coord,coord2)
            if min_d == None or d < min_d:
                min_d = d
                min_atom = atomnr
                min_atom2 = atomnr2

    if min_d == None:
        print c_map,res,chain,res2,chain2,het,het2

    return min_d,min_atom,min_atom2

#called by serializedPipeline
def structuralAnalysis(sub_info_map,pdb_id,stored_annotations,pdb_path,dssp_path,rin_db_path,anno_anno=False,dssp=True,calculate_interaction_profiles=True):
    #print pdb_id,sub_info_map

    page,chain_id_map = pdb.standardParsePDB(pdb_id,pdb_path,return_ori_chain_map=True)

    if page == '':
        print "Error while parsing: ",pdb_id
        return {},{}

    coordinate_map,siss_coord_map,res_contig_map,ligands = parsePDB(page,chain_id_map)

    errorlist = []

    if dssp:
        #write temp pdb file
        tmp_path = 'tmp_%s.pdb' % pdb_id
        f = open(tmp_path,'w')
        f.write(page)
        f.close()

        #call DSSP  -structure of dssp_dict: {Chain:{Res_id:(racc,ssa)}}
        dssp_dict,errlist = calcDSSP(tmp_path,dssp_path)
        errorlist += errlist
        #remove tmp_file
        try:
            os.remove(tmp_path)
        except:
            pass
        if dssp_dict == {}:
            dssp = False

    #t_m_amount = 0

    #if not dssp:
    #surface-calculation with SISI
    target_residues_map = {}
    target_residues = []
    for template_id in sub_info_map:

        (sub_infos,target_chain,oligos) = sub_info_map[template_id]
        #t_m_amount += len(sub_infos)
        if not target_chain in coordinate_map:
            continue
        for aac_base in sub_infos:
            sub_info = sub_infos[aac_base]
            substitution_residue_template = sub_info[1]
            if substitution_residue_template != "-":
                target_res_id = sub_info[0]
                substitution_residue_template = sub_info[1]
                if substitution_residue_template == "-" or not target_res_id in coordinate_map[target_chain][0]:
                    continue

                target_residues.append((chain_id_map[target_chain],target_res_id))

                if target_chain not in target_residues_map:
                    target_residues_map[target_chain] = []
                target_residues_map[target_chain].append(target_res_id)
                #if target_res_id not in siss_coord_map:
                #    print template,target_res_id,siss_coord_map

    siss_map = {}
    for target_chain in target_residues_map:
        #siss_coord_map: {Chain:{String:[String,{String:(String,float,float,float)}]}} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.
        try:
            centroid_map = siss.calcCentroidMap(siss_coord_map[target_chain],target_residues_map[target_chain],False)

            dist_matrix = siss.calcDistMatrix(siss_coord_map[target_chain],centroid_map,target_residues_map[target_chain],False)

        
            siss_map[target_chain] = siss.calculateSiss(siss_coord_map[target_chain],centroid_map,dist_matrix,target_residues_map[target_chain],False)
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            errorlist.append(("Siss error: %s,%s\n%s\n" % (pdb_id,target_chain,str(target_residues_map[target_chain])),e,f,g))

    profiles = {}
    error = None
    if calculate_interaction_profiles:
        profiles,error = rin.lookup(pdb_id,target_residues,ligands,res_contig_map,rin_db_path)
    if error != None:
        errorlist.append((error,'','',''))

    annotation_dict = {}

    #anno_amount = 0
    
    for template_id in sub_info_map:
        (sub_infos,target_chain,oligos) = sub_info_map[template_id]
        if not target_chain in coordinate_map:
            continue
        annotation_dict[template_id] = [pdb_id,{}] #to save the pdb_id here is useful later for the database
        for aac_base in sub_infos:
            if template_id in stored_annotations:
                if aac_base in stored_annotations[template_id]:
                    continue #Do not recalculate, whats already in the database
            sub_info = sub_infos[aac_base]
            target_res_id = sub_info[0]

            substitution_residue_template = sub_info[1]
            if substitution_residue_template == "-" or not target_res_id in coordinate_map[target_chain][0]:
                annotation_dict[template_id][1][aac_base] = 0
                continue
            min_chain_dists = {}
            min_dists = {}
            
            

            #Distance Calculations
            for chain in coordinate_map:
                #Sub - Chain - Calculations
                if chain != target_chain and len(coordinate_map[chain][0]) > 0:
                    min_chain_dist,min_res,atom_sub,atom_chain = getMinSubDist(coordinate_map,target_res_id,target_chain,chain)
                    min_chain_dists[chain] = (min_chain_dist,(atom_chain,atom_sub),min_res)
                #Sub - Lig - Calculations
                for hetres in coordinate_map[chain][1]:
                    min_d,atom,atom2 = getMinDist(coordinate_map,target_res_id,target_chain,hetres,chain,het2=1)
                    abr = coordinate_map[chain][1][hetres][0]
                    ligand_identifier = "%s_%s_%s" % (abr,hetres,chain)
                    min_dists[ligand_identifier] = (min_d,(atom,atom2))
            
            #coordinate_map: {Chain:[{Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]},{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}    
            #Homomer Sub - Sub Calculations:
            homomer_map = {}
            residue_type = coordinate_map[target_chain][0][target_res_id][0]
            for homo_chain in oligos:
                if homo_chain == target_chain:
                    continue
                if not homo_chain in coordinate_map:
                    continue
                if not target_res_id in coordinate_map[homo_chain][0]:
                    continue #Sanity Check - residue-id must be the same in the homo-chain
                if coordinate_map[homo_chain][0][target_res_id][0] != residue_type:
                    continue #Sanity Check - residue type of the homomer residue must be the same as the residue type of the target residue
                min_d,atom,atom2 = getMinDist(coordinate_map,target_res_id,target_chain,target_res_id,homo_chain)                
                homomer_map[homo_chain] = (min_d,atom,atom2)        

            
            if dssp:
                if target_chain in dssp_dict:
                    if target_res_id in dssp_dict[target_chain]:
                        (raac,ssa) = dssp_dict[target_chain][target_res_id]
                    else:
                        if target_res_id in siss_map:
                            siss_value = siss_map[target_res_id]
                        else:
                            siss_value = None
                        ssa = None 
                        raac  = siss_value
                else:
                    [e,f,g] = sys.exc_info()
                    g = traceback.format_exc(g)
                    errorlist.append(("dssp error: chain not in dssp_dict; %s; %s" % (pdb_id,target_chain),e,f,g))
                    if target_res_id in siss_map:
                        siss_value = siss_map[target_res_id]
                    else:
                        siss_value = None
                    ssa = None 
                    raac  = siss_value
            else:
                if target_res_id in siss_map:
                    siss_value = siss_map[target_res_id]
                else:
                    siss_value = None
                ssa = None 
                raac  = siss_value
            if (chain_id_map[target_chain],target_res_id) in profiles:
                profile = profiles[(chain_id_map[target_chain],target_res_id)]
            else:
                profile = None
            annotation_dict[template_id][1][aac_base] = (min_dists,min_chain_dists,raac,target_res_id,ssa,homomer_map,profile)
            #anno_amount += 1

    

    anno_anno_dict = {}
    if anno_anno:
        #compute the distances between all mutations for a whole structure
        used_template_ids = set([])
        for template_id in sub_info_map:
            (sub_infos,target_chain,oligos) = sub_info_map[template_id]
            if not target_chain in coordinate_map:
                continue
            for template_id_2 in sub_info_map:
                if template_id_2 in used_template_ids:
                    continue
                
                (sub_infos_2,target_chain_2,oligos_2) = sub_info_map[template_id_2]
                if not target_chain_2 in coordinate_map:
                    continue
                for aac_base in sub_infos:

                    sub_info = sub_infos[aac_base]
                    res = sub_info[0]
                    residue_template = sub_info[1]
                    if residue_template == '-' or not res in coordinate_map[target_chain][0]:
                        continue

                    for aac_base_2 in sub_infos_2:
                        if template_id in stored_annotations and template_id_2 in stored_annotations:
                            if aac_base in stored_annotations[template_id] and aac_base_2 in stored_annotations[template_id_2]:
                                continue
                        
                        sub_info_2 = sub_infos_2[aac_base_2]
                        res_2 = sub_info_2[0]
                        residue_template_2 = sub_info_2[1]
                        if residue_template_2 == '-' or not res_2 in coordinate_map[target_chain_2][0]:
                            continue

                        if target_chain == target_chain_2 and res == res_2:
                            continue #No Calculations for the same annotation, this can be eventually be need in the future. If this is the case, just comment the clause.

                        #search in homomer chains, if there is a homomer mutation with shorter distance
                        min_min_d = None
                        for homo_chain in oligos_2:
                            if not homo_chain in coordinate_map:
                                continue
                            if not res_2 in coordinate_map[homo_chain][0]:
                                continue #Sanity Check - residue-id must be the same in the homo-chain
                            if coordinate_map[homo_chain][0][res_2][0] != coordinate_map[target_chain_2][0][res_2][0]:
                                continue #Sanity Check - residue type of the homomer residue must be the same as the residue type of the target residue
                            min_d,atom,atom2 = getMinDist(coordinate_map,res,target_chain,res_2,homo_chain)
                            if min_min_d == None or min_d < min_min_d:
                                min_min_d = min_d
                                min_atom = atom
                                min_atom2 = atom2
                                min_chain = homo_chain
                        if min_min_d != None:
                            anno_anno_dict[(template_id,template_id_2,aac_base,aac_base_2)] = (min_min_d,(min_atom,min_atom2),target_chain,min_chain)
                      

        used_template_ids.add(template_id)

    #print pdb_id,t_m_amount,anno_amount

    #print annotation_dict
    return annotation_dict,anno_anno_dict,errorlist


#Return True, if the template is not good
def weakCriteria(seq,res,rel_aln_len,seq_thresh,res_thresh,rel_aln_len_thresh):
    #print seq,seq_thresh
    #print res,res_thresh
    if float(seq) <= 1.0:
        seq = float(seq)*100.0
    if float(seq) < seq_thresh:
        return True
    if float(res) > res_thresh:
        return True
    else:
        return False

#template-structure: [pdb_id,seq_id,chain,aln_length,resolution,ligands,r-value,templateSelectionScore]
def good(template,seq_threshold,resolution_threshold,cov_threshold):
    seq_id = template[1]
    resolution = template[4]
    rel_aln_len = template[3]
    if weakCriteria(seq_id,resolution,rel_aln_len,seq_threshold,resolution_threshold,cov_threshold):
        return False
    return True

#called by serializedPipeline
def filterTemplates(templates,seq_threshold,resolution_threshold,cov_threshold):  
    filtered_templates = []   
    for template in templates:
        if good(template,seq_threshold,resolution_threshold,cov_threshold):
            filtered_templates.append(template)
    return filtered_templates

