import os
import sys
import traceback
import subprocess
import math

from Bio.PDB import *
from Bio.SubsMat import MatrixInfo
from operator import itemgetter

import pdbParser as pdb
import spherecon as siss

import rin
import database

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
def qualityScore(resolution,coverage,seq_id,resolution_wf=0.25,coverage_wf=0.5,seq_id_wf=1.0):
    seq_id = float(seq_id)
    resolution = float(resolution)
    coverage = float(coverage)
    if seq_id > 1.0:
        seq_id = seq_id/100.0
    if coverage > 1.0:
        coverage = coverage/100.0
    
    #project the criteria to [0,1] via logistic regression
    resolution_value = (1+math.exp((1.5*resolution)-4))**(-1)

    seq_value = (1+math.exp(10*(0.4-seq_id)))**(-1)

    ws = sum((resolution_wf,coverage_wf,seq_id_wf))
    #print resolution_value,coverage,seq_id,r_value
    quality = sum((resolution_value*resolution_wf,coverage*coverage_wf,seq_value*seq_id_wf))/ws
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

    #print([DSSP,path])
    p = subprocess.Popen([DSSP,path], universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    out, err = p.communicate() 

    # Alert user for errors 
    if err.strip():
        if not out.strip():
            #errorlist.append(('DSSP failed to produce an output\n%s\n%s\n' % (path,err),'','',''))
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
            if aa_type_one_letter == 'X':
                macc = 220.0
            else:
                aa = Polypeptide.one_to_three(aa_type_one_letter)
                macc = residue_max_acc['Sander'][aa]
            racc = acc/macc
        except:
            errorlist.append(('DSSP failed for: \n%s\nChain: %s, Residue: %s\n' % (path,chain,res),'','',''))
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

def parsePDB(input_page):
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

    chain_type_map = {}

    box_map = {}

    ligands = set()
    metals = set()
    ions = set()
    b_factors = {}

    firstAltLoc = None

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
            
                altLoc = line[16]
                if firstAltLoc == None and altLoc != ' ':
                    firstAltLoc = altLoc #The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                if altLoc != ' ' and altLoc != firstAltLoc: #Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                    continue

            if len(line) > 60:
                try:
                    b_factor = float(line[61:67].replace(" ",""))
                except:
                    b_factor = 0.0
            if record_name == "ATOM" or record_name == 'MODRES' or record_name == "HETATM":
                if chain_id not in chain_type_map:
                    if record_name == "ATOM":
                        if len(res_name) == 1:
                            chain_type = "RNA"
                        elif len(res_name) == 2:
                            chain_type = "DNA"
                        elif len(res_name) == 3:
                            chain_type = "Protein"
                        chain_type_map[chain_id] = chain_type
                        
                    elif record_name == 'HETATM':
                        if res_name in pdb.threeToOne or not pdb.boring(res_name):#For a hetero peptide 'boring' hetero amino acids are allowed as well as other non boring molecules not in threeToOne, which are hopefully some kind of anormal amino acids
                            chain_type = 'Peptide'
                            chain_type_map[chain_id] = chain_type
                elif record_name == "ATOM" and chain_type_map[chain_id] == 'Peptide':
                    if len(res_name) == 1:
                        chain_type = "RNA"
                    elif len(res_name) == 2:
                        chain_type = "DNA"
                    elif len(res_name) == 3:
                        chain_type = "Protein"
                    chain_type_map[chain_id] = chain_type

            if record_name == "ATOM":
                if len(line) > 50 and not atom_name[0] in ('H','D'):
                    x = float(line[30:38].replace(" ",""))
                    y = float(line[38:46].replace(" ",""))
                    z = float(line[46:54].replace(" ",""))
                    if not chain_id in coordinate_map:
                        coordinate_map[chain_id] = [{},{}]
                        box_map[chain_id] = [x,x,y,y,z,z]
                    if res_nr not in coordinate_map[chain_id][0]:
                        coordinate_map[chain_id][0][res_nr] = [res_name,{}]
                    coordinate_map[chain_id][0][res_nr][1][atom_nr] = (atom_name,x,y,z)

                    if not chain_id in b_factors:
                        b_factors[chain_id] = {}
                    if not res_nr in b_factors[chain_id]:
                        b_factors[chain_id][res_nr] = []
                    b_factors[chain_id][res_nr].append(b_factor)

                    if x < box_map[chain_id][0]:
                        box_map[chain_id][0] = x
                    if x > box_map[chain_id][1]:
                        box_map[chain_id][1] = x
                    if y < box_map[chain_id][2]:
                        box_map[chain_id][2] = y
                    if y > box_map[chain_id][3]:
                        box_map[chain_id][3] = y
                    if z < box_map[chain_id][4]:
                        box_map[chain_id][4] = z
                    if z > box_map[chain_id][1]:
                        box_map[chain_id][4] = z

                    if not chain_id in siss_map:
                        siss_map[chain_id] = {}
                    if not res_nr in siss_map[chain_id]:
                        siss_map[chain_id][res_nr] = [res_name,{}]
                    siss_map[chain_id][res_nr][1][atom_nr] = (atom_name,x,y,z)
                    
                    if not chain_id in res_contig_map:
                        res_contig_map[chain_id] = {res_nr:[1,res_name]}
                        contig_help_map[chain_id] = 1
                    elif not res_nr in res_contig_map[chain_id]:
                        contig_help_map[chain_id] += 1
                        res_contig_map[chain_id][res_nr] = [contig_help_map[chain_id],res_name]

            if record_name == "HETATM":
                if len(line) > 50:
                    x = float(line[30:38].replace(" ",""))
                    y = float(line[38:46].replace(" ",""))
                    z = float(line[46:54].replace(" ",""))
                    if not chain_id in coordinate_map:
                        coordinate_map[chain_id] = [{},{}]
                        box_map[chain_id] = [x,x,y,y,z,z]
                    
                    if (chain_id,res_nr) in modres_map or res_name in pdb.threeToOne: #If it is a modified residue, than add it to the normal residues...
                        if atom_name[0] in ('H','D'):
                            continue
                        if not (chain_id,res_nr) in modres_map:
                            modres_map[(chain_id,res_nr)] = res_name
                        if res_nr not in coordinate_map[chain_id][0]:
                            coordinate_map[chain_id][0][res_nr] = [res_name,{}]
                        coordinate_map[chain_id][0][res_nr][1][atom_nr] = (atom_name,x,y,z)

                        if not chain_id in b_factors:
                            b_factors[chain_id] = {}
                        if not res_nr in b_factors[chain_id]:
                            b_factors[chain_id][res_nr] = []
                        b_factors[chain_id][res_nr].append(b_factor)

                        if x < box_map[chain_id][0]:
                            box_map[chain_id][0] = x
                        if x > box_map[chain_id][1]:
                            box_map[chain_id][1] = x
                        if y < box_map[chain_id][2]:
                            box_map[chain_id][2] = y
                        if y > box_map[chain_id][3]:
                            box_map[chain_id][3] = y
                        if z < box_map[chain_id][4]:
                            box_map[chain_id][4] = z
                        if z > box_map[chain_id][1]:
                            box_map[chain_id][4] = z

                        if pdb.threeToOne[res_name][0] in pdb.oneToThree:
                            siss_het_res_name = pdb.oneToThree[pdb.threeToOne[res_name][0]]
                        else:
                            siss_het_res_name = 'UNK'
                        if not chain_id in siss_map:
                            siss_map[chain_id] = {}
                        if not res_nr in siss_map[chain_id]:
                            siss_map[chain_id][res_nr] = [siss_het_res_name,{}]
                        siss_map[chain_id][res_nr][1][atom_nr] = (atom_name,x,y,z)

                        if not chain_id in res_contig_map:
                            res_contig_map[chain_id] = {res_nr:[1,res_name]}
                            contig_help_map[chain_id] = 1
                        elif not res_nr in res_contig_map[chain_id]:
                            contig_help_map[chain_id] += 1
                            res_contig_map[chain_id][res_nr] = [contig_help_map[chain_id],res_name]
                    else:    
                        if res_nr not in coordinate_map[chain_id][1]:#If not, then add it to the ligands
                            coordinate_map[chain_id][1][res_nr] = [res_name,{}]
                        coordinate_map[chain_id][1][res_nr][1][atom_nr] = (atom_name,x,y,z)
                        if res_name in database.metal_atoms:
                            metals.add((chain_id,res_nr))
                        elif res_name in database.ion_atoms:
                            ions.add((chain_id,res_nr))
                        else:
                            ligands.add((chain_id,res_nr))

    return coordinate_map,siss_map,res_contig_map,ligands,metals,ions,box_map,chain_type_map,b_factors,modres_map

def dist(coord1,coord2):
    diff = [coord1[0]-coord2[0],coord1[1]-coord2[1],coord1[2]-coord2[2]]
    d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
    return d

def getMinSubDist(c_map,fuzzy_dist_matrix,target_res_id,res_chain,chain,distance_threshold=10.0):

    top20 = []
    top = 20
    for res in c_map[chain][0]:
        if chain == res_chain and res == target_res_id:
            continue
        if (res_chain,chain,target_res_id,res) in fuzzy_dist_matrix:
            d = fuzzy_dist_matrix[(res_chain,chain,target_res_id,res)]
        else:
            center_dist = fuzzy_dist_matrix[(res_chain,chain)]
            return center_dist,'None','None','None',{}

        if len(top20) < top:
            top20.append((res,d))
        elif len(top20) == top:
            top20.append((res,d))
            top20 = sorted(top20,key = lambda x:x[1])
        elif d < top20[top][1]:
            top20[top] = (res,d)
            top20 = sorted(top20,key = lambda x:x[1])

    min_sub_d = None
    min_res = None
    min_atom_sub = None
    min_atom_chain = None
    atomlist = c_map[res_chain][0][target_res_id][1]
    inside_sphere = {}
    for (res,d) in top20:
        for atomnr in c_map[chain][0][res][1]:
            #exclude Hydrogens
            if c_map[chain][0][res][1][atomnr][0][0] != 'H':
                coord1 = c_map[chain][0][res][1][atomnr][1:]
                for atomnr2 in atomlist:
                    if atomlist[atomnr2][0][0] != 'H':
                        coord2 = atomlist[atomnr2][1:]
                        d = dist(coord1,coord2)
                        if d < 1.2:
                            continue
                        if d < distance_threshold:
                            if not res in inside_sphere:
                                inside_sphere[res] = d
                            elif d < inside_sphere[res]:
                                inside_sphere[res] = d
                        if min_sub_d == None or d < min_sub_d:
                            min_sub_d = d
                            min_res = res
                            min_atom_sub = atomnr2
                            min_atom_chain = atomnr

    return min_sub_d,min_res,min_atom_sub,min_atom_chain,inside_sphere

def box_check(box_1,box_2,distance_threshold=5.0):
    [min_x_1,max_x_1,min_y_1,max_y_1,min_z_1,max_z_1] = box_1
    [min_x_2,max_x_2,min_y_2,max_y_2,min_z_2,max_z_2] = box_2

    center_1 = [(min_x_1 + max_x_1)/2.,(min_y_1 + max_y_1)/2.,(min_z_1 + max_z_1)/2.]
    center_2 = [(min_x_2 + max_x_2)/2.,(min_y_2 + max_y_2)/2.,(min_z_2 + max_z_2)/2.]

    center_dist = dist(center_1,center_2)

    if center_dist < (max(max_x_1-min_x_1,max_y_1-min_y_1,max_z_1-min_z_1)/2.)+(max(max_x_2-min_x_2,max_y_2-min_y_2,max_z_2-min_z_2)/2.)+distance_threshold:
        return True,center_dist
    else:
        return False,center_dist

def calcFuzzyDM(coordinate_map,box_map):
    #coordinate_map: {Chain:[{Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]},{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}
    fuzzy_dm = {}
    processed_chains = set()
    for chain in coordinate_map:
        for chain_2 in coordinate_map:
            if chain_2 in processed_chains:
                continue
            neighbors,center_dist = box_check(box_map[chain],box_map[chain_2])
            if not neighbors:
                fuzzy_dm[(chain,chain_2)] = center_dist
                fuzzy_dm[(chain_2,chain)] = center_dist
                #print chain,chain_2
                continue

            for res in coordinate_map[chain][0]:
                test_coord = list(coordinate_map[chain][0][res][1].values())[0][1:]
                for res_2 in coordinate_map[chain_2][0]:
                    test_coord_2 = list(coordinate_map[chain_2][0][res_2][1].values())[0][1:]
                    d = dist(test_coord,test_coord_2)
                    fuzzy_dm[(chain,chain_2,res,res_2)] = d
                    fuzzy_dm[(chain_2,chain,res_2,res)] = d
        processed_chains.add(chain)
    return fuzzy_dm

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
        print(c_map,res,chain,res2,chain2,het,het2)

    return min_d,min_atom,min_atom2

#called by lite_pipeline
def liteAnalysis(pdb_id,chain_structure_map,pdb_path,dssp_path,rin_db_path,neighborhood_calculation=False,dssp=True,calculate_interaction_profiles=True,distance_threshold=5.0):
    page = pdb.standardParsePDB(pdb_id,pdb_path)
    resolution,homomer_dict = pdb.getInfo(pdb_id,pdb_path)

    if page == '':
        print("Error while parsing: ",pdb_id)
        return {},{},[]

    coordinate_map,siss_coord_map,res_contig_map,ligands,metals,ions,box_map,chain_type_map,b_factors,modres_map = parsePDB(page)

    fuzzy_dist_matrix = calcFuzzyDM(coordinate_map,box_map)
    errorlist = []
    annotations = {}

    for t_chain in homomer_dict:#Happens for repeated chains in asymetric units, which do not occur in the biological assembly
        irregular_homo_chains = []
        for i,homo_chain in enumerate(homomer_dict[t_chain]):
            if not homo_chain in chain_type_map:
                irregular_homo_chains.append(i)

        irregular_homo_chains = sorted(irregular_homo_chains,reverse=True)
        for i in irregular_homo_chains:
            del homomer_dict[t_chain][i]

    residue_residue_dict = {}

    if dssp:
        #write temp pdb file
        tmp_path = 'tmp_%s.pdb' % (pdb_id)
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

    profiles = {}
    ligand_profiles = {}
    metal_profiles = {}
    ion_profiles = {}
    chain_chain_profiles = {}
    if calculate_interaction_profiles:
        #print(rin.lookup)                                                               lookup(pdb_id,inp_residues,chains,ligands,metals,ions,res_contig_map,base_path,chain_type_map):
        profiles,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles = rin.lookup(pdb_id,None,None,ligands,metals,ions,res_contig_map,rin_db_path,chain_type_map)

    for target_chain in res_contig_map:
        if chain_type_map[target_chain] != 'Protein':
            continue
        
        target_residues = []
        if target_chain in chain_structure_map:
            structure,sub_infos = chain_structure_map[target_chain]
            oligos = structure['Oligo']
            for aacbase in sub_infos:
                res_id = sub_infos[aacbase][0]
                target_residues.append(res_id)
        else:
            continue

        #siss_coord_map: {Chain:{String:[String,{String:(String,float,float,float)}]}} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.
        try:
            centroid_map = siss.calcCentroidMap(siss_coord_map[target_chain],target_residues,False)

            dist_matrix = siss.calcDistMatrix(siss_coord_map[target_chain],centroid_map,target_residues,False)

            siss_map = siss.calculateSiss(siss_coord_map[target_chain],centroid_map,dist_matrix,target_residues,False)
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc()
            errorlist.append(("Siss error: %s,%s\n" % (pdb_id,target_chain),e,f,g))

        

        annotation_dict = {}

        
        #anno_amount = 0
        residue_residue_dict[target_chain] = {}
        for target_res_id in target_residues:

            residue_residue_dict[target_chain][target_res_id] = {}
            three_letter = res_contig_map[target_chain][target_res_id][1]
            if three_letter in pdb.threeToOne:
                if pdb.threeToOne[three_letter][0] in pdb.oneToThree:
                    one_letter = pdb.threeToOne[three_letter][0]
                else:
                    one_letter = 'X'
            else:
                one_letter = 'X'
            min_dists = {}
            min_chain_dists = {}
            #Distance Calculations
            for chain in coordinate_map:

                #Sub - Chain - Calculations
                if chain != target_chain and len(coordinate_map[chain][0]) > 0:
                    min_chain_dist,min_res,atom_sub,atom_chain,inside_sphere = getMinSubDist(coordinate_map,fuzzy_dist_matrix,target_res_id,target_chain,chain,distance_threshold=distance_threshold)
                    min_chain_dists[chain] = (min_chain_dist,(atom_chain,atom_sub),min_res)
                    if chain_type_map[chain] == 'Protein' and neighborhood_calculation: #Filter DNA and RNA chains
                        residue_residue_dict[target_chain][target_res_id][chain] = inside_sphere


                #Residue-Residue Calculations
                elif chain == target_chain and len(coordinate_map[chain][0]) > 0 and neighborhood_calculation:
                    if chain_type_map[chain] == 'Protein': #Filter DNA and RNA chains
                        min_chain_dist,min_res,atom_sub,atom_chain,inside_sphere = getMinSubDist(coordinate_map,fuzzy_dist_matrix,target_res_id,target_chain,chain,distance_threshold=distance_threshold)
                        residue_residue_dict[target_chain][target_res_id][chain] = inside_sphere

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
                homomer_map[homo_chain] = min_d        

            
            if dssp:
                if target_chain in dssp_dict:
                    if target_res_id in dssp_dict[target_chain]:
                        (raac,ssa) = dssp_dict[target_chain][target_res_id]
                    else:
                        #print(target_res_id)
                        #print(siss_map)
                        siss_value = None
                        if target_res_id in siss_map:
                            siss_value = siss_map[target_res_id]

                        ssa = None 
                        raac  = siss_value
                else:
                    [e,f,g] = sys.exc_info()
                    g = traceback.format_exc()
                    errorlist.append(("dssp error: chain not in dssp_dict; %s; %s" % (pdb_id,target_chain),e,f,g))
                    dssp = False
                    if target_res_id in siss_map:
                        siss_value = siss_map[target_res_id]
                    else:
                        siss_value = None
                    ssa = None 
                    raac  = siss_value
            else:
                siss_value = None
                if target_res_id in siss_map:
                    siss_value = siss_map[target_res_id]
                        
                ssa = None 
                raac  = siss_value
            if target_chain in profiles:
                if target_res_id in profiles[target_chain]:
                    profile,centrality_scores = profiles[target_chain][target_res_id]
                else:
                    profile = None
                    centrality_scores = (None,None,None,None,None,None)
            else:
                profile = None
                centrality_scores = (None,None,None,None,None,None)

            avg_b_factor = sum(b_factors[target_chain][target_res_id])/float(len(b_factors[target_chain][target_res_id]))

            if (target_chain,target_res_id) in modres_map:
                modres = True
            else:
                modres = False

            annotation_dict[target_res_id] = (one_letter,min_dists,min_chain_dists,raac,ssa,homomer_map,profile,centrality_scores,avg_b_factor,modres)

        annotations[target_chain] = annotation_dict

    return annotations,residue_residue_dict,errorlist,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles

#called by serializedPipeline
def structuralAnalysis(pdb_id,chain_structure_map,pdb_path,dssp_path,rin_db_path,neighborhood_calculation=False,dssp=True,calculate_interaction_profiles=True,distance_threshold=5.0):

    page = pdb.standardParsePDB(pdb_id,pdb_path)
    resolution,homomer_dict = pdb.getInfo(pdb_id,pdb_path)

    if page == '':
        print("Error while parsing: ",pdb_id)
        return {},{},[]

    coordinate_map,siss_coord_map,res_contig_map,ligands,metals,ions,box_map,chain_type_map,b_factors,modres_map = parsePDB(page)

    fuzzy_dist_matrix = calcFuzzyDM(coordinate_map,box_map)
    errorlist = []
    annotations = {}

    example_structure = chain_structure_map[list(chain_structure_map.keys())[0]][1]

    for t_chain in homomer_dict:#Happens for repeated chains in asymetric units, which do not occur in the biological assembly
        irregular_homo_chains = []
        for i,homo_chain in enumerate(homomer_dict[t_chain]):
            if not homo_chain in chain_type_map:
                irregular_homo_chains.append(i)

        irregular_homo_chains = sorted(irregular_homo_chains,reverse=True)
        for i in irregular_homo_chains:
            del homomer_dict[t_chain][i]

    interacting_chain_map = {}

    residue_residue_dict = {}

    if dssp:
        #write temp pdb file
        tmp_path = 'tmp_%s.pdb' % (pdb_id)
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

    profiles = {}
    ligand_profiles = {}
    metal_profiles = {}
    ion_profiles = {}
    chain_chain_profiles = {}
    if calculate_interaction_profiles:
        #print(rin.lookup)                                                               lookup(pdb_id,inp_residues,chains,ligands,metals,ions,res_contig_map,base_path,chain_type_map):
        profiles,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles = rin.lookup(pdb_id,None,None,ligands,metals,ions,res_contig_map,rin_db_path,chain_type_map)

    for target_chain in res_contig_map:
        if chain_type_map[target_chain] != 'Protein':
            continue
        

        if target_chain in chain_structure_map:
            s_id,structure = chain_structure_map[target_chain]
            #print structure
            oligos = structure['Oligo']
        else:
            if target_chain in homomer_dict:
                oligos = homomer_dict[target_chain]
            else:
                oligos = target_chain
            interacting_chain_map[target_chain] = (resolution,example_structure['IAP'],oligos)

        target_residues = list(res_contig_map[target_chain].keys())

        #if not dssp:
        #surface-calculation with SISI

        #siss_coord_map: {Chain:{String:[String,{String:(String,float,float,float)}]}} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.
        try:
            centroid_map = siss.calcCentroidMap(siss_coord_map[target_chain],target_residues,False)

            dist_matrix = siss.calcDistMatrix(siss_coord_map[target_chain],centroid_map,target_residues,False)

            siss_map = siss.calculateSiss(siss_coord_map[target_chain],centroid_map,dist_matrix,target_residues,False)
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc()
            errorlist.append(("Siss error: %s,%s\n" % (pdb_id,target_chain),e,f,g))

        

        annotation_dict = {}

        
        #anno_amount = 0
        residue_residue_dict[target_chain] = {}
        for target_res_id in res_contig_map[target_chain]:

            residue_residue_dict[target_chain][target_res_id] = {}
            three_letter = res_contig_map[target_chain][target_res_id][1]
            if three_letter in pdb.threeToOne:
                if pdb.threeToOne[three_letter][0] in pdb.oneToThree:
                    one_letter = pdb.threeToOne[three_letter][0]
                else:
                    one_letter = 'X'
            else:
                one_letter = 'X'
            min_dists = {}
            min_chain_dists = {}
            #Distance Calculations
            for chain in coordinate_map:

                #Sub - Chain - Calculations
                if chain != target_chain and len(coordinate_map[chain][0]) > 0:
                    min_chain_dist,min_res,atom_sub,atom_chain,inside_sphere = getMinSubDist(coordinate_map,fuzzy_dist_matrix,target_res_id,target_chain,chain,distance_threshold=distance_threshold)
                    min_chain_dists[chain] = (min_chain_dist,(atom_chain,atom_sub),min_res)
                    if chain_type_map[chain] == 'Protein' and neighborhood_calculation: #Filter DNA and RNA chains
                        residue_residue_dict[target_chain][target_res_id][chain] = inside_sphere


                #Residue-Residue Calculations
                elif chain == target_chain and len(coordinate_map[chain][0]) > 0 and neighborhood_calculation:
                    if chain_type_map[chain] == 'Protein': #Filter DNA and RNA chains
                        min_chain_dist,min_res,atom_sub,atom_chain,inside_sphere = getMinSubDist(coordinate_map,fuzzy_dist_matrix,target_res_id,target_chain,chain,distance_threshold=distance_threshold)
                        residue_residue_dict[target_chain][target_res_id][chain] = inside_sphere

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
                homomer_map[homo_chain] = min_d        

            
            if dssp:
                if target_chain in dssp_dict:
                    if target_res_id in dssp_dict[target_chain]:
                        (raac,ssa) = dssp_dict[target_chain][target_res_id]
                    else:
                        #print(target_res_id)
                        #print(siss_map)
                        siss_value = None
                        if target_res_id in siss_map:
                            siss_value = siss_map[target_res_id]

                        ssa = None 
                        raac  = siss_value
                else:
                    [e,f,g] = sys.exc_info()
                    g = traceback.format_exc()
                    errorlist.append(("dssp error: chain not in dssp_dict; %s; %s" % (pdb_id,target_chain),e,f,g))
                    dssp = False
                    if target_res_id in siss_map:
                        siss_value = siss_map[target_res_id]
                    else:
                        siss_value = None
                    ssa = None 
                    raac  = siss_value
            else:
                siss_value = None
                if target_res_id in siss_map:
                    siss_value = siss_map[target_res_id]
                        
                ssa = None 
                raac  = siss_value
            if target_chain in profiles:
                if target_res_id in profiles[target_chain]:
                    profile,centrality_scores = profiles[target_chain][target_res_id]
                else:
                    profile = None
                    centrality_scores = (None,None,None,None,None,None)
            else:
                profile = None
                centrality_scores = (None,None,None,None,None,None)

            avg_b_factor = sum(b_factors[target_chain][target_res_id])/float(len(b_factors[target_chain][target_res_id]))

            if (target_chain,target_res_id) in modres_map:
                modres = True
            else:
                modres = False

            annotation_dict[target_res_id] = (one_letter,min_dists,min_chain_dists,raac,ssa,homomer_map,profile,centrality_scores,avg_b_factor,modres)

        annotations[target_chain] = annotation_dict

    return annotations,interacting_chain_map,residue_residue_dict,errorlist,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles


#Return True, if the template is not good
def weakCriteria(seq,res,rel_aln_len,seq_thresh,res_thresh,rel_aln_len_thresh):

    if float(seq) <= 1.0:
        seq = float(seq)*100.0
    if float(seq) < seq_thresh:
        return True
    if float(res) > res_thresh:
        return True
    else:
        return False

#template-structure: [pdb_id,seq_id,chain,aln_length,resolution,ligands,r-value,templateSelectionScore]
def good(structure,seq_threshold,resolution_threshold,cov_threshold):
    seq_id = structure['Seq_Id']
    resolution = structure['Resolution']
    cov = structure['Coverage']
    if weakCriteria(seq_id,resolution,cov,seq_threshold,resolution_threshold,cov_threshold):
        return False
    return True

#called by serializedPipeline
def filterTemplates(structures,seq_threshold,resolution_threshold,cov_threshold):
    #print structures
    filtered_structures = {}
    for (pdb_id,chain) in structures:
        if good(structures[(pdb_id,chain)],seq_threshold,resolution_threshold,cov_threshold):
            filtered_structures[(pdb_id,chain)] = structures[(pdb_id,chain)]

    return filtered_structures

