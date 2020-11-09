import os
import sys
import traceback
import subprocess
import math
import time

from Bio.PDB import *
from Bio.SubsMat import MatrixInfo
from operator import itemgetter

import pdbParser
import spherecon

import rin
import database
import sdsc
import ray
import psutil

try:
    from memory_profiler import profile
except:
    pass

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
hydropathy = {'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,'W':-0.9,'Y':-1.3,
                'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5,'X':0.0,'B':0.0,'Z':0.0}

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

def calcDSSP(path,DSSP,angles=False,verbosity_level=0):

    dssp_dict = {}
    errorlist = []

    #print([DSSP,path])
    try:
        p = subprocess.Popen([DSSP,path], universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
        out, err = p.communicate() 
    except:
        out = ''
        err = '%s %s' % (str(sys.exc_info()[0]),str(sys.exc_info()[1]))
    # Alert user for errors 
    if err.strip():
        if not out.strip():
            if verbosity_level >= 4:
                print('DSSP failed to produce an output\n%s\n%s\n' % (path,err))
            errorlist.append('DSSP failed to produce an output\n%s\n%s\n' % (path,err))
            return dssp_dict,errorlist

    lines = out.split('\n')
    
    i = 1
    for line in lines:
        words = line.split()
        if words[0] == '#':
            break
        i += 1


    for line in lines[i:]:
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
            elif aa_type_one_letter.islower(): #This happens for SS-bond cysteines
                aa = Polypeptide.one_to_three('C')
                macc = residue_max_acc['Sander'][aa]
            else:
                aa = Polypeptide.one_to_three(aa_type_one_letter)
                macc = residue_max_acc['Sander'][aa]
            racc = acc/macc
        except:
            errorlist.append('DSSP failed for: \n%s\nChain: %s, Residue: %s,%s\n' % (path,chain,res,aa_type_one_letter))
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

    ssbond_map = {}
    link_map = {}
    cis_conformation_map = {}
    cis_follower_map = {}

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

                if record_name == 'SSBOND':
                    chain_1 = line[15]
                    res_nr_1 = line[17:22].replace(" ","")
                    chain_2 = line[29]
                    res_nr_2 = line[31:36].replace(" ","")
                    ssbond_len = float(line[73:78].replace(' ',''))
                    ssbond_map[(chain_1,res_nr_1)] = (chain_2,res_nr_2,ssbond_len)
                    ssbond_map[(chain_2,res_nr_2)] = (chain_1,res_nr_1,ssbond_len)

                if record_name == 'LINK  ':
                    altLoc_1 = line[16]
                    altLoc_2 = line[46]
                    if firstAltLoc == None and altLoc_1 != ' ':
                        firstAltLoc = altLoc_1 #The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                    if altLoc_1 != ' ' and altLoc_1 != firstAltLoc: #Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                        continue
                    if altLoc_2 != ' ' and altLoc_2 != firstAltLoc: #Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                        continue

                    atom_1 = line[12:16].replace(" ","")
                    res_name_1 = line[17:20].replace(" ","")
                    res_nr_1 = line[22:27].replace(" ","")
                    chain_1 = line[21]

                    atom_2 = line[42:46].replace(" ","")
                    res_name_2 = line[47:50].replace(" ","")
                    res_nr_2 = line[52:57].replace(" ","")
                    chain_2 = line[51]

                    link_dist = float(line[73:78].replace(" ",""))

                    link_map[(chain_1,res_nr_1)] = (atom_1,res_name_1,atom_2,res_name_2,res_nr_2,chain_2,link_dist)
                    link_map[(chain_2,res_nr_2)] = (atom_2,res_name_2,atom_1,res_name_1,res_nr_1,chain_1,link_dist)

                if record_name == 'CISPEP':
                    res_name_1 = line[11:14].replace(" ","")
                    chain_1 = line[15]
                    res_nr_1 = line[17:22].replace(" ","")
                    
                    res_name_2 = line[25:28].replace(" ","")
                    chain_2 = line[29]
                    res_nr_2 = line[31:36].replace(" ","")

                    angle = float(line[53:59].replace(" ",""))

                    cis_conformation_map[(chain_1,res_nr_1)] = (res_name_2,chain_2,res_nr_2,angle)
                    cis_follower_map[(chain_2,res_nr_2)] = (res_name_1,chain_1,res_nr_1,angle)

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
                        if res_name in pdbParser.threeToOne or not pdbParser.boring(res_name):#For a hetero peptide 'boring' hetero amino acids are allowed as well as other non boring molecules not in threeToOne, which are hopefully some kind of anormal amino acids
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
                    
                    if (chain_id,res_nr) in modres_map or res_name in pdbParser.threeToOne: #If it is a modified residue, than add it to the normal residues...
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
                            
                            ### LOOK HERE

                        if pdbParser.threeToOne[res_name][0] in pdbParser.oneToThree:
                            siss_het_res_name = pdbParser.oneToThree[pdbParser.threeToOne[res_name][0]]
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
                        if res_name in sdsc.metal_atoms:
                            metals.add((chain_id,res_nr))
                        elif res_name in sdsc.ion_atoms:
                            ions.add((chain_id,res_nr))
                        else:
                            ligands.add((chain_id,res_nr))

    return (coordinate_map,siss_map,res_contig_map,ligands,metals,ions,box_map,
            chain_type_map,b_factors,modres_map,ssbond_map,link_map,cis_conformation_map,cis_follower_map)

def dist(coord1,coord2):
    diff = [coord1[0]-coord2[0],coord1[1]-coord2[1],coord1[2]-coord2[2]]
    d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
    return d

def getMinSubDist(c_map,fuzzy_dist_matrix,target_res_id,res_chain,chain,config,distance_threshold=10.0):

    top20 = []
    top = 20
    for res in c_map[chain][0]:
        if chain == res_chain and res == target_res_id:
            continue
        if (res_chain,chain,target_res_id,res) in fuzzy_dist_matrix:
            d = fuzzy_dist_matrix[(res_chain,chain,target_res_id,res)]
        else:
            continue

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

def calcFuzzyDM(coordinate_map,box_map,config,calc_exact_distances = False):
    #coordinate_map: {Chain:[{Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]},{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}
    fuzzy_dm = {}
    processed_chains = set()
    for chain in coordinate_map:
        for chain_2 in coordinate_map:
            if chain_2 in processed_chains:
                continue
            neighbors,center_dist = box_check(box_map[chain],box_map[chain_2],distance_threshold = config.short_distance_threshold)
            if not neighbors:
                #fuzzy_dm[(chain,chain_2)] = center_dist
                #fuzzy_dm[(chain_2,chain)] = center_dist
                #print chain,chain_2
                continue

            for res in coordinate_map[chain][0]:
                test_coord = list(coordinate_map[chain][0][res][1].values())[0][1:]
                for res_2 in coordinate_map[chain_2][0]:
                    if chain == chain_2 and res == res_2:
                        continue
                    test_coord_2 = list(coordinate_map[chain_2][0][res_2][1].values())[0][1:]
                    d = dist(test_coord,test_coord_2)
                    if d > 2.*config.short_distance_threshold and d > 2.*config.milieu_threshold:
                        continue
                    if calc_exact_distances:
                        d,atom,atom2 = getMinDist(coordinate_map,res,chain,res_2,chain_2)
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

@ray.remote(num_cpus=1)
def analysis_chain_remote_wrapper(target_chain,analysis_dump):
    #hack proposed by the devs of ray to prevent too many processes being spawned
    resources = ray.ray.get_resource_ids() 
    cpus = [v[0] for v in resources['CPU']]
    psutil.Process().cpu_affinity(cpus)

    return analysis_chain(target_chain,analysis_dump)

def analysis_chain(target_chain,analysis_dump):
    (pdb_id,config,target_residues,siss_coord_map,centroid_map,
    res_contig_map,coordinate_map,fuzzy_dist_matrix,chain_type_map,
    b_factors,modres_map,ssbond_map,link_map,cis_conformation_map,cis_follower_map,
    profiles,milieu_dict,dssp,dssp_dict) = analysis_dump
    errorlist = []

    siss_map = {}
    try:
        siss_targets = {}
        siss_targets[target_chain] = list(target_residues[target_chain].keys())
        dist_matrix,res_res_angle_map = spherecon.calcDistMatrix(siss_coord_map,centroid_map,siss_targets,False)

        siss_map[target_chain] = spherecon.calculateSiss(siss_coord_map,centroid_map,dist_matrix,res_res_angle_map,siss_targets,False)[target_chain]

        del dist_matrix
    except:

        errorlist.append("Siss error: %s,%s" % (pdb_id,target_chain))

    structural_analysis_dict = {}

    milieu_dict[target_chain] = {}
    for target_res_id in target_residues[target_chain]:

        milieu_dict[target_chain][target_res_id] = {}
        three_letter = res_contig_map[target_chain][target_res_id][1]
        if three_letter in pdbParser.threeToOne:
            if pdbParser.threeToOne[three_letter][0] in pdbParser.oneToThree:
                one_letter = pdbParser.threeToOne[three_letter][0]
            else:
                one_letter = 'X'
        else:
            one_letter = 'X'
        lig_dists = {}
        min_chain_dists = {}

        inter_chain_median_kd = None
        inter_chain_dist_weighted_kd = None
        inter_chain_median_rsa = None
        inter_chain_dist_weighted_rsa = None
        intra_chain_median_kd = None
        intra_chain_dist_weighted_kd = None
        intra_chain_median_rsa = None
        intra_chain_dist_weighted_rsa = None

        #Distance Calculations
        for chain in coordinate_map:

            #Sub - Chain - Calculations
            if chain != target_chain and len(coordinate_map[chain][0]) > 0:
                min_chain_dist,min_res,atom_sub,atom_chain,inside_sphere = getMinSubDist(coordinate_map,fuzzy_dist_matrix,target_res_id,
                                                                                        target_chain,chain,config,distance_threshold=config.milieu_threshold)
                min_chain_dists[chain] = (min_chain_dist,(atom_chain,atom_sub),min_res)
                if chain_type_map[chain] == 'Protein' and config.neighborhood_calculation: #Filter DNA and RNA chains
                    inter_chain_kds = []
                    inter_chain_dist_weighted_kds = []
                    inter_chain_rsas = []
                    inter_chain_dist_weighted_rsas = []
                    total_dist_weights = 0.
                    for inter_chain_res in inside_sphere:
                        if not chain in dssp_dict:
                            continue
                        if not inter_chain_res in dssp_dict[chain]:
                            continue
                        dist = inside_sphere[inter_chain_res]
                        total_dist_weights += 1/dist
                        inter_chain_res_one_letter = pdbParser.threeToOne[res_contig_map[chain][inter_chain_res][1]][0]
                        kd = hydropathy[inter_chain_res_one_letter]
                        inter_chain_kds.append(kd)
                        inter_chain_dist_weighted_kds.append(kd/dist)
                        rsa = dssp_dict[chain][inter_chain_res][0]
                        inter_chain_rsas.append(rsa)
                        inter_chain_dist_weighted_rsas.append(rsa/dist)
                    if total_dist_weights > 0:
                        inter_chain_median_kd = database.median(inter_chain_kds)
                        inter_chain_dist_weighted_kd = sum(inter_chain_dist_weighted_kds)/total_dist_weights
                        inter_chain_median_rsa = database.median(inter_chain_rsas)
                        inter_chain_dist_weighted_rsa = sum(inter_chain_dist_weighted_rsas)/total_dist_weights

            #Residue-Residue Calculations
            elif chain == target_chain and len(coordinate_map[chain][0]) > 0 and config.neighborhood_calculation:
                if chain_type_map[chain] == 'Protein': #Filter DNA and RNA chains
                    if not chain in dssp_dict:
                        continue
                    min_chain_dist,min_res,atom_sub,atom_chain,inside_sphere = getMinSubDist(coordinate_map,fuzzy_dist_matrix,target_res_id,
                                                                                            target_chain,chain,config,distance_threshold=config.milieu_threshold)
                    intra_chain_kds = []
                    intra_chain_dist_weighted_kds = []
                    intra_chain_rsas = []
                    intra_chain_dist_weighted_rsas = []
                    total_dist_weights = 0.

                    if config.verbosity >= 5:
                        print('Residue-Residue calc:',pdb_id,chain,target_res_id,'inside the sphere:',len(inside_sphere))

                    for intra_chain_res in inside_sphere:
                        if not intra_chain_res in dssp_dict[chain]:
                            continue
                        dist = inside_sphere[intra_chain_res]
                        total_dist_weights += 1/dist
                        intra_chain_res_one_letter = pdbParser.threeToOne[res_contig_map[chain][intra_chain_res][1]][0]
                        kd = hydropathy[intra_chain_res_one_letter]
                        intra_chain_kds.append(kd)
                        intra_chain_dist_weighted_kds.append(kd/dist)
                        rsa = dssp_dict[chain][intra_chain_res][0]
                        intra_chain_rsas.append(rsa)
                        intra_chain_dist_weighted_rsas.append(rsa/dist)
                    if total_dist_weights > 0:
                        intra_chain_median_kd = database.median(intra_chain_kds)
                        intra_chain_dist_weighted_kd = sum(intra_chain_dist_weighted_kds)/total_dist_weights
                        intra_chain_median_rsa = database.median(intra_chain_rsas)
                        intra_chain_dist_weighted_rsa = sum(intra_chain_dist_weighted_rsas)/total_dist_weights

            #Sub - Lig - Calculations
            for hetres in coordinate_map[chain][1]:
                min_d,atom,atom2 = getMinDist(coordinate_map,target_res_id,target_chain,hetres,chain,het2=1)
                #only return ligand distances that are inside a given threshold
                if min_d > config.ligand_interest_sphere:
                    continue
                abr = coordinate_map[chain][1][hetres][0]
                ligand_identifier = "%s_%s_%s" % (abr,hetres,chain)
                lig_dists[ligand_identifier] = (min_d,(atom,atom2))
        
        #coordinate_map: {Chain:[{Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]},{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}
        
        #Homomer Sub - Sub Calculations:
        homomer_map = {}
        '''
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
        '''
        
        if dssp:
            if target_chain in dssp_dict:
                if target_res_id in dssp_dict[target_chain]:
                    (rsa,ssa,phi,psi) = dssp_dict[target_chain][target_res_id]
                else:
                    #print(target_res_id)
                    #print(siss_map)
                    siss_value = None
                    if target_res_id in siss_map:
                        siss_value = siss_map[target_chain][target_res_id]

                    ssa = None 
                    rsa  = siss_value
                    phi = None
                    psi = None
            else:
                errorlist.append(("dssp error: chain not in dssp_dict; %s; %s" % (pdb_id,target_chain)))
                dssp = False
                if target_res_id in siss_map:
                    siss_value = siss_map[target_chain][target_res_id]
                else:
                    siss_value = None
                ssa = None 
                rsa  = siss_value
                phi = None
                psi = None
        else:
            siss_value = None
            if target_res_id in siss_map:
                siss_value = siss_map[target_chain][target_res_id]
                    
            ssa = None
            phi = None
            psi = None
            rsa  = siss_value
        if target_chain in profiles:
            if target_res_id in profiles[target_chain]:
                profile,centrality_scores = profiles[target_chain][target_res_id]
            else:
                profile = None
                centrality_scores = None
        else:
            profile = None
            centrality_scores = None

        avg_b_factor = sum(b_factors[target_chain][target_res_id])/float(len(b_factors[target_chain][target_res_id]))

        if (target_chain,target_res_id) in modres_map:
            modres = True
        else:
            modres = False

        if (target_chain,target_res_id) in ssbond_map:
            (chain_2,res_nr_2,ssbond_len) = ssbond_map[(target_chain,target_res_id)]
            intra_ssbond = target_chain == chain_2
        else:
            intra_ssbond = None
            ssbond_len = None

        if (target_chain,target_res_id) in link_map:
            (atom_1,res_name_1,atom_2,res_name_2,res_nr_2,chain_2,link_dist) = link_map[(target_chain,target_res_id)]
            intra_link = target_chain == chain_2
        else:
            intra_link = None
            link_dist = None

        if (target_chain,target_res_id) in cis_conformation_map:
            (res_name_2,chain_2,res_nr_2,angle) = cis_conformation_map[(target_chain,target_res_id)]
            cis_conformation = angle
        else:
            cis_conformation = None

        if (target_chain,target_res_id) in cis_follower_map:
            (res_name_2,chain_2,res_nr_2,angle) = cis_follower_map[(target_chain,target_res_id)]
            cis_follower = angle
        else:
            cis_follower = None

        residue = sdsc.Residue(target_res_id,aa = one_letter,lig_dists = lig_dists,chain_distances = min_chain_dists,RSA = rsa,
                    SSA = ssa,homomer_distances = homomer_map,interaction_profile = profile,centralities = centrality_scores,
                    modres = modres,b_factor = avg_b_factor,phi = phi, psi = psi, intra_ssbond = intra_ssbond, ssbond_length = ssbond_len,
                    intra_link = intra_link, link_length = link_dist, cis_conformation = cis_conformation, cis_follower = cis_follower,
                    inter_chain_median_kd = inter_chain_median_kd, inter_chain_dist_weighted_kd = inter_chain_dist_weighted_kd,
                    inter_chain_median_rsa = inter_chain_median_rsa, inter_chain_dist_weighted_rsa = inter_chain_dist_weighted_rsa,
                    intra_chain_median_kd = intra_chain_median_kd, intra_chain_dist_weighted_kd = intra_chain_dist_weighted_kd,
                    intra_chain_median_rsa = intra_chain_median_rsa, intra_chain_dist_weighted_rsa = intra_chain_dist_weighted_rsa)

        structural_analysis_dict[target_res_id] = residue

    return target_chain,structural_analysis_dict,errorlist

#called by serializedPipeline
#@profile
def structuralAnalysis(pdb_id,config,target_dict = None):

    calculate_interaction_profiles = config.calculate_interaction_profiles
    dssp = config.dssp
    dssp_path = config.dssp_path
    pdb_path = config.pdb_path
    rin_db_path = config.rin_db_path
    verbosity = config.verbosity

    if verbosity >= 4:
        t0 = time.time()
        print('Start structuralAnalysis of:',pdb_id)

    page,fixed_10k_bug,path = pdbParser.standardParsePDB(pdb_id,pdb_path,return_10k_bool = True,get_is_local = True)

    if verbosity >= 4:
        t1 = time.time()
        print('Time for structuralAnalysis Part 1:',t1-t0)

    if page == '':
        print("Error while parsing: ",pdb_id)
        return {},{},[]

    (coordinate_map,siss_coord_map,res_contig_map,ligands,metals,ions,box_map,chain_type_map,
        b_factors,modres_map,ssbond_map,link_map,cis_conformation_map,cis_follower_map) = parsePDB(page)

    if verbosity >= 4:
        t2 = time.time()
        print('Time for structuralAnalysis Part 2:',t2-t1)

    if target_dict == None:
        target_residues = {}
        for chain in siss_coord_map:
            target_residues[chain] = list(siss_coord_map[chain].keys())
    else:
        target_residues = target_dict

    centroid_map = spherecon.calcCentroidMap(siss_coord_map,target_residues,False)

    if verbosity >= 4:
        t3 = time.time()
        print('Time for structuralAnalysis Part 3:',t3-t2)

    fuzzy_dist_matrix = calcFuzzyDM(coordinate_map,box_map,config)

    if verbosity >= 4:
        t4 = time.time()
        print('Time for structuralAnalysis Part 4:',t4-t3)

    errorlist = []
    structural_analysis_dict = {}

    milieu_dict = {}

    if dssp:
        #write temp pdb file only if we had to fix the 10k atom bug or the file is not stored locally
        if fixed_10k_bug or path == None:
            tmp_path = '%s/tmp_%s.pdb' % (config.temp_folder,pdb_id)
            f = open(tmp_path,'w')
            f.write(page)
            f.close()
        else:
            tmp_path = path

        #call DSSP  -structure of dssp_dict: {Chain:{Res_id:(racc,ssa)}}
        dssp_dict,errlist = calcDSSP(tmp_path,dssp_path,angles = True,verbosity_level=verbosity)
        errorlist += errlist
        #remove tmp_file
        if fixed_10k_bug:
            try:
                os.remove(tmp_path)
            except:
                pass
        if dssp_dict == {}:
            dssp = False

    if verbosity >= 4:
        t5 = time.time()
        print('Time for structuralAnalysis Part 5:',t5-t4)

    profiles = {}
    ligand_profiles = {}
    metal_profiles = {}
    ion_profiles = {}
    chain_chain_profiles = {}

    if calculate_interaction_profiles:
        (profiles,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles) = rin.lookup(pdb_id,page,config,None,None,ligands,
                                                                                                    metals,ions,res_contig_map,
                                                                                                    rin_db_path,chain_type_map,encoded = False)

    if verbosity >= 4:
        t6 = time.time()
        print('Time for structuralAnalysis Part 6:',t6-t5)

    if target_dict == None:
        target_residues = res_contig_map

    parent_dir = '/wibicom/SHARED_DATA/agress/structman/lib'
    os.environ["PYTHONPATH"] = parent_dir + ":" + os.environ.get("PYTHONPATH", "")

    parent_dir = '/wibicom/SHARED_DATA/agress/structman'
    os.environ["PYTHONPATH"] = parent_dir + ":" + os.environ.get("PYTHONPATH", "")

    if len(target_residues) > 10:

        analysis_dump = ray.put((pdb_id,config,target_residues,siss_coord_map,centroid_map,
                        res_contig_map,coordinate_map,fuzzy_dist_matrix,chain_type_map,
                        b_factors,modres_map,ssbond_map,link_map,cis_conformation_map,cis_follower_map,
                        profiles,milieu_dict,dssp,dssp_dict))

        core_results = []

        for target_chain in target_residues:
            if chain_type_map[target_chain] != 'Protein':
                continue

            core_results.append(analysis_chain_remote_wrapper.remote(target_chain,analysis_dump))

        core_outs = ray.get(core_results)
        for target_chain,chain_structural_analysis_dict,chain_errorlist in core_outs:
            structural_analysis_dict[target_chain] = chain_structural_analysis_dict
            errorlist += chain_errorlist
    else:
        analysis_dump = (pdb_id,config,target_residues,siss_coord_map,centroid_map,
                        res_contig_map,coordinate_map,fuzzy_dist_matrix,chain_type_map,
                        b_factors,modres_map,ssbond_map,link_map,cis_conformation_map,cis_follower_map,
                        profiles,milieu_dict,dssp,dssp_dict)
        for target_chain in target_residues:
            if chain_type_map[target_chain] != 'Protein':
                continue

            target_chain,chain_structural_analysis_dict,chain_errorlist = analysis_chain(target_chain,analysis_dump)
            structural_analysis_dict[target_chain] = chain_structural_analysis_dict
            errorlist += chain_errorlist

    if verbosity >= 4:
        t7 = time.time()
        print('Time for structuralAnalysis Part 7:',t7-t6)

    return structural_analysis_dict,errorlist,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles


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

