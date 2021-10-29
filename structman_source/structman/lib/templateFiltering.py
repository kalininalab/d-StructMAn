import math
import os
import subprocess
import sys
import time
import traceback
import gc
import psutil
import ray
import pickle

from structman.lib import database, pdbParser, rin, sdsc, spherecon, serializedPipeline
from structman.utils import median, distance, pack, unpack, decompress
from structman.lib.sdsc.consts.residues import CORRECT_COUNT, THREE_TO_ONE
try:
    from memory_profiler import profile
except:
    pass

# called by templateSelection
def qualityScore(resolution, coverage, seq_id, resolution_wf=0.25, coverage_wf=0.5, seq_id_wf=1.0):
    seq_id = float(seq_id)
    resolution = float(resolution)
    coverage = float(coverage)
    if seq_id > 1.0:
        seq_id = seq_id / 100.0
    if coverage > 1.0:
        coverage = coverage / 100.0

    # project the criteria to [0,1] via logistic regression
    resolution_value = (1 + math.exp((1.5 * resolution) - 4))**(-1)

    seq_value = (1 + math.exp(10 * (0.4 - seq_id)))**(-1)

    ws = sum((resolution_wf, coverage_wf, seq_id_wf))
    # print resolution_value,coverage,seq_id,r_value
    quality = sum((resolution_value * resolution_wf, coverage * coverage_wf, seq_value * seq_id_wf)) / ws
    return quality


# called by database
def candidateScore(lig_sub_dist, chain_sub_dist, lig_wf=1.0, chain_wf=1.0, useBlosum=False, aac="", blosum_wf=0.5):
    if useBlosum:
        if aac == "":
            raise NameError("If useBlosum is True, an aac is needed")

        try:
            blosum_value = 0.6 - float(sdsc.BLOSUM62[(aac[0], aac[-1])]) / 10
        except:
            blosum_value = 0.6 - float(sdsc.BLOSUM62[(aac[-1], aac[0])]) / 10
        if blosum_value < 0.0:
            blosum_value = 0.0
        # print blosum_value
    else:
        blosum_value = 1.0
        blosum_wf = 0.0

    # project the criteria to [0,1] via logistic regression
    lig_value = (1 + math.exp(lig_sub_dist - 10))**(-1)
    chain_value = (1 + math.exp(chain_sub_dist - 10))**(-1)

    if lig_sub_dist == -1:
        lig_value = 0.0
    if chain_sub_dist == -1:
        chain_value = 0.0

    ws = sum((lig_wf, chain_wf, blosum_wf))

    candidate = sum((lig_value * lig_wf, chain_value * chain_wf, blosum_value * blosum_wf)) / ws
    return candidate


def calcDSSP(path, DSSP, angles=False, verbosity_level=0):

    dssp_dict = {}
    errorlist = []

    # print([DSSP,path])
    try:
        p = subprocess.Popen([DSSP, path], universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
    except:
        out = ''
        err = '%s %s' % (str(sys.exc_info()[0]), str(sys.exc_info()[1]))
    # Alert user for errors
    if err.strip():
        if not out.strip():
            if verbosity_level >= 4:
                print('DSSP failed to produce an output\n%s\n%s\n' % (path, err))
            errorlist.append('DSSP failed to produce an output\n%s\n%s\n' % (path, err))
            return dssp_dict, errorlist

    lines = out.split('\n')

    i = 1
    for line in lines:
        if line == '':
            continue
        words = line.split()
        if words[0] == '#':
            break
        i += 1

    for line in lines[i:]:
        if line == '':
            continue
        if len(line) < 46:
            continue
        if line[9] == " ":
            continue

        res = (line[5:11]).replace(" ", "")  # this includes the insertion code
        insertion_code = line[10]
        chain = line[11]
        aa_type_one_letter = line[13]
        ssa = line[16]
        try:
            acc = float(line[34:38])

            main_chain_acc = float(line[38:42])

            side_chain_acc = float(line[42:46])
        except:
            errorlist.append('DSSP float conversion failed for: \n%s\nChain: %s, Residue: %s,%s\n%s\n\n' % (path, chain, res, aa_type_one_letter, line))
            racc = None
            relative_main_chain_acc = None
            relative_side_chain_acc = None

        try:
            if aa_type_one_letter == 'X':
                macc = 220.0
            elif aa_type_one_letter.islower():  # This happens for SS-bond cysteines
                aa = sdsc.ONE_TO_THREE['C']
                macc = sdsc.RESIDUE_MAX_ACC['Sander'][aa]
            else:
                aa = sdsc.ONE_TO_THREE[aa_type_one_letter]
                macc = sdsc.RESIDUE_MAX_ACC['Sander'][aa]
            racc = acc / macc
            relative_main_chain_acc = main_chain_acc / macc
            relative_side_chain_acc = side_chain_acc / macc
        except:
            errorlist.append('DSSP (%s) failed for: \n%s\nChain: %s, Residue: %s,%s\n%s\n\n' % (DSSP, path, chain, res, aa_type_one_letter, line))
            racc = None
            relative_main_chain_acc = None
            relative_side_chain_acc = None

        try:
            if len(line) > 123:
                phi = float(line[111:117].replace(" ", ""))
                psi = float(line[117:123].replace(" ", ""))
        except:
            errorlist.append('DSSP failed for (angles): \n%s\nChain: %s, Residue: %s,%s\n%s\n\n' % (path, chain, res, aa_type_one_letter, line))
            phi = None
            psi = None

        if chain not in dssp_dict:
            dssp_dict[chain] = {}
        dssp_dict[chain][res] = (racc, relative_main_chain_acc, relative_side_chain_acc, ssa)
        if angles:
            dssp_dict[chain][res] = (racc, relative_main_chain_acc, relative_side_chain_acc, ssa, phi, psi)
    return dssp_dict, errorlist


def parsePDB(input_page):
    """
    Parses a PDB-file and takes all atomic coordinates.

    Input:
    input_page: String ; content of a pdb file

    Output:
    coordinate_map: {Chain:[{Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]},{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}
    """
    # siss_map: {String:[String,{String:(String,float,float,float)}]} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.

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
    chain_type = '-'
    chainlist = []

    peptide_count = {}

    box_map = {}

    ligands = set()
    metals = set()
    ions = set()
    rare_residues = set()
    b_factors = {}

    firstAltLoc = None

    for line in lines:
        if len(line) > 5:
            record_name = line[0:6].replace(" ", "")
            if record_name == "ENDMDL":
                # print len(coordinate_map)
                break
            elif record_name == 'SEQRES':
                for tlc in line[19:].split():
                    tlc = tlc.strip()
                    if len(tlc) != 3:
                        continue
                    if tlc not in sdsc.THREE_TO_ONE:
                        rare_residues.add(tlc)
        # ignore short lines
        if len(line) > 20:
            atom_nr = line[6:11].replace(" ", "")
            if record_name.count('ATOM') > 0 and record_name != 'ATOM':  # 100k atom bug fix
                atom_nr = '%s%s' % (record_name[4:], atom_nr)
                record_name = 'ATOM'
            atom_name = line[12:16].replace(" ", "")
            res_name = line[17:20].replace(" ", "")

            if len(line) > 21:
                if record_name == 'MODRES':
                    chain_id = line[16]
                    res_nr = line[18:23].replace(" ", "")
                    res_name = line[24:27].replace(" ", "")
                    if not (chain_id, res_nr) in modres_map:
                        modres_map[(chain_id, res_nr)] = res_name

                if record_name == 'SSBOND':
                    chain_1 = line[15]
                    res_nr_1 = line[17:22].replace(" ", "")
                    chain_2 = line[29]
                    res_nr_2 = line[31:36].replace(" ", "")
                    ssbond_len = float(line[73:78].replace(' ', ''))
                    ssbond_map[(chain_1, res_nr_1)] = (chain_2, res_nr_2, ssbond_len)
                    ssbond_map[(chain_2, res_nr_2)] = (chain_1, res_nr_1, ssbond_len)

                if record_name == 'LINK  ':
                    atom_1 = line[12:16].replace(" ", "")
                    res_name_1 = line[17:20].replace(" ", "")
                    res_nr_1 = line[22:27].replace(" ", "")
                    chain_1 = line[21]

                    atom_2 = line[42:46].replace(" ", "")
                    res_name_2 = line[47:50].replace(" ", "")
                    res_nr_2 = line[52:57].replace(" ", "")
                    chain_2 = line[51]

                    altLoc_1 = line[16]
                    altLoc_2 = line[46]
                    if firstAltLoc is None and altLoc_1 != ' ':
                        firstAltLoc = altLoc_1  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                    if altLoc_1 != ' ' and altLoc_1 != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                        continue
                    if altLoc_2 != ' ' and altLoc_2 != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                        continue

                    link_dist = float(line[73:78].replace(" ", ""))

                    link_map[(chain_1, res_nr_1)] = (atom_1, res_name_1, atom_2, res_name_2, res_nr_2, chain_2, link_dist)
                    link_map[(chain_2, res_nr_2)] = (atom_2, res_name_2, atom_1, res_name_1, res_nr_1, chain_1, link_dist)

                if record_name == 'CISPEP':
                    res_name_1 = line[11:14].replace(" ", "")
                    chain_1 = line[15]
                    res_nr_1 = line[17:22].replace(" ", "")

                    res_name_2 = line[25:28].replace(" ", "")
                    chain_2 = line[29]
                    res_nr_2 = line[31:36].replace(" ", "")

                    angle = float(line[53:59].replace(" ", ""))

                    cis_conformation_map[(chain_1, res_nr_1)] = (res_name_2, chain_2, res_nr_2, angle)
                    cis_follower_map[(chain_2, res_nr_2)] = (res_name_1, chain_1, res_nr_1, angle)

                chain_id = line[21]
                res_nr = line[22:27].replace(" ", "")  # [22:27] includes the insertion_code
                if record_name == "ATOM" or record_name == "HETATM":
                    altLoc = line[16]
                    if firstAltLoc is None and altLoc != ' ':
                        firstAltLoc = altLoc  # The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                    if altLoc != ' ' and altLoc != firstAltLoc:  # Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                        continue

            if len(line) > 60:
                try:
                    b_factor = float(line[61:67].replace(" ", ""))
                except:
                    b_factor = 0.0
            if record_name == "ATOM" or record_name == 'MODRES' or record_name == "HETATM":
                if chain_id not in chain_type_map:
                    chainlist.append(chain_id)
                    if record_name == "ATOM" or record_name == 'MODRES':
                        if len(res_name) == 1:
                            chain_type = "RNA"
                        elif len(res_name) == 2:
                            chain_type = "DNA"
                        elif len(res_name) == 3:
                            chain_type = "Protein"
                        chain_type_map[chain_id] = chain_type
                        peptide_count[chain_id] = [0, 0]
                    elif record_name == 'HETATM':
                        if res_name in sdsc.THREE_TO_ONE or not sdsc.boring(res_name) or (res_name in rare_residues):
                            chain_type_map[chain_id] = "Protein"
                            peptide_count[chain_id] = [0, 0]
                        elif len(res_name) == 1:
                            chain_type_map[chain_id] = "RNA"
                            peptide_count[chain_id] = [0, 0]
                        elif len(res_name) == 2:
                            chain_type_map[chain_id] = "DNA"
                            peptide_count[chain_id] = [0, 0]
                        else:
                            chain_type_map[chain_id] = 'MISC'


                if record_name == 'HETATM':
                    if res_name in sdsc.THREE_TO_ONE or not sdsc.boring(res_name) or (res_name in rare_residues):  # For a hetero peptide 'boring' hetero amino acids are allowed as well as other non boring molecules not in THREE_TO_ONE, which are hopefully some kind of anormal amino acids
                        peptide_count[chain_id][1] += 1
                    elif len(res_name) < 3:
                        peptide_count[chain_id][0] += 1

                elif record_name == "ATOM" or record_name == 'MODRES':
                    peptide_count[chain_id][0] += 1
                    if len(res_name) == 1:
                        chain_type_map[chain_id] = "RNA"
                    elif len(res_name) == 2:
                        chain_type_map[chain_id] = "DNA"

            if record_name == "ATOM":
                if len(line) > 50 and not atom_name[0] in ('H', 'D'):
                    x = float(line[30:38].replace(" ", ""))
                    y = float(line[38:46].replace(" ", ""))
                    z = float(line[46:54].replace(" ", ""))
                    if chain_id not in coordinate_map:
                        coordinate_map[chain_id] = [{}, {}]
                        box_map[chain_id] = [x, x, y, y, z, z]
                    if res_nr not in coordinate_map[chain_id][0]:
                        coordinate_map[chain_id][0][res_nr] = [res_name, {}]
                    coordinate_map[chain_id][0][res_nr][1][atom_nr] = (atom_name, x, y, z)

                    if chain_id not in b_factors:
                        b_factors[chain_id] = {}
                    if res_nr not in b_factors[chain_id]:
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

                    if chain_id not in siss_map:
                        siss_map[chain_id] = {}
                    if res_nr not in siss_map[chain_id]:
                        siss_map[chain_id][res_nr] = [res_name, {}]
                    siss_map[chain_id][res_nr][1][atom_nr] = (atom_name, x, y, z)

                    if chain_id not in res_contig_map:
                        res_contig_map[chain_id] = {res_nr: [1, res_name]}
                        contig_help_map[chain_id] = 1
                    elif res_nr not in res_contig_map[chain_id]:
                        contig_help_map[chain_id] += 1
                        res_contig_map[chain_id][res_nr] = [contig_help_map[chain_id], res_name]

            if record_name == "HETATM":
                if len(line) > 50:
                    x = float(line[30:38].replace(" ", ""))
                    y = float(line[38:46].replace(" ", ""))
                    z = float(line[46:54].replace(" ", ""))
                    if chain_id not in coordinate_map:
                        coordinate_map[chain_id] = [{}, {}]
                        box_map[chain_id] = [x, x, y, y, z, z]

                    if ((chain_id, res_nr) in modres_map) or (res_name in sdsc.THREE_TO_ONE) or (res_name in rare_residues):  # If it is a modified residue, than add it to the normal residues...
                        if atom_name[0] in ('H', 'D'):
                            continue
                        if not (chain_id, res_nr) in modres_map:
                            modres_map[(chain_id, res_nr)] = res_name
                        if res_nr not in coordinate_map[chain_id][0]:
                            coordinate_map[chain_id][0][res_nr] = [res_name, {}]
                        coordinate_map[chain_id][0][res_nr][1][atom_nr] = (atom_name, x, y, z)

                        if chain_id not in b_factors:
                            b_factors[chain_id] = {}
                        if res_nr not in b_factors[chain_id]:
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

                        if not res_name in sdsc.THREE_TO_ONE:
                            siss_het_res_name = 'UNK'
                        elif sdsc.THREE_TO_ONE[res_name][0] in sdsc.ONE_TO_THREE:
                            siss_het_res_name = sdsc.ONE_TO_THREE[sdsc.THREE_TO_ONE[res_name][0]]
                        else:
                            siss_het_res_name = 'UNK'
                        if chain_id not in siss_map:
                            siss_map[chain_id] = {}
                        if res_nr not in siss_map[chain_id]:
                            siss_map[chain_id][res_nr] = [siss_het_res_name, {}]
                        siss_map[chain_id][res_nr][1][atom_nr] = (atom_name, x, y, z)

                        if chain_id not in res_contig_map:
                            res_contig_map[chain_id] = {res_nr: [1, res_name]}
                            contig_help_map[chain_id] = 1
                        elif res_nr not in res_contig_map[chain_id]:
                            contig_help_map[chain_id] += 1
                            res_contig_map[chain_id][res_nr] = [contig_help_map[chain_id], res_name]
                    else:
                        if res_nr not in coordinate_map[chain_id][1]:  # If not, then add it to the ligands
                            coordinate_map[chain_id][1][res_nr] = [res_name, {}]
                        coordinate_map[chain_id][1][res_nr][1][atom_nr] = (atom_name, x, y, z)
                        if res_name in sdsc.METAL_ATOMS:
                            metals.add((chain_id, res_nr))
                        elif res_name in sdsc.ION_ATOMS:
                            ions.add((chain_id, res_nr))
                        else:
                            ligands.add((chain_id, res_nr))

    # New Peptide detection counts irregular amino acids
    for chain_id in peptide_count:
        if chain_type_map[chain_id] == 'DNA' or chain_type_map[chain_id] == 'RNA':
            continue
        if peptide_count[chain_id][0] <= peptide_count[chain_id][1] * 2:
            chain_type_map[chain_id] = 'Peptide'
        elif (peptide_count[chain_id][0] + peptide_count[chain_id][1]) < 150:  # Total number of atoms
            chain_type_map[chain_id] = 'Peptide'

    return (coordinate_map, siss_map, res_contig_map, ligands, metals, ions, box_map,
            chain_type_map, chainlist, b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map)

def get_atomic_coverage(coordinate_map, chain):
    expected_number_of_atoms = 0
    resolved_number_of_atoms = 0
    for res_nr in coordinate_map[chain][0]:
        res_name, atoms = coordinate_map[chain][0][res_nr]
        one_letter = THREE_TO_ONE[res_name]
        if one_letter == 'X':
            expected_number_of_atoms += 1
            continue
        if one_letter not in CORRECT_COUNT:
            expected_number_of_atoms += 1
            continue
        expected_number_of_atoms += CORRECT_COUNT[one_letter]
        resolved_number_of_atoms += len(atoms)
    if expected_number_of_atoms == 0:
        return 0
    return (resolved_number_of_atoms / expected_number_of_atoms)

def getMinSubDist(c_map, fuzzy_dist_matrix, target_res_id, res_chain, chain, config, distance_threshold=10.0):
    top20 = []
    top = 20
    for res in c_map[chain][0]:
        if chain == res_chain and res == target_res_id:
            continue
        if (res_chain, chain, target_res_id, res) in fuzzy_dist_matrix:
            d = fuzzy_dist_matrix[(res_chain, chain, target_res_id, res)]
        else:
            continue

        if len(top20) < top:
            top20.append((res, d))
        elif len(top20) == top:
            top20.append((res, d))
            top20 = sorted(top20, key=lambda x: x[1])
        elif d < top20[top][1]:
            top20[top] = (res, d)
            top20 = sorted(top20, key=lambda x: x[1])

    min_sub_d = None
    min_res = None
    min_atom_sub = None
    min_atom_chain = None
    atomlist = c_map[res_chain][0][target_res_id][1]
    inside_sphere = {}
    for (res, d) in top20:
        for atomnr in c_map[chain][0][res][1]:
            # exclude Hydrogens
            if c_map[chain][0][res][1][atomnr][0][0] != 'H':
                coord1 = c_map[chain][0][res][1][atomnr][1:]
                for atomnr2 in atomlist:
                    if atomlist[atomnr2][0][0] != 'H':
                        coord2 = atomlist[atomnr2][1:]
                        d = distance(coord1, coord2)
                        if d < 1.2:
                            continue
                        if d < distance_threshold:
                            if res not in inside_sphere:
                                inside_sphere[res] = d
                            elif d < inside_sphere[res]:
                                inside_sphere[res] = d
                        if min_sub_d is None or d < min_sub_d:
                            min_sub_d = d
                            min_res = res
                            min_atom_sub = atomnr2
                            min_atom_chain = atomnr

    return min_sub_d, min_res, min_atom_sub, min_atom_chain, inside_sphere


def box_check(box_1, box_2, distance_threshold=5.0):
    [min_x_1, max_x_1, min_y_1, max_y_1, min_z_1, max_z_1] = box_1
    [min_x_2, max_x_2, min_y_2, max_y_2, min_z_2, max_z_2] = box_2

    center_1 = [(min_x_1 + max_x_1) / 2., (min_y_1 + max_y_1) / 2., (min_z_1 + max_z_1) / 2.]
    center_2 = [(min_x_2 + max_x_2) / 2., (min_y_2 + max_y_2) / 2., (min_z_2 + max_z_2) / 2.]

    center_dist = distance(center_1, center_2)

    if center_dist < (max(max_x_1 - min_x_1, max_y_1 - min_y_1, max_z_1 - min_z_1) / 2.) + (max(max_x_2 - min_x_2, max_y_2 - min_y_2, max_z_2 - min_z_2) / 2.) + distance_threshold:
        return True, center_dist
    else:
        return False, center_dist


def calcFuzzyDM(coordinate_map, box_map, config, calc_exact_distances=False, target_chains=None):
    # coordinate_map: {Chain:[{Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]},{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}
    fuzzy_dm = {}
    processed_chains = set()
    for chain in coordinate_map:
        for chain_2 in coordinate_map:
            if chain_2 in processed_chains:
                continue
            if target_chains is not None:
                if not (chain in target_chains or chain_2 in target_chains):
                    continue
            neighbors, center_dist = box_check(box_map[chain], box_map[chain_2], distance_threshold=config.short_distance_threshold)
            if not neighbors:
                #fuzzy_dm[(chain,chain_2)] = center_dist
                #fuzzy_dm[(chain_2,chain)] = center_dist
                # print chain,chain_2
                continue

            for res in coordinate_map[chain][0]:
                test_coord = list(coordinate_map[chain][0][res][1].values())[0][1:]
                for res_2 in coordinate_map[chain_2][0]:
                    if chain == chain_2 and res == res_2:
                        continue
                    test_coord_2 = list(coordinate_map[chain_2][0][res_2][1].values())[0][1:]
                    d = distance(test_coord, test_coord_2)
                    if d > 2. * config.short_distance_threshold and d > 2. * config.milieu_threshold:
                        continue
                    if calc_exact_distances:
                        d, atom, atom2 = getMinDist(coordinate_map, res, chain, res_2, chain_2)
                    fuzzy_dm[(chain, chain_2, res, res_2)] = d
                    fuzzy_dm[(chain_2, chain, res_2, res)] = d
        processed_chains.add(chain)
    return fuzzy_dm
# coordinate_map: {Chain:[{Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]},{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}


def getMinDist(c_map, res, chain, res2, chain2, het=0, het2=0):
    min_d = None
    if res not in c_map[chain][het]:
        return None, None, None
    if res2 not in c_map[chain2][het2]:
        return None, None, None

    for atomnr in c_map[chain][het][res][1]:
        coord = c_map[chain][het][res][1][atomnr][1:]
        for atomnr2 in c_map[chain2][het2][res2][1]:
            coord2 = c_map[chain2][het2][res2][1][atomnr2][1:]
            d = distance(coord, coord2)
            if min_d is None or d < min_d:
                min_d = d
                min_atom = atomnr
                min_atom2 = atomnr2

    if min_d is None:
        print(c_map, res, chain, res2, chain2, het, het2)

    return min_d, min_atom, min_atom2

@ray.remote(max_retries = -1, max_calls = 1)
def analysis_chain_remote_wrapper(target_chain, analysis_dump):
    serializedPipeline.ray_hack()

    result = [analysis_chain(target_chain, analysis_dump)]

    #print('Finished nested chain:', target_chain)

    result = pack(result)

    #print('Finished packaging of nested chain:', target_chain)

    return result

@ray.remote(max_retries = -1, max_calls = 1)
def analysis_chain_package_wrapper(package, analysis_dump):
    serializedPipeline.ray_hack()
    results = []
    for target_chain in package:
        results.append(analysis_chain(target_chain, analysis_dump))
    return pack(results)

@ray.remote(max_retries = -1, max_calls = 1)
def annotate(config, chunk):
    serializedPipeline.ray_hack()

    outputs = []
    nested_outputs = []
    for pdb_id, target_dict, nested_cores_limiter, nested_call in chunk:
        t0 = time.time()
        (structural_analysis_dict, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, nested_processes) = structuralAnalysis(pdb_id, config, target_dict=target_dict, nested_cores_limiter = nested_cores_limiter, nested_call = nested_call)
        t1 = time.time()

        if config.verbosity >= 4:
            print('Time for structural analysis of', pdb_id, ':', t1 - t0)

        if nested_call:
            nested_outputs.append((pdb_id, nested_processes, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, t1 - t0))
        else:
            outputs.append((pdb_id, structural_analysis_dict, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, t1 - t0))

    outputs = pack(outputs)

    return (outputs, nested_outputs)

def analysis_chain(target_chain, analysis_dump):
    (pdb_id, config, target_residues, siss_coord_map, centroid_map,
     res_contig_map, coordinate_map, fuzzy_dist_matrix, chain_type_map,
     b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map,
     profiles, milieu_dict, dssp, dssp_dict, page_altered) = analysis_dump
    errorlist = []

    siss_map = {}
    try:
        siss_targets = {}
        siss_targets[target_chain] = list(target_residues[target_chain].keys())
        dist_matrix, res_res_angle_map = spherecon.calcDistMatrix(siss_coord_map, centroid_map, siss_targets, False)

        siss_map[target_chain] = spherecon.calculateSiss(siss_coord_map, centroid_map, dist_matrix, res_res_angle_map, siss_targets, False)[target_chain]

        del dist_matrix
    except:

        errorlist.append("Siss error: %s,%s" % (pdb_id, target_chain))

    structural_analysis_dict = {}

    milieu_dict[target_chain] = {}
    sub_loop_broken = False

    for target_res_id in target_residues[target_chain]:
        if sub_loop_broken:
            break
        milieu_dict[target_chain][target_res_id] = {}
        three_letter = res_contig_map[target_chain][target_res_id][1]
        if three_letter in sdsc.THREE_TO_ONE:
            if sdsc.THREE_TO_ONE[three_letter][0] in sdsc.ONE_TO_THREE:
                one_letter = sdsc.THREE_TO_ONE[three_letter][0]
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

        inter_chain_interactions_median = None
        inter_chain_interactions_dist_weighted = None
        intra_chain_interactions_median = None
        intra_chain_interactions_dist_weighted = None

        inter_chain_kds = []
        inter_chain_dist_weighted_kds = []
        inter_chain_rsas = []
        inter_chain_dist_weighted_rsas = []
        total_inter_dist_weights = 0.

        inter_chain_interactions = []
        inter_chain_interactions_dw = []

        # Distance Calculations
        for chain in coordinate_map:

            #Sub - Chain - Calculations
            if chain != target_chain and len(coordinate_map[chain][0]) > 0:
                min_chain_dist, min_res, atom_sub, atom_chain, inside_sphere = getMinSubDist(coordinate_map, fuzzy_dist_matrix, target_res_id,
                                                                                             target_chain, chain, config, distance_threshold=config.milieu_threshold)
                min_chain_dists[chain] = (min_chain_dist, (atom_chain, atom_sub), min_res)
                if chain_type_map[chain] == 'Protein' and config.neighborhood_calculation:  # Filter DNA and RNA chains

                    for inter_chain_res in inside_sphere:
                        if chain not in dssp_dict:
                            continue
                        if inter_chain_res not in dssp_dict[chain]:
                            continue
                        dist = inside_sphere[inter_chain_res]
                        total_inter_dist_weights += 1 / dist
                        if not inter_chain_res in res_contig_map[chain]:
                            continue
                        if not res_contig_map[chain][inter_chain_res][1] in sdsc.THREE_TO_ONE:
                            continue
                        inter_chain_res_one_letter = sdsc.THREE_TO_ONE[res_contig_map[chain][inter_chain_res][1]][0]
                        kd = sdsc.HYDROPATHY[inter_chain_res_one_letter]
                        inter_chain_kds.append(kd)
                        inter_chain_dist_weighted_kds.append(kd / dist)
                        rsa = dssp_dict[chain][inter_chain_res][0]
                        if rsa is not None:
                            inter_chain_rsas.append(rsa)
                            inter_chain_dist_weighted_rsas.append(rsa / dist)
                        if chain not in profiles:
                            errorlist.append('(1) target_chain not in profiles: %s %s %s' % (pdb_id, chain, list(profiles.keys())))
                            sub_loop_broken = True
                            break
                        else:
                            interface_score = profiles[chain][inter_chain_res][0].getTotalInterfaceInteractionScore()
                            inter_chain_interactions.append(interface_score)
                            inter_chain_interactions_dw.append(interface_score / dist)

            # Residue-Residue Calculations
            elif chain == target_chain and len(coordinate_map[chain][0]) > 0 and config.neighborhood_calculation:
                if chain_type_map[chain] == 'Protein':  # Filter DNA and RNA chains
                    if chain not in dssp_dict:
                        continue
                    min_chain_dist, min_res, atom_sub, atom_chain, inside_sphere = getMinSubDist(coordinate_map, fuzzy_dist_matrix,
                                                                                                 target_res_id, target_chain, chain,
                                                                                                 config, distance_threshold=config.milieu_threshold)
                    intra_chain_kds = []
                    intra_chain_dist_weighted_kds = []
                    intra_chain_rsas = []
                    intra_chain_dist_weighted_rsas = []
                    total_dist_weights = 0.


                    if config.verbosity >= 5:
                        print('Residue-Residue calc:', pdb_id, chain, target_res_id, 'inside the sphere:', len(inside_sphere))

                    for intra_chain_res in inside_sphere:
                        if intra_chain_res not in dssp_dict[chain]:
                            continue
                        if not res_contig_map[chain][intra_chain_res][1] in sdsc.THREE_TO_ONE:
                            continue
                        dist = inside_sphere[intra_chain_res]
                        total_dist_weights += 1 / dist
                        intra_chain_res_one_letter = sdsc.THREE_TO_ONE[res_contig_map[chain][intra_chain_res][1]][0]
                        kd = sdsc.HYDROPATHY[intra_chain_res_one_letter]
                        intra_chain_kds.append(kd)
                        intra_chain_dist_weighted_kds.append(kd / dist)
                        rsa = dssp_dict[chain][intra_chain_res][0]
                        if rsa is not None:
                            intra_chain_rsas.append(rsa)
                            intra_chain_dist_weighted_rsas.append(rsa / dist)
                    if total_dist_weights > 0:
                        intra_chain_median_kd = median(intra_chain_kds)
                        intra_chain_dist_weighted_kd = sum(intra_chain_dist_weighted_kds) / total_dist_weights
                        intra_chain_median_rsa = median(intra_chain_rsas)
                        intra_chain_dist_weighted_rsa = sum(intra_chain_dist_weighted_rsas) / total_dist_weights

                    min_chain_dist, min_res, atom_sub, atom_chain, inside_sphere = getMinSubDist(coordinate_map, fuzzy_dist_matrix,
                                                                                                 target_res_id, target_chain, chain,
                                                                                                 config, distance_threshold=config.intra_milieu_threshold)

                    intra_chain_interactions = []
                    intra_chain_interactions_dw = []
                    total_dist_weights = 0.

                    for intra_chain_res in inside_sphere:
                        dist = inside_sphere[intra_chain_res]
                        total_dist_weights += 1 / dist
                        if not chain in profiles:
                            errorlist.append('(2) target_chain not in profiles: %s %s %s' % (pdb_id, chain, list(profiles.keys())))
                            sub_loop_broken = True
                            break
                        elif intra_chain_res in profiles[chain]:
                            interface_score = profiles[chain][intra_chain_res][0].getTotalInterfaceInteractionScore()
                            intra_chain_interactions.append(interface_score)
                            intra_chain_interactions_dw.append(interface_score / dist)

                    if total_dist_weights > 0:
                        intra_chain_interactions_median = median(intra_chain_interactions)
                        intra_chain_interactions_dist_weighted = sum(intra_chain_interactions_dw) / total_dist_weights

            #Sub - Lig - Calculations
            for hetres in coordinate_map[chain][1]:
                min_d, atom, atom2 = getMinDist(coordinate_map, target_res_id, target_chain, hetres, chain, het2=1)
                # only return ligand distances that are inside a given threshold
                if min_d > config.ligand_interest_sphere:
                    continue
                abr = coordinate_map[chain][1][hetres][0]
                ligand_identifier = "%s_%s_%s" % (abr, hetres, chain)
                lig_dists[ligand_identifier] = (min_d, (atom, atom2))

        if total_inter_dist_weights > 0:
            inter_chain_median_kd = median(inter_chain_kds)
            inter_chain_dist_weighted_kd = sum(inter_chain_dist_weighted_kds) / total_inter_dist_weights
            inter_chain_median_rsa = median(inter_chain_rsas)
            inter_chain_dist_weighted_rsa = sum(inter_chain_dist_weighted_rsas) / total_inter_dist_weights
            inter_chain_interactions_median = median(inter_chain_interactions)
            inter_chain_interactions_dist_weighted = sum(inter_chain_interactions_dw) / total_inter_dist_weights

        # coordinate_map: {Chain:[{Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]},{Hetatm-Residue (inlcuding insertion_code):[Resname,{atomnr:(atomname,x,y,z)}]}]}

        # Homomer Sub - Sub Calculations:
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
                    (rsa, relative_main_chain_acc, relative_side_chain_acc, ssa, phi, psi) = dssp_dict[target_chain][target_res_id]
                else:
                    # print(target_res_id)
                    # print(siss_map)
                    siss_value = None
                    if target_res_id in siss_map:
                        siss_value = siss_map[target_chain][target_res_id]

                    ssa = None
                    rsa = siss_value
                    relative_main_chain_acc = None
                    relative_side_chain_acc = None
                    phi = None
                    psi = None
            else:
                atomic_coverage = get_atomic_coverage(coordinate_map, chain)
                if atomic_coverage > 0.5 and not page_altered: #DSSP does not work for chains with missing atomic coordinates or chains over the 10K atoms mark, no need for a warning here
                    errorlist.append(("dssp error: chain not in dssp_dict; %s; %s" % (pdb_id, target_chain)))
                dssp = False
                if target_res_id in siss_map:
                    siss_value = siss_map[target_chain][target_res_id]
                else:
                    siss_value = None
                ssa = None
                rsa = siss_value
                relative_main_chain_acc = None
                relative_side_chain_acc = None
                phi = None
                psi = None
        else:
            siss_value = None
            if target_res_id in siss_map:
                siss_value = siss_map[target_chain][target_res_id]

            ssa = None
            phi = None
            psi = None
            relative_main_chain_acc = None
            relative_side_chain_acc = None
            rsa = siss_value
        if target_chain in profiles:
            if target_res_id in profiles[target_chain]:
                profile, centrality_scores = profiles[target_chain][target_res_id]

                #if profile is not None:
                #    profile = profile.encode()
                #if centrality_scores is not None:
                #    centrality_scores = centrality_scores.str_encode()

                if profile is None:
                    errorlist.append('profile is None: %s %s %s' % (pdb_id, target_chain, target_res_id))
            else:
                #errorlist.append('target_res not in profiles')
                profile = None
                centrality_scores = None
        else:
            #errorlist.append('target_chain not in profiles: %s %s' % (target_chain,list(profiles.keys())))
            profile = None
            centrality_scores = None

        avg_b_factor = sum(b_factors[target_chain][target_res_id]) / float(len(b_factors[target_chain][target_res_id]))

        if (target_chain, target_res_id) in modres_map:
            modres = True
        else:
            modres = False

        if (target_chain, target_res_id) in ssbond_map:
            (chain_2, res_nr_2, ssbond_len) = ssbond_map[(target_chain, target_res_id)]
            intra_ssbond = target_chain == chain_2
        else:
            intra_ssbond = None
            ssbond_len = None

        if (target_chain, target_res_id) in link_map:
            (atom_1, res_name_1, atom_2, res_name_2, res_nr_2, chain_2, link_dist) = link_map[(target_chain, target_res_id)]
            intra_link = target_chain == chain_2
        else:
            intra_link = None
            link_dist = None

        if (target_chain, target_res_id) in cis_conformation_map:
            (res_name_2, chain_2, res_nr_2, angle) = cis_conformation_map[(target_chain, target_res_id)]
            cis_conformation = angle
        else:
            cis_conformation = None

        if (target_chain, target_res_id) in cis_follower_map:
            (res_name_2, chain_2, res_nr_2, angle) = cis_follower_map[(target_chain, target_res_id)]
            cis_follower = angle
        else:
            cis_follower = None

        residue = sdsc.Residue(target_res_id, aa=one_letter, lig_dists=lig_dists, chain_distances=min_chain_dists, RSA=rsa,
                           relative_main_chain_acc=relative_main_chain_acc, relative_side_chain_acc=relative_side_chain_acc,
                           SSA=ssa, homomer_distances=homomer_map, interaction_profile=profile, centralities=centrality_scores,
                           modres=modres, b_factor=avg_b_factor, phi=phi, psi=psi, intra_ssbond=intra_ssbond, ssbond_length=ssbond_len,
                           intra_link=intra_link, link_length=link_dist, cis_conformation=cis_conformation, cis_follower=cis_follower,
                           inter_chain_median_kd=inter_chain_median_kd, inter_chain_dist_weighted_kd=inter_chain_dist_weighted_kd,
                           inter_chain_median_rsa=inter_chain_median_rsa, inter_chain_dist_weighted_rsa=inter_chain_dist_weighted_rsa,
                           intra_chain_median_kd=intra_chain_median_kd, intra_chain_dist_weighted_kd=intra_chain_dist_weighted_kd,
                           intra_chain_median_rsa=intra_chain_median_rsa, intra_chain_dist_weighted_rsa=intra_chain_dist_weighted_rsa,
                           inter_chain_interactions_median=inter_chain_interactions_median,
                           inter_chain_interactions_dist_weighted=inter_chain_interactions_dist_weighted,
                           intra_chain_interactions_median=intra_chain_interactions_median,
                           intra_chain_interactions_dist_weighted=intra_chain_interactions_dist_weighted)
        residue.convert_centrality_str()
        residue.convert_interaction_profile_str()
        residue.generate_interacting_chains_str()
        residue.generate_interacting_ligands_str()
        structural_analysis_dict[target_res_id] = residue

    return target_chain, structural_analysis_dict, errorlist


# called by serializedPipeline
# @profile
def structuralAnalysis(pdb_id, config, target_dict=None, model_path=None, keep_rin_files=True, nested_cores_limiter = None, nested_call = False):

    dssp = config.dssp
    dssp_path = config.dssp_path
    pdb_path = config.pdb_path
    rin_db_path = config.rin_db_path
    verbosity = config.verbosity

    if verbosity >= 4:
        t0 = time.time()
        print('Start structuralAnalysis of:', pdb_id)

    page, page_altered, path, atom_count = pdbParser.standardParsePDB(pdb_id, pdb_path, return_10k_bool=True, get_is_local=True, model_path=model_path)

    if verbosity >= 4:
        t1 = time.time()
        print('Time for structuralAnalysis Part 1:', t1 - t0)

    if verbosity >= 5:
        print('\n',pdb_id,'\n\n',page)

    errorlist = []

    if page == '':
        errorlist.append("Error while parsing: %s" % pdb_id)
        return {}, errorlist, {}, {}, {}, {}, {}, []

    (coordinate_map, siss_coord_map, res_contig_map, ligands, metals, ions, box_map, chain_type_map, chainlist,
        b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map) = parsePDB(page)

    if verbosity >= 4:
        t2 = time.time()
        print('Time for structuralAnalysis Part 2:', t2 - t1)

    if target_dict is None:
        target_residues = {}
        for chain in siss_coord_map:
            target_residues[chain] = list(siss_coord_map[chain].keys())
    elif isinstance(target_dict, list):
        target_residues = {}
        for chain in target_dict:
            target_residues[chain] = list(siss_coord_map[chain].keys())
    else:
        target_residues = target_dict

    centroid_map = spherecon.calcCentroidMap(siss_coord_map, target_residues, False)

    if verbosity >= 4:
        t3 = time.time()
        print('Time for structuralAnalysis Part 3:', t3 - t2)

    fuzzy_dist_matrix = calcFuzzyDM(coordinate_map, box_map, config, target_chains=set(target_residues.keys()))

    if verbosity >= 4:
        t4 = time.time()
        print('Time for structuralAnalysis Part 4:', t4 - t3)

    structural_analysis_dict = {}

    milieu_dict = {}

    if dssp:
        # write temp pdb file only if we had to fix the 10k atom bug or the file is not stored locally
        if page_altered or path is None:
            tmp_path = '%s/tmp_%s.pdb' % (config.temp_folder, pdb_id)
            f = open(tmp_path, 'w')
            f.write(page)
            f.close()
        else:
            tmp_path = path

        # call DSSP  -structure of dssp_dict: {Chain:{Res_id:(racc,ssa)}}
        dssp_dict, errlist = calcDSSP(tmp_path, dssp_path, angles=True, verbosity_level=verbosity)
        errorlist += errlist
        # remove tmp_file
        if page_altered:
            try:
                os.remove(tmp_path)
            except:
                pass
        if dssp_dict == {}:
            dssp = False

    if verbosity >= 4:
        t5 = time.time()
        print('Time for structuralAnalysis Part 5:', t5 - t4)

    profiles = {}
    ligand_profiles = {}
    metal_profiles = {}
    ion_profiles = {}
    chain_chain_profiles = {}

    try:
        lookup_out = rin.lookup(pdb_id, page, config, None, None, ligands,
                                metals, ions, res_contig_map,
                                rin_db_path, chain_type_map, encoded=False, model_path=model_path, keep_tmp_files=keep_rin_files)
        if isinstance(lookup_out, str):
            errorlist.append(lookup_out)
        else:
            (profiles, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles) = lookup_out
            if len(profiles) == 0:
                errorlist.append('Empty profiles %s %s %s' % (pdb_id, str(len(page)), str(len(res_contig_map))))
            elif verbosity >= 4:
                print('Chains in profiles:', list(profiles.keys()))
    except:
        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        errortext = 'RIN lookup Error:%s\nRIN_db_path:%s\n' % (pdb_id, rin_db_path) + '\n'.join([str(e), str(f), str(g)])
        errorlist.append(errortext)

    if verbosity >= 4:
        t6 = time.time()
        print('Time for structuralAnalysis Part 6:', t6 - t5)

    if target_dict is None:
        target_residues = res_contig_map
    elif isinstance(target_dict, list):
        target_residues = {}
        for chain in target_dict:
            target_residues[chain] = res_contig_map[chain]
    else:
        target_residues = target_dict

    #parent_dir = '/wibicom/SHARED_DATA/agress/structman/lib'
    #os.environ["PYTHONPATH"] = parent_dir + ":" + os.environ.get("PYTHONPATH", "")

    #parent_dir = '/wibicom/SHARED_DATA/agress/structman'
    #os.environ["PYTHONPATH"] = parent_dir + ":" + os.environ.get("PYTHONPATH", "")

    number_protein_chains = 0
    for target_chain in target_residues:
        if chain_type_map[target_chain] != 'Protein':
            continue
        number_protein_chains += 1

    if nested_cores_limiter is None:
        nested_cores_limiter = config.proc_n - 1

    if nested_call and not config.low_mem_system:

        analysis_dump = ray.put((pdb_id, config, target_residues, siss_coord_map, centroid_map,
                                 res_contig_map, coordinate_map, fuzzy_dist_matrix, chain_type_map,
                                 b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map,
                                 profiles, milieu_dict, dssp, dssp_dict, page_altered))

        core_results = []

        nested_cores = min([config.proc_n - 1, nested_cores_limiter])

        if number_protein_chains < nested_cores:
            packaged = False
            for target_chain in target_residues:
                if chain_type_map[target_chain] != 'Protein':
                    continue

                core_results.append(analysis_chain_remote_wrapper.remote(target_chain, analysis_dump))

        else:
            if verbosity >= 2:
                print('Packaging a big structure:',pdb_id,'with',number_protein_chains,'chains')
            packaged = True
            max_package_size = number_protein_chains // nested_cores
            if (number_protein_chains % nested_cores) != 0:
                max_package_size += 1

            package = []
            for target_chain in target_residues:
                if chain_type_map[target_chain] != 'Protein':
                    continue
                package.append(target_chain)
                if len(package) == max_package_size:
                    core_results.append(analysis_chain_package_wrapper.remote(package, analysis_dump))
                    package = []

            if len(package) != 0:
                core_results.append(analysis_chain_package_wrapper.remote(package, analysis_dump))

        if verbosity >= 3:
            print('Reached nested remotes:', pdb_id, number_protein_chains)

        return {}, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, core_results

    else:
        analysis_dump = (pdb_id, config, target_residues, siss_coord_map, centroid_map,
                         res_contig_map, coordinate_map, fuzzy_dist_matrix, chain_type_map,
                         b_factors, modres_map, ssbond_map, link_map, cis_conformation_map, cis_follower_map,
                         profiles, milieu_dict, dssp, dssp_dict, page_altered)
        for target_chain in target_residues:
            if chain_type_map[target_chain] != 'Protein':
                continue

            target_chain, chain_structural_analysis_dict, chain_errorlist = analysis_chain(target_chain, analysis_dump)
            structural_analysis_dict[target_chain] = chain_structural_analysis_dict
            errorlist += ['%s - %s' % (x, pdb_id) for x in chain_errorlist[:5]]

    if verbosity >= 4:
        t7 = time.time()
        print('Time for structuralAnalysis Part 7:', t7 - t6)

    return structural_analysis_dict, errorlist, ligand_profiles, metal_profiles, ion_profiles, chain_chain_profiles, chain_type_map, chainlist, []

def paraAnnotate(config, proteins, lite=False, indel_analysis_follow_up=False):

    annotation_processes = config.proc_n

    t0 = time.time()

    # new structure for paralell annotation: give each annotation process a pdb and dict of chain and structure information,
    # all chains without structure information shall go into the database as fully annotated structures, which are not aligned (yet)
    # these are unmapped chains, or protein interaction partners, their presence in the database is necesary for the computation of the structural neighborhood
    # a new structure-structure entry has to be created and has to be inserted into the database before the residue insertion
    # homooligomer information is not present for such cases, this information has to be updated, if the chain is mapped at a later iteration or run

    size_map = {}

    total_structural_analysis = {}
    interaction_structures = set()
    background_insert_residues_process = None
    structure_list = proteins.get_structure_list()

    complex_list = proteins.get_complex_list()

    n_of_stored_complexes = 0
    n_of_comps = 0
    n_of_chains_to_analyze = 0
    n_of_small_comps = 0

    for pdb_id in complex_list:
        if proteins.is_complex_stored(pdb_id):
            n_of_stored_complexes += 1
            continue
        s = len(proteins.get_complex_chains(pdb_id, only_protein=True))
        ac = proteins.get_atom_count(pdb_id)
        if ac is None:
            config.errorlog.add_warning('Complex with None as atom count: %s' % pdb_id)
        if s not in size_map:
            size_map[s] = {}
        size_map[s][pdb_id] = ac
        n_of_comps += 1
        n_of_chains_to_analyze += s
        if s < config.n_of_chain_thresh:
            n_of_small_comps += 1

    chunksize = min([8,max([n_of_small_comps // (4 * config.proc_n), 1])])

    if config.verbosity >= 2:
        print('Starting structural analysis with', n_of_comps, 'complexes.', n_of_stored_complexes, ' complexes are already stored.')
        print('Total amount in structure_list:', len(structure_list), '. Amount of structures to analyze:', n_of_chains_to_analyze)

    sorted_sizes = sorted(size_map.keys(), reverse=True)

    anno_result_ids = []

    package_cost = 0
    assigned_costs = {}
    temporarily_freed = {}
    conf_dump = ray.put(config)
    cost_function_constant = 1280000000

    n_started = 0

    current_chunk = []

    nested_packages = set()

    for s in sorted_sizes:
        if s < config.n_of_chain_thresh:
            cost = 1
        else:
            cost = min([s, config.proc_n])
        '''
        if config.low_mem_system and s >= 15: #Skip large structure for low mem systems
            del size_map[s]
            continue
        '''
        del_list = []
        for pdb_id in size_map[s]:
            ac = size_map[s][pdb_id]
            if config.low_mem_system:
                cost = max([1, min([config.proc_n, (((ac**2) / cost_function_constant) * config.proc_n) // config.gigs_of_ram])])

            #if (cost + package_cost) <= config.proc_n:
            if (1 + package_cost) <= config.proc_n:
                if not lite:
                    target_dict = None
                else:
                    target_dict = proteins.get_target_dict(pdb_id)

                if s < config.n_of_chain_thresh:
                    current_chunk.append((pdb_id, target_dict, None, False))
                    if len(current_chunk) >= chunksize:
                        anno_result_ids.append(annotate.remote(conf_dump, current_chunk))
                        current_chunk = []
                        package_cost += 1  # chunks cost 1
                else:
                    anno_result_ids.append(annotate.remote(conf_dump, [(pdb_id, target_dict, None, True)]))
                    package_cost += cost  # big structures with nested remotes cost the amount of nested calls
                    nested_packages.add(pdb_id)
                    if config.verbosity >= 3:
                        print('started nested package:', pdb_id, cost)
                n_started += 1
                del_list.append(pdb_id)
                assigned_costs[pdb_id] = cost
                temporarily_freed[pdb_id] = 0
            else:
                break

        for pdb_id in del_list:
            del size_map[s][pdb_id]
        if len(size_map[s]) == 0:
            del size_map[s]

    sorted_sizes = sorted(size_map.keys(), reverse=True)

    if config.verbosity >= 2:
        t11 = time.time()
        print("Annotation Part 1.1: %s" % (str(t11 - t0)))

    database.getStoredResidues(proteins, config)  # has to be called after insertStructures
    #anno_results = ray.get(anno_result_ids)

    if config.verbosity >= 2:
        t1 = time.time()
        print("Annotation Part 1.2: %s" % (str(t1 - t11)))
        print("Annotation Part 1: %s" % (str(t1 - t0)))

    max_comp_time = 0.
    total_comp_time = 0.
    max_comp_time_structure = None
    amount_of_structures = 0
    amount_of_chains_in_analysis_dict = 0

    total_computed = 0

    interacting_structure_ids = {}
    no_result_run = 1

    core_not_ready_dict = {}

    loop_counter = 0

    #gc.disable()

    while True:
        loop_counter += 1
        chunksize = min([8,max([(n_of_small_comps - n_started) // (4 * config.proc_n), 1])])  # dynamically adjust chunksize
        finished = 0
        del_list = []
        # this loop collects the results from the nested remotes
        for ret_pdb_id in core_not_ready_dict:
            while True:
                core_ready, core_not_ready = ray.wait(core_not_ready_dict[ret_pdb_id], timeout=0.1)
                if config.verbosity >= 3:
                    print('Nested wait:', ret_pdb_id, len(core_ready), len(core_not_ready))
                if len(core_ready) > 0:
                    errorlist = []
                    structural_analysis_dict = {}
                    core_outs = ray.get(core_ready)

                    for core_package in core_outs:
                        core_package = unpack(core_package)
                        for target_chain, chain_structural_analysis_dict, chain_errorlist in core_package:
                            structural_analysis_dict[target_chain] = chain_structural_analysis_dict
                            errorlist += chain_errorlist
                            package_cost -= 1  # relieve cost for each finished nested remote (individual chains)
                            temporarily_freed[ret_pdb_id] += 1
                    if len(errorlist) > 0:
                        for error_text in errorlist:
                            config.errorlog.add_warning(error_text)
                    for chain in structural_analysis_dict:
                        if lite:
                            for res_id in structural_analysis_dict[chain]:
                                residue = structural_analysis_dict[chain][res_id]
                                proteins.add_residue(ret_pdb_id, chain, res_id, residue)
                        else:
                            total_structural_analysis[(ret_pdb_id, chain)] = structural_analysis_dict[chain]
                            amount_of_chains_in_analysis_dict += 1

                            if not (ret_pdb_id, chain) in structure_list:
                                interaction_structures.add((ret_pdb_id, chain))
                    #for core_ray_id in core_ready:
                    #    ray.cancel(core_ray_id, force=True)
                if len(core_not_ready) == 0:
                    del_list.append(ret_pdb_id)
                    finished += 1
                    break
                else:
                    core_not_ready_dict[ret_pdb_id] = core_not_ready
                if len(core_ready) == 0:
                    break


        for ret_pdb_id in del_list:
            del core_not_ready_dict[ret_pdb_id]
            nested_packages.remove(ret_pdb_id)
            package_cost = (package_cost - assigned_costs[ret_pdb_id]) + temporarily_freed[ret_pdb_id] # correct package cost, when the number of protein chains got estimated incorrectly
        
        not_ready = []
        if len(anno_result_ids) > 0:
            if config.verbosity >= 3:
                print('Before wait, loop counter:', loop_counter)

            '''
            ready, not_ready = ray.wait(anno_result_ids, timeout = 0.00001)
            if len(ready) == 0:
                t_gc_0 = time.time()
                gc.collect()
                t_gc_1 = time.time()
                if config.verbosity >= 3:
                    print('Manual gc happend:', (t_gc_1 - t_gc_0))
            '''
            ready, not_ready = ray.wait(anno_result_ids)

            if config.verbosity >= 3:
                print('Wait happened:', no_result_run)

            t_0 = time.time()

            if len(ready) > 0:
                chunk_struct_outs = ray.get(ready)
                no_result_run = 1
            else:
                chunk_struct_outs = []
                no_result_run = min([10, no_result_run + 1])

            t_1 = time.time()
            if config.verbosity >= 3:
                print('Time for loop part 1:', (t_1-t_0))

            t_integrate_0 = 0.
            t_integrate_1 = 0.
            t_integrate_2 = 0.
            t_integrate_3 = 0.
            t_integrate_4 = 0.
            t_pickle_0 = 0.
            t_pickle_1 = 0.
            t_pickle_2 = 0.

            returned_pdbs = []

            for chunk_struct_out, nested_struct_outs in chunk_struct_outs:
                package_cost -= 1  # finished chunks free 1 cost
                t_pickle_0 += time.time()
                chunk_struct_out = unpack(chunk_struct_out)
                t_pickle_1 += time.time()
                for struct_out in chunk_struct_out:
                    t_integrate_0 += time.time()
                    (ret_pdb_id, structural_analysis_dict, errorlist, ligand_profiles, metal_profiles,
                        ion_profiles, chain_chain_profiles, chain_type_map, chainlist, comp_time) = struct_out

                    returned_pdbs.append(ret_pdb_id)

                    finished += 1

                    if comp_time > max_comp_time:
                        max_comp_time = comp_time
                        max_comp_time_structure = ret_pdb_id
                    amount_of_structures += 1
                    total_comp_time += comp_time

                    #This part adds the results from the structural analysis into the central datastructure
                    t_integrate_1 += time.time()
                    if len(errorlist) > 0:
                        for error_text in errorlist:
                            config.errorlog.add_warning(error_text)

                    for chain in structural_analysis_dict:
                        if lite:
                            for res_id in structural_analysis_dict[chain]:
                                residue = structural_analysis_dict[chain][res_id]
                                proteins.add_residue(ret_pdb_id, chain, res_id, residue)
                        else:
                            total_structural_analysis[(ret_pdb_id, chain)] = structural_analysis_dict[chain]
                            amount_of_chains_in_analysis_dict += 1

                            if not (ret_pdb_id, chain) in structure_list:
                                interaction_structures.add((ret_pdb_id, chain))

                    t_integrate_2 += time.time()

                    proteins.set_lig_profile(ret_pdb_id, ligand_profiles)
                    proteins.set_ion_profile(ret_pdb_id, ion_profiles)
                    proteins.set_metal_profile(ret_pdb_id, metal_profiles)
                    proteins.set_chain_chain_profile(ret_pdb_id, chain_chain_profiles)

                    t_integrate_3 += time.time()

                    proteins.set_chain_type_map(ret_pdb_id, chain_type_map, chainlist)

                    if config.low_mem_system:
                        if amount_of_chains_in_analysis_dict > 2 * (config.chunksize):
                            if background_insert_residues_process is not None:
                                background_insert_residues_process.join()
                            interacting_structure_ids = database.insertInteractingChains(interaction_structures, proteins, config)
                            interaction_structures = set()

                            background_insert_residues_process = database.insertResidues(total_structural_analysis, interacting_structure_ids, proteins, config)
                            total_structural_analysis = {}
                            amount_of_chains_in_analysis_dict = 0
                            gc.collect()
                    t_integrate_4 += time.time()
                #print('Amount of nested outs:',len(nested_struct_outs))
                for struct_out in nested_struct_outs:
                    package_cost += 1 #repay one cost, since this a nested out
                    t_integrate_0 += time.time()
                    (ret_pdb_id, nested_processes, errorlist, ligand_profiles, metal_profiles,
                        ion_profiles, chain_chain_profiles, chain_type_map, chainlist, comp_time) = struct_out

                    returned_pdbs.append(ret_pdb_id)

                    core_ready, core_not_ready = ray.wait(nested_processes, timeout=0.01)
                    structural_analysis_dict = {}
                    if len(core_ready) > 0:
                        core_outs = ray.get(core_ready)

                        for core_package in core_outs:
                            core_package = unpack(core_package)
                            for target_chain, chain_structural_analysis_dict, chain_errorlist in core_package:
                                structural_analysis_dict[target_chain] = chain_structural_analysis_dict
                                errorlist += chain_errorlist
                                package_cost -= 1  # relieve cost for each finished nested remote (individual chains)
                                temporarily_freed[ret_pdb_id] += 1

                    if len(core_not_ready) > 0:
                        core_not_ready_dict[ret_pdb_id] = core_not_ready
                    else:
                        finished += 1
                        if ret_pdb_id in core_not_ready_dict:
                            del core_not_ready_dict[ret_pdb_id]
                        nested_packages.remove(ret_pdb_id)
                        package_cost = (package_cost - assigned_costs[ret_pdb_id]) + temporarily_freed[ret_pdb_id] # correct package cost, when the number of protein chains got estimated incorrectly

                    if comp_time > max_comp_time:
                        max_comp_time = comp_time
                        max_comp_time_structure = ret_pdb_id
                    amount_of_structures += 1
                    total_comp_time += comp_time

                    #This part adds the results from the structural analysis into the central datastructure
                    t_integrate_1 += time.time()
                    if len(errorlist) > 0:
                        for error_text in errorlist:
                            config.errorlog.add_warning(error_text)

                    for chain in structural_analysis_dict:
                        if lite:
                            for res_id in structural_analysis_dict[chain]:
                                residue = structural_analysis_dict[chain][res_id]
                                proteins.add_residue(ret_pdb_id, chain, res_id, residue)
                        else:
                            total_structural_analysis[(ret_pdb_id, chain)] = structural_analysis_dict[chain]
                            amount_of_chains_in_analysis_dict += 1

                            if not (ret_pdb_id, chain) in structure_list:
                                interaction_structures.add((ret_pdb_id, chain))

                    t_integrate_2 += time.time()

                    proteins.set_lig_profile(ret_pdb_id, ligand_profiles)
                    proteins.set_ion_profile(ret_pdb_id, ion_profiles)
                    proteins.set_metal_profile(ret_pdb_id, metal_profiles)
                    proteins.set_chain_chain_profile(ret_pdb_id, chain_chain_profiles)

                    t_integrate_3 += time.time()

                    proteins.set_chain_type_map(ret_pdb_id, chain_type_map, chainlist)

                    if config.low_mem_system:
                        if amount_of_chains_in_analysis_dict > 2 * (config.chunksize):
                            if background_insert_residues_process is not None:
                                background_insert_residues_process.join()
                            interacting_structure_ids = database.insertInteractingChains(interaction_structures, proteins, config)
                            interaction_structures = set()

                            background_insert_residues_process = database.insertResidues(total_structural_analysis, interacting_structure_ids, proteins, config)
                            total_structural_analysis = {}
                            amount_of_chains_in_analysis_dict = 0
                            gc.collect()
                    t_integrate_4 += time.time()
            t_2 = time.time()

            if config.verbosity >= 3:
                print('Time for result integration:', (t_integrate_4 - t_integrate_0), 'for amount of analyzed chunks:', len(chunk_struct_outs),
                      'individual times:', (t_integrate_1 - t_integrate_0),
                                           (t_integrate_2 - t_integrate_1),
                                           (t_integrate_3 - t_integrate_2),
                                           (t_integrate_4 - t_integrate_3))
                print('Time for unpacking:', (t_pickle_1 - t_pickle_0))

                if (t_pickle_2 - t_pickle_1) > 60:
                    print('Long depickling:', returned_pdbs)

            if config.verbosity >= 3:
                print('Time for loop part 2:', (t_2-t_1))

        new_anno_result_ids = []
        if len(size_map) > 0:
            # Start new jobs regarding the freed resources

            for s in sorted_sizes:
                if s < config.n_of_chain_thresh:
                    cost = 1
                else:
                    cost = min([s, config.proc_n])

                del_list = []
                for pdb_id in size_map[s]:
                    ac = size_map[s][pdb_id]
                    if config.low_mem_system:
                        cost = max([1, min([config.proc_n, (((ac**2) / cost_function_constant) * config.proc_n) // config.gigs_of_ram])])
                    #if (cost + package_cost) <= config.proc_n:
                    if (1 + package_cost) <= config.proc_n:
                        if not lite:
                            target_dict = None
                        else:
                            target_dict = proteins.get_target_dict(pdb_id)

                        if s < config.n_of_chain_thresh:
                            current_chunk.append((pdb_id, target_dict, None, False))
                            if len(current_chunk) >= chunksize:
                                new_anno_result_ids.append(annotate.remote(conf_dump, current_chunk))
                                current_chunk = []
                                package_cost += 1
                        else:
                            new_anno_result_ids.append(annotate.remote(conf_dump, [(pdb_id, target_dict, None, True)]))
                            package_cost += cost
                            nested_packages.add(pdb_id)
                            if config.verbosity >= 3:
                                print('started nested package:', pdb_id, cost)
                        n_started += 1

                        del_list.append(pdb_id)
                        assigned_costs[pdb_id] = cost
                        temporarily_freed[pdb_id] = 0
                    elif s >= config.n_of_chain_thresh and package_cost <= ((3*config.proc_n)//4): #if there are some free resource, use them
                        cost = min([s, config.proc_n - package_cost])
                        new_anno_result_ids.append(annotate.remote(conf_dump, [(pdb_id, target_dict, cost, True)]))
                        package_cost += cost
                        nested_packages.add(pdb_id)

                        n_started += 1

                        del_list.append(pdb_id)
                        assigned_costs[pdb_id] = cost
                        temporarily_freed[pdb_id] = 0

                        if config.verbosity >= 3:
                            print('started nested package:', pdb_id, cost)
                    else:
                        break

                for pdb_id in del_list:
                    del size_map[s][pdb_id]
                if len(size_map[s]) == 0:
                    del size_map[s]
            sorted_sizes = sorted(size_map.keys(), reverse=True)
        elif len(current_chunk) > 0:  # starting the last chunk
            new_anno_result_ids.append(annotate.remote(conf_dump, current_chunk))
            current_chunk = []
            package_cost += 1
            n_started += 1

        if len(anno_result_ids) > 0 and config.verbosity >= 3:
            t_3 = time.time()
            print('Time for loop part 3:', (t_3-t_2))

        #for ray_id in ready:
        #    ray.cancel(ray_id, force=True)

        anno_result_ids = new_anno_result_ids + not_ready

        if finished > 0:
            total_computed += finished
            if config.verbosity >= 3:
                print('Newly finished structural analysis of', finished, '. Total:', total_computed, 'len of size map:', len(size_map), 'Current package:', len(anno_result_ids), 'with cost:', package_cost, 'and nested packages:', str(nested_packages), ', number of pending nested packaged chains:', len(core_not_ready_dict), 'Amount of new started processes:', len(new_anno_result_ids), 'Amount of pending processes:', len(not_ready))
        if len(anno_result_ids) == 0 and len(size_map) == 0 and len(core_not_ready_dict) == 0:
            break

    t_1_5 = time.time()
    #gc.enable()

    if config.verbosity >= 2:
        t2 = time.time()
        print('Time for reenabling gc:', (t2- t_1_5))
        print("Annotation Part 2: %s" % (str(t2 - t1)))
        print('Longest computation for:', max_comp_time_structure, 'with:', max_comp_time, 'In total', amount_of_structures, 'structures', 'Accumulated time:', total_comp_time)

    if background_insert_residues_process is not None:
        background_insert_residues_process.join()

    if not lite:
        interacting_structure_ids = database.insertInteractingChains(interaction_structures, proteins, config)

        if config.verbosity >= 2:
            t32 = time.time()
            print('Annotation Part 3.1:', t32 - t2)

        database.insertComplexes(proteins, config)

        if config.verbosity >= 2:
            t33 = time.time()
            print('Annotation Part 3.2:', t33 - t32)

        if len(total_structural_analysis) > 0:
            background_insert_residues_process = database.insertResidues(total_structural_analysis, interacting_structure_ids, proteins, config)

        if config.verbosity >= 2:
            t34 = time.time()
            print('Annotation Part 3.3:', t34 - t33)
    else:
        background_insert_residues_process = None

    if config.verbosity >= 2:
        t3 = time.time()
        print("Annotation Part 3: %s" % (str(t3 - t2)))

    serializedPipeline.classification(proteins, config, indel_analysis_follow_up=indel_analysis_follow_up)

    if config.verbosity >= 2:
        t4 = time.time()
        print("Annotation Part 4: %s" % (str(t4 - t3)))

    if not lite:
        database.insertClassifications(proteins, config)

    if config.verbosity >= 2:
        t5 = time.time()
        print("Annotation Part 5: %s" % (str(t5 - t4)))

    gc.collect()

    if config.verbosity >= 2:
        t6 = time.time()
        print("Time for garbage collection: %s" % (str(t6 - t5)))

    return background_insert_residues_process, amount_of_structures


# Return True, if the template is not good
def weakCriteria(seq, res, rel_aln_len, seq_thresh, res_thresh, rel_aln_len_thresh):

    if float(seq) <= 1.0:
        seq = float(seq) * 100.0
    if float(seq) < seq_thresh:
        return True
    if float(res) > res_thresh:
        return True
    else:
        return False


# template-structure: [pdb_id,seq_id,chain,aln_length,resolution,ligands,r-value,templateSelectionScore]


def good(structure, seq_threshold, resolution_threshold, cov_threshold):
    seq_id = structure['Seq_Id']
    resolution = structure['Resolution']
    cov = structure['Coverage']
    if weakCriteria(seq_id, resolution, cov, seq_threshold, resolution_threshold, cov_threshold):
        return False
    return True


# called by serializedPipeline
def filterTemplates(structures, seq_threshold, resolution_threshold, cov_threshold):
    # print structures
    filtered_structures = {}
    for (pdb_id, chain) in structures:
        if good(structures[(pdb_id, chain)], seq_threshold, resolution_threshold, cov_threshold):
            filtered_structures[(pdb_id, chain)] = structures[(pdb_id, chain)]

    return filtered_structures
