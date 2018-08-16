import os
import sys
import getopt
import traceback
import math
import numpy
import time

vdw_radius = {"C":1.7,"O":1.52,"N":1.55,"F":1.47,"P":1.8,"S":1.8,"X":1.7}
#threshs = [3.5,3.75,4.0,4.25,4.5,4.75,5.0,5.25,5.5,5.75,6.0,6.25,6.5,6.75,7.0]
#threshs = [3.5,3.75,4.0,4.25,4.5,4.75,5.0,5.25,5.5,5.75,6.0,6.25,6.5,6.75,7.0,7.25,7.5,7.75,8.0,8.25,8.5,8.75,
#9.0,9.25,9.5,9.75,10.0,10.25,10.5,10.75,11.0,11.25,11.75,12.0,12.25,12.5,12.75,
#13.0,13.25,13.4,13.75,14.0,14.25,14.5,14.75,15.0,15.25,15.0,15.25,15.5,15.75,16.0,16.25,16.5,16.75]

cone_cut_volumes = {'CYS': 542.8672105403161, 'ILE': 859.5397500221673, 'SER': 592.3734747731354, 'GLN': 859.5397500221673, 'LYS': 662.0644718052689, 'ASN': 662.0644718052689, 'PRO': 662.0644718052689, 'THR': 592.3734747731354, 'PHE': 1092.8291844899893, 'ALA': 487.836979224935, 'HIS': 696.9099703213358, 'GLY': 471.23889803846896, 'ASP': 662.0644718052689, 'SEC': 542.8672105403161, 'LEU': 814.3008158104743, 'ARG': 814.3008158104743, 'TRP': 1364.9172882296452, 'VAL': 814.3008158104743, 'GLU': 662.0644718052689, 'TYR': 1092.8291844899893, 'MET': 859.5397500221673}

threshs = numpy.linspace(2,14,25)
angle_range = numpy.linspace(-1.0,0.5,16)
#threshs = angle_range
max_sphere = 40.0
radii_map = {'CYS': 3.029584442714009, 'SEC': 3.029584442714009,'ILE': 3.3236043438778844, 'SER': 2.942862799880502, 'GLN': 3.392265173640708, 'LYS': 3.4323769248759275, 'TRP': 4.020778501571494, 'PRO': 3.1681967871564782, 'THR': 3.1209640826827605, 'PHE': 3.7193696421024685, 'ALA': 2.8009640812798096, 'GLY': 2.5743877263600927, 'HIS': 3.534660601593296, 'ASN': 3.2435252870320834, 'LEU': 3.3236043438778844, 'ARG': 3.631344304999469, 'ASP': 3.236792124221129, 'VAL': 3.1681967871564782, 'GLU': 3.3861111330098956, 'TYR': 3.8021338661518342, 'MET': 3.351107764635289}

avg_cent_side = {'ILE': (-0.90613294951389123, -1.2696649328274769, 0.75787678394555402), 'GLN': (-1.2338776755913703, -1.6271318863182545, 0.98235396198659397), 'GLY': (-0.0047299661865784631, -0.0088845691075302401, 9.8567539791531887e-05), 'GLU': (-1.2149070583116599, -1.5273224890561412, 1.0878295037481032), 'CYS': (-0.58968190337597537, -0.89958163307781258, 0.62776825003201231), 'HIS': (-0.9776890826599115, -1.6575223413455686, 0.78038430859923047), 'SER': (-0.59170713881015236, -0.65044117400421442, 0.79631448621746626), 'LYS': (-1.4201311703581181, -1.8190188812624906, 1.0949225828431293), 'PRO': (0.21186055970214385, -0.81932530473289122, 1.0393498184079177), 'SEC': (-0.83941591976390406, -0.71433578942505138, 0.85735910033472262), 'ASN': (-0.84857508137952831, -1.3074529598975138, 0.7776964830504921), 'VAL': (-0.79991957561390237, -0.94684856115170413, 0.62104228622222646), 'THR': (-0.65195583478889552, -0.90724937534582129, 0.82850637927897564), 'ASP': (-0.90583763237409776, -1.1848068294814498, 0.82225350991223289), 'TRP': (-1.528171801920162, -1.7410619994461729, 0.89554186538385805), 'UNK': (-0.3674831972563582, -0.32646621880413212, 0.49726141188550765), 'PHE': (-1.1835851942281017, -1.7644149833489573, 0.692845774001975), 'ALA': (-0.40699060927186675, -0.36885898636123338, 0.49092580988001461), 'MET': (-1.1860699223797722, -1.5396420763646954, 0.84853591931547168), 'LEU': (-1.0557131262899473, -1.4510720809764079, 0.67350078632184163), 'ARG': (-1.7179619489063336, -2.0937800414713896, 1.3885442733060744), 'TYR': (-1.2761868751958905, -1.9287097734708729, 0.69640859762760676)}

unknown_avg_cent_side = (-0.866137516594, -1.16682367833, 0.769774990962) #average vector over all residue types
unknown_thresh = 7.5
unknown_angle_thresh = -1.0
unknown_rad = 3.23342916549


#Residue distribution:  {'CYS': 0.012580215522460346, 'ASP': 0.0612362271461436, 'SER': 0.059280784598619685, 'VAL': 0.07643782540259111, 'GLN': 0.034986075796101225, 'LYS': 0.05867538442910764, 'ILE': 0.06109093110546071, 'PRO': 0.04416999636759898, 'THR': 0.05295435282721879, 'PHE': 0.04113694151834363, 'ALA': 0.08321225329943092, 'GLY': 0.07504540501271341, 'HIS': 0.022236348226177503, 'GLU': 0.06778060297856883, 'LEU': 0.09070710739799007, 'ARG': 0.05026032207289018, 'TRP': 0.011520765225814264, 'ASN': 0.042135851798038505, 'TYR': 0.03321225329943092, 'MET': 0.021340355975299673}

avg_cent = {'ILE': (-0.51939923132192789, -0.42924346469309471, 0.57478473070875846), 'GLN': (-0.80834986146513732, -0.76972077063003586, 0.68335543707029034), 'GLY': (0.003710108790125123, 0.74850253699086222, -0.011107678632392972), 'GLU': (-0.78586343234711975, -0.68596611605471958, 0.72074313020193725), 'CYS': (-0.25453284884232952, 0.025290152329503414, 0.46347585732429047), 'HIS': (-0.67516717655691016, -0.87046086894289898, 0.61423651641973831), 'SER': (-0.25186295200626596, 0.15435437507216301, 0.50821742795023583), 'LYS': (-0.93157547397808582, -0.90270366121918721, 0.7564363100362167), 'PRO': (0.19092956143389894, -0.054932017506649618, 0.7190732935459736), 'SEC': (-0.4041486429034104, 0.19276231247130715, 0.31022206884216508), 'ASN': (-0.52337772135489558, -0.46464652365549852, 0.5839222978918297), 'VAL': (-0.40625381747595429, -0.12555581649592371, 0.49308718269365215), 'THR': (-0.3458279994581282, -0.10156321384665247, 0.56845475732264672), 'ASP': (-0.5395255947207902, -0.37285940531331668, 0.57492539764983719), 'TRP': (-1.1801359800054108, -1.1546490468404469, 0.74010431975141766), 'UNK': (-0.16216122024915333, 0.48606836524146568, 0.13849611887759408), 'PHE': (-0.83581175433622168, -1.0098420145081217, 0.57161519965229146), 'ALA': (-0.11519817930399062, 0.43429992078978386, 0.21401964289355344), 'MET': (-0.70928160313770094, -0.59184549532823161, 0.58483387893500882), 'LEU': (-0.62691092485685773, -0.54401570724965409, 0.4676475183702658), 'ARG': (-1.2457704954369544, -1.2737121830387164, 1.037835446531018), 'TYR': (-0.9381954559286525, -1.1994693553734199, 0.57799294311076543)}

average_correction = {'CYS': -0.022946150704117624, 'SEC': -0.022946150704117624, 'GLN': -0.02961937816659793, 'ASP': -0.013910148202608552, 'SER': 0.013356648242644154, 'VAL': -0.03795902095389797, 'LYS': -0.03942889003072603, 'ILE': -0.05977847074136369, 'PRO': -0.004513042646024423, 'THR': -0.009849324564857875, 'PHE': -0.04798429467867554, 'ALA': 0.009567509980448752, 'GLY': 0.01743352671606324, 'HIS': -0.022794305142026405, 'GLU': -0.044020993580337175, 'LEU': -0.0641172463977473, 'ARG': -0.0848773838260547, 'TRP': -0.10482640485221231, 'ASN': 0.012512801252194339, 'TYR': -0.04262395441484643, 'MET': -0.059906004162288945}

average_correction_alpha = {'CYS': -0.08178375570565236, 'SEC': -0.08178375570565236, 'ILE': -0.0763942769641649, 'SER': -0.056276272825953266, 'GLN': -0.05432425352766346, 'LYS': -0.04347620974426458, 'TRP': -0.13120029100479658, 'PRO': -0.02893434811025425, 'THR': -0.06906454564366421, 'PHE': -0.08442494055881009, 'ALA': -0.05209043431185281, 'GLY': -0.029668853492853364, 'HIS': -0.07807760374070327, 'ASN': -0.04205405440700721, 'LEU': -0.11715899801205235, 'ARG': -0.08525007619343955, 'ASP': -0.05334347102978171, 'VAL': -0.09920832263668614, 'GLU': -0.05777287475315668, 'TYR': -0.08999231904470734, 'MET': -0.10386081073354511}


#values learned on old training set
#threshs = {'CYS' : 7.25, 'SEC': 7.25, 'GLN' : 7.25, 'ASP' : 6.75, 'SER' : 6.5, 'VAL' : 7.5, 'LYS' : 7.0, 'ASN' : 7.0, 'PRO' : 7.25, 'THR' : 7.0, 'PHE' : 8.25, 'ALA' : 6.75, 'HIS' : 7.25, 'GLY' : 6.75, 'ILE' : 7.5 , 'LEU' : 7.75, 'ARG' : 7.75, 'TRP' : 8.25, 'GLU' : 7.0, 'TYR' : 8.0, 'MET' : 7.5}

#angle_threshs = {'CYS' : -0.85, 'SEC' : -0.85, 'GLN' : -1.0, 'ASP' : -0.95, 'SER' : -0.9, 'VAL' : -0.85, 'LYS' : -1.0, 'ASN' : -0.9, 'PRO' : -0.9, 'THR' : -0.9, 'PHE' : -0.8, 'ALA' : -0.85, 'HIS' : -0.75, 'GLY' : -1.0, 'ILE' : -0.9, 'LEU' : -0.9, 'ARG' : -1.0, 'TRP' : -0.95, 'GLU' : -1.0, 'TYR' : -0.75, 'MET' : -0.95}

#threshs_alpha = {'CYS' : 7.25, 'SEC' : 7.25,'ILE' : 7.75, 'SER' : 6.75, 'GLN' : 8.0, 'LYS' : 8.0, 'PRO' : 7.25, 'ASP' : 7.0, 'THR' : 7.0, 'PHE' : 8.5, 'ALA' : 7.0, 'GLY' : 6.75, 'HIS' : 8.0, 'GLU' : 7.75, 'LEU' : 7.75, 'ARG' : 8.75, 'TRP' : 8.75, 'VAL' : 7.5, 'ASN' : 7.25, 'TYR' : 8.75, 'MET' : 8.0}

#angle_threshs_alpha = {'CYS' : -0.9, 'SEC' : -0.9, 'GLN' : -1.0, 'ASP' : -1.0, 'SER' : -1.0, 'VAL' : -1.0, 'LYS' : -1.0, 'ASN' : -1.0, 'PRO' : -1.0, 'THR' : -1.0, 'PHE' : -0.85, 'ALA' : -0.9, 'HIS' : -1.0, 'GLY' : -1.0, 'ILE' : -0.85, 'LEU' : -1.0, 'ARG' : -1.0, 'TRP' : -1.0, 'GLU' : -1.0, 'TYR' : -0.95, 'MET' : -1.0}

#learned on golden train set without included target residue
#"""
threshs_alpha ={
'CYS' : 7.25, 
'ASP' : 7.25,
'SER' : 7.0,
'VAL' : 7.5,
'GLN' : 7.75,
'LYS' : 8.0,
'PRO' : 7.5,
'THR' : 7.0,
'PHE' : 8.25,
'ALA' : 7.0,
'HIS' : 8.0,
'GLY' : 6.75,
'ILE' : 7.75,
'GLU' : 7.5,
'LEU' : 7.75,
'ARG' : 8.5,
'TRP' : 9.0,
'ASN' : 7.25,
'TYR' : 8.75,
'MET' : 8.0
}

angle_threshs_alpha = {
'CYS' : -0.85,
'ASP' : -0.90,
'SER' : -0.85,
'VAL' : -0.90,
'GLN' : -1.0,
'LYS' : -1.0,
'PRO' : -1.0,
'THR' : -1.0,
'PHE' : -0.90,
'ALA' : -0.85,
'HIS' : -0.95,
'GLY' : -1.0,
'ILE' : -0.85,
'GLU' : -1.0,
'LEU' : -1.0,
'ARG' : -1.0,
'TRP' : -0.6,
'ASN' : -1.0,
'TYR' : -0.85,
'MET' : -1.0
}
#"""
"""
#learned on golden train with including intersecting residue spheres
threshs_alpha ={
'CYS' : 8.5,
'ASP' : 9.0,
'SER' : 9.0,
'VAL' : 9.0,
'GLN' : 9.5,
'LYS' : 9.5,
'PRO' : 9.0,
'THR' : 9.0,
'PHE' : 10.0,
'ALA' : 9.0,
'HIS' : 9.5,
'GLY' : 8.5,
'ILE' : 10.0,
'GLU' : 9.5,
'LEU' : 9.5,
'ARG' : 10.0,
'TRP' : 10.5,
'ASN' : 9.0,
'TYR' : 10.5,
'MET' : 9.5
}
angle_threshs_alpha = {
'CYS' : -0.9,
'ASP' : -0.9,
'SER' : -0.8,
'VAL' : -0.9,
'GLN' : -1.0,
'LYS' : -0.9,
'PRO' : -1.0,
'THR' : -0.9,
'PHE' : -0.4,
'ALA' : -0.8,
'HIS' : -0.9,
'GLY' : -1.0,
'ILE' : -0.8,
'GLU' : -0.9,
'LEU' : -0.9,
'ARG' : -0.9,
'TRP' : -0.4,
'ASN' : -0.9,
'TYR' : -0.4,
'MET' : -1.0
}
"""
"""
#learned on golden train set with target residue included
threshs_alpha ={
'CYS' : 7.25,
'ASP' : 7.25,
'SER' : 7.0,
'VAL' : 7.5,
'GLN' : 7.75,
'LYS' : 8.0,
'PRO' : 7.5,
'THR' : 7.0,
'PHE' : 8.25,
'ALA' : 7.0,
'HIS' : 8.0,
'GLY' : 6.75,
'ILE' : 7.75,
'GLU' : 7.5,
'LEU' : 7.75,
'ARG' : 8.5,
'TRP' : 9.0,
'ASN' : 7.25,
'TYR' : 8.75,
'MET' : 8.0}

angle_threshs_alpha = {
'CYS' : -0.85,
'ASP' : -0.90,
'SER' : -0.85,
'VAL' : -0.90,
'GLN' : -1.0,
'LYS' : -1.0,
'PRO' : -1.0,
'THR' : -1.0,
'PHE' : -0.90,
'ALA' : -0.85,
'HIS' : -0.95,
'GLY' : -1.0,
'ILE' : -0.85,
'GLU' : -1.0,
'LEU' : -1.0,
'ARG' : -1.0,
'TRP' : -0.60,
'ASN' : -1.0,
'TYR' : -0.85,
'MET' : -1.0
}
"""
#Best threshs for c_alpha and unknown target residue: 7.5,-1.0

#learned values from golden train set
threshs = {'CYS' : 7.0, 'SEC': 7.0, 'GLN' : 7.25, 'ASP' : 6.75, 'ASX': 6.75, 'SER' : 6.5, 'VAL' : 7.5, 'LYS' : 7.5, 'ASN' : 6.75, 'PRO' : 7.0, 'THR' : 6.75, 'PHE' : 7.75, 'ALA' : 7.0, 'HIS' : 7.25, 'GLY' : 6.75, 'ILE' : 7.75 , 'LEU' : 7.75, 'ARG' : 7.75, 'TRP' : 8.0, 'GLU' : 7.0, 'GLX': 7.0, 'TYR' : 7.75, 'MET' : 7.5, 'UNK': unknown_thresh}

angle_threshs = {'CYS' : -0.85, 'SEC' : -0.85, 'GLN' : -0.95, 'ASP' : -0.95, 'ASX' : -0.95, 'SER' : -0.9, 'VAL' : -0.85, 'LYS' : -0.9, 'ASN' : -0.95, 'PRO' : -0.85, 'THR' : -0.85, 'PHE' : -0.85, 'ALA' : -0.85, 'HIS' : -0.8, 'GLY' : -1.0, 'ILE' : -0.9, 'LEU' : -0.9, 'ARG' : -0.95, 'TRP' : -0.95, 'GLU' : -0.95, 'GLX' : -0.95, 'TYR' : -0.85, 'MET' : -1.0, 'UNK': unknown_angle_thresh}

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
        'THR': 142.0, 'TRP': 227.0, 'TYR': 222.0, 'VAL': 142.0,
        'SEC': 135.0 
    } 
} 

def calcVol(r,cos):
    vol = (2.0/3.0)*math.pi*(r**3.0)*(1.0-cos)
    return vol

def sissCorrectionVol(siss,r,cos):
    v = calcVol(r,cos)
    return (v-siss)/v

def sissCorrection(siss,res):
    v = cone_cut_volumes[res]
    max_acc = residue_max_acc['Sander'][res]
    #corrected_siss = siss/(v*max_acc)
    #corrected_siss = siss/(v)
    corrected_siss = (v-siss)/v
    #corrected_siss = siss/(max_acc)
    #corrected_siss = siss/(v+max_acc)
    return corrected_siss

def parsePDB(input_file,chain,c_alpha):
    """
    Parses a PDB-file and takes all atomic coordinates of a specified chain.

    Input:
    input_file: String ; Path to a PDB file
    chain: String or None ; Chain identifier, if None is given, the first Chain found in the file is taken
    c_alpha: Boolean ; If True, then only C alpha atoms are taken

    Output:
    coordinate_map: {String:[String,{String:(String,float,float,float)}]} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.   
    """
    f = open(input_file,'r')
    lines = f.readlines()
    f.close()
    coordinate_map = {}

    x_total = 0.0
    y_total = 0.0
    z_total = 0.0
    n = 0.0

    #print lines

    for line in lines:
        if len(line) > 5:        
            record_name = line[0:6].replace(" ","")
            if record_name == "ENDMDL":
                #print len(coordinate_map)
                break
        #ignore short lines
        if len(line) > 20:
            atom_nr = line[6:11].replace(" ","")
            atom_name = line[12:16].replace(" ","")
            res_name = line[17:20].replace(" ","")
            
            #len(res_name) == 3 is only true for amino acid chains
            if len(res_name) == 3:
                if len(line) > 21:
                    chain_id = line[21]
                    res_nr = int(line[22:27].replace(" ",""))
                    #insertion_code = line[26]
                
                #consider only lines with record name ATOM
                if record_name == "ATOM":
                    #if chain not given, take first one found
                    if chain == None:
                        chain = chain_id
                    if len(line) > 50:
                        #print chain,chain_id
                        #consider only lines with the correct chain id
                        if chain_id == chain:
                            #if c_alpha is true, then take only C alpha atoms
                            if res_name != 'UNK':
                                if (not c_alpha) or atom_name == 'CA':
                                    if atom_name[0] != 'H' and atom_name[0] != 'D':
                                        x = float(line[30:38].replace(" ",""))
                                        y = float(line[38:46].replace(" ",""))
                                        z = float(line[46:54].replace(" ",""))
                                        if res_nr not in coordinate_map:
                                            coordinate_map[res_nr]=[res_name,{}]
                                        coordinate_map[res_nr][1][atom_nr] = (atom_name,x,y,z)
                                        x_total += x
                                        y_total += y
                                        z_total += z
                                        n += 1.0
    if n > 0.0:                                        
        protein_centroid = numpy.array([x_total/n,y_total/n,z_total/n])
    else:
        protein_centroid = numpy.array([0.0,0.0,0.0])
    #print coordinate_map
    return coordinate_map,protein_centroid

def calcCentroidMap(coordinate_map,target_residues,c_alpha,double_unknown_mode = False):
    #print coordinate_map
    centroid_map = {}
    if c_alpha:
        c1 = None
        c2 = None
        res_2 = None
        res_1 = None
        res_name_1 = None
        res_name_2 = None
        for res in coordinate_map:
            c0 = c1
            c1 = c2
            res_0 = res_1
            res_1 = res_2
            res_2 = res
            res_name_1 = res_name_2
            res_name_2 = coordinate_map[res_2][0]
            atomlist = coordinate_map[res_2][1]
            #if len(atomlist) > 1:
            #    print atomlist
            #    raise NameError('More than one atom in C alpha only case.')
            (atomname,x,y,z) = atomlist.values()[0]
            if atomname != ('CA'):
                raise NameError('Only atom is not C alpha')
            c2 = numpy.array([x,y,z])           
            if res_0 != None and res_1 != None:
                if res_2 - res_1 == 1 and res_1 - res_0 == 1:
                    avg_cent_side_vec = avg_cent_side[res_name_1]
                    avg_cent_vec = avg_cent[res_name_1]
                    if double_unknown_mode: #in double unknown mode, the residue types are not needed!
                        avg_cent_side_vec = unknown_avg_cent_side
                    side_centroid = predict_centroid(c0,c1,c2,avg_cent_side_vec)
                    centroid = predict_centroid(c0,c1,c2,avg_cent_vec)
                    
                    centroid_map[res_1] = (side_centroid,centroid)

                else:
                    centroid_map[res_1] = (c1,c1)
            elif res_1 != None:
                centroid_map[res_1] = (c1,c1)
        if len(coordinate_map) > 0:
            centroid_map[res_2] = (c2,c2)
    else:
        for res in target_residues:
            atomlist = coordinate_map[res][1]
            centroid = getCentroid(atomlist)
            centroid_map[res] = centroid

    return centroid_map

"""
def calcDistMatrix(coordinate_map,target_residues):
    
    Calculates distances matrix for all atoms in a set of target residues and all atoms in a coordinate map. Distances between atoms of the same residue are not computed
    
    Input:
    coordinate_map: {String:[String,{String:(String,float,float,float)}]} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.   
    target_residues: [String] or None; contains the residue-ids for all residues, for which the distances are needed. If None is given, compute all distances inbetween the coordinate_map.
    
    Output:
    dist_matrix: {String:{String:Float}} ; Matrix of [Atom-ID][Atom-ID] -> euclidean distance. IDs are not ordered!
    
    dist_matrix= {}

    #If None is given, take all residues from the coordinate_map
    if target_residues == None:
        target_residues = coordinate_map.keys()
    for res_id_1 in target_residues:
        for atom_id_1 in coordinate_map[res_id_1][1]:
            (atom_name,x,y,z) = coordinate_map[res_id_1][1][atom_id_1]
            vec1 = numpy.array([x,y,z])
            if not atom_id_1 in dist_matrix:
                dist_matrix[atom_id_1 = {}
            for res_id_2 in coordinate_map:
                #ignore if same residue
                if res_id_2 != res_id_1:
                    for atom_id_2 in coordinate_map[res_id_2][1]:
                        #check if distance already computed
                        if not atom_id_2 in dist_matrix[atom_id_1]:
                            if not atom_id_2 in dist_matrix:
                                (atom_name,x,y,z) = coordinate_map[res_id_2][1][atom_id_2]
                                vec2 = numpy.array([x,y,z])
                                diff = vec2 - vec1
                                d = numpy.sqrt(numpy.dot(diff, diff))
                                dist_matrix[atom_id_1][atom_id_2] = d
                            elif not atom_id_1 in dist_matrix[atom_id_2]:
                                (atom_name,x,y,z) = coordinate_map[res_id_2][1][atom_id_2]
                                vec2 = numpy.array([x,y,z])
                                diff = vec2 - vec1
                                d = numpy.sqrt(numpy.dot(diff, diff))
                                dist_matrix[atom_id_1][atom_id_2] = d

    return dist_matrix 
"""

def calcDistMatrix(coordinate_map,centroid_map,target_residues,c_alpha):
    dist_matrix = {}
    if c_alpha:
        for res_1 in centroid_map:
            #dist_matrix[res_1] = {}
            centroid_1 = centroid_map[res_1][0]
            for res_2 in centroid_map:
                if not (res_1,res_2) in dist_matrix:
                    centroid_2 = centroid_map[res_2][0]
                    #if centroid_2 == None:
                    #    print centroid_map
                    diff = centroid_2 - centroid_1
                    d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                    dist_matrix[(res_1,res_2)] = d
                    dist_matrix[(res_2,res_1)] = d
    else:
        for res_1 in target_residues:
            centroid = centroid_map[res_1]
            res_name = coordinate_map[res_1][0]
            thresh = threshs[res_name]
            dist_matrix[res_1] = {}
            for res_2 in coordinate_map:
                #if res_2 != res_1: #if this is commented, the atoms of the target residue are included in the measure calculation
                dist_matrix[res_1][res_2] = {}
                atomlist = coordinate_map[res_2][1]
                for atom in atomlist:
                    (atomname,x,y,z) = atomlist[atom]
                    diff = centroid - numpy.array([x,y,z])
                    #d = numpy.sqrt(numpy.dot(diff, diff))
                    d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                    #If one atom is further than 10A+threshhold, then ignore the whole residue
                    if d - 10.0 > thresh:
                        break

                    dist_matrix[res_1][res_2][atom] = (d,atomname,x,y,z)
    return dist_matrix

def get_gly_cb_vector(residue): 
    """ 
    Return a pseudo CB vector for a Gly residue. 
    The pseudoCB vector is centered at the origin. 

    CB coord=N coord rotated over -120 degrees 
    along the CA-C axis. 
    """ 
    try: 
      n_v = residue["N"].get_vector() 
      c_v = residue["C"].get_vector() 
      ca_v = residue["CA"].get_vector() 
    except: 
      return None 
    # center at origin 
    n_v = n_v - ca_v 
    c_v = c_v - ca_v 
    # rotation around c-ca over -120 deg 
    rot = rotaxis(-math.pi * 120.0 / 180.0, c_v) 
    cb_at_origin_v = n_v.left_multiply(rot)
    vec = numpy.array([cb_at_origin_v[0],cb_at_origin_v[1],cb_at_origin_v[2]])
    return vec


def sphere_intersection(R,r,d):
    #print R,r,d,calcVol(R,-1.0),calcVol(r,-1.0)
    if R < r:
        a = r
        r = R
        R = a
    if r+R <= d:
        return 0.0
    if d+r <= R:
        return (4.0/3.0)*math.pi*r**3.0
    sum1 = (R+r-d)**2
    sum2 = (d**2.0+2.0*d*r-3.0*r**2.0+2.0*d*R+6.0*r*R-3.0*R**2.0)
    si = math.pi*sum1*sum2/(12.0*d)
    #print si
    return si

def createRotMatrix(axis,cos):
    sin = (1.0-cos**2.0)**(1.0/2.0)
    if axis == "x":
        Rot = [[1.0,0.0,0.0],[0.0,cos,-sin],[0.0,sin,cos]]
    if axis == "y":
        Rot = [[cos,0.0,sin],[0.0,1.0,0.0],[-sin,0.0,cos]]
    if axis == "z":
        Rot = [[cos,-sin,0.0],[sin,cos,0.0],[0.0,0.0,1.0]]
    return numpy.matrix(Rot)

def createPlaneRotMatrix(plane,vec):
    x = vec[0]
    y = vec[1]
    z = vec[2]
    if plane == 'xy':
        if z == 0.0:
            cos = 1.0
        elif y == 0.0:
            cos = 0.0
        else:
            cos = (((z**2.0/y**2.0))+1.0)**(-1.0/2.0)
        #sin = (1.0-cos**2.0)**(1.0/2.0)
        rot = createRotMatrix('x',cos)
        if (vec*rot).A1[2]**2 > 0.0000001:
            rot = rot.T
    if plane == 'xz':
        if y == 0.0:
            cos = 1.0
        elif x == 0.0:
            cos = 0.0
        else:
            cos = (((y**2.0/x**2.0))+1.0)**(-1.0/2.0)
        #sin = (1.0-cos**2.0)**(1.0/2.0)
        rot = createRotMatrix('z',cos)
        if (vec*rot).A1[1]**2 > 0.0000001:
            rot = rot.T
    if plane == 'yz':
        if x == 0.0:
            cos = 1.0
        elif z == 0.0:
            cos = 0.0
        else:
            cos = (((x**2.0/z**2.0))+1.0)**(-1.0/2.0)
        #sin = (1.0-cos**2.0)**(1.0/2.0)
        rot = createRotMatrix('y',cos)
        if (vec*rot).A1[0]**2 > 0.0000001:
            rot = rot.T
    return rot

def createAxisRotMatrix(axis,vec):
    if axis == 'x':
        rot1 = createPlaneRotMatrix('xz',vec)
        ivec = (vec*rot1).A1
        rot2 = createRotMatrix('y',getCos('x',ivec))
        if (ivec*rot2).A1[2]**2 > 0.0000001:
            rot2 = rot2.T
    if axis == 'y':
        rot1 = createPlaneRotMatrix('xy',vec)
        ivec = (vec*rot1).A1
        rot2 = createRotMatrix('z',getCos('y',ivec))
        if (ivec*rot2).A1[0]**2 > 0.0000001:
            rot2 = rot2.T
    if axis == 'z':
        rot1 = createPlaneRotMatrix('yz',vec)
        ivec = (vec*rot1).A1
        rot2 = createRotMatrix('x',getCos('z',ivec))
        if (ivec*rot2).A1[1]**2 > 0.0000001:
            rot2 = rot2.T
    #print vec,'vec'
    #print (vec*rot1.T).A1
    #print ivec,'ivec'
    return rot1*rot2

def createRotAxisMatrix(axis,cos):
    sin = (1.0-cos**2.0)**(1.0/2.0)
    axis = axis/numpy.linalg.norm(axis)
    x = axis[0]
    y = axis[1]
    z = axis[2]
    m1 = [cos+x**2.0*(1.0-cos),x*y*(1.0-cos)-z*sin,x*z*(1.0-cos)+y*sin]
    m2 = [y*x*(1.0-cos)+z*sin,cos+y**2.0*(1.0-cos),y*z*(1.0-cos)-x*sin]
    m3 = [z*x*(1.0-cos)-y*sin,z*y*(1.0-cos)+x*sin,cos+z**2.0*(1.0-cos)]
    Rot = [m1,m2,m3]
    return numpy.matrix(Rot)

def gly_vector(n_v,c_v,ca_v):
    n_v = n_v - ca_v 
    c_v = c_v - ca_v
    rot = createRotAxisMatrix(c_v,-0.5)
    vec = (n_v*rot).A1
    return vec
    
def getCosAngle(vec1,vec2):
    n1 = numpy.linalg.norm(vec1)
    n2 = numpy.linalg.norm(vec2)
    norm = n1*n2
    dot = numpy.dot(vec1,vec2)
    if norm != 0.0:
        c = dot / norm
    else:
        return None
    #print c 
    # Take care of roundoff errors 
    c = min(c, 1.0) 
    c = max(-1.0, c)
    return c

def getCos(axis,vec):
    if axis == "x":
        other = numpy.array([1.0,0.0,0.0])
    if axis == "y":
        other = numpy.array([0.0,1.0,0.0])
    if axis == "z":
        other = numpy.array([0.0,0.0,1.0])
    n1 = numpy.linalg.norm(vec)
    n2 = numpy.linalg.norm(other)
    #print n1,n2
    c = ((numpy.dot(vec,other)) / (n1 * n2))
    #print c 
    # Take care of roundoff errors 
    c = min(c, 1.0) 
    c = max(-1.0, c)
    return c

def nullTest(vec):
    if vec[0] == 0.0 and vec[1] == 0.0 and vec[2] == 0.0:
        return False
    return True

def predict_centroid(c0,c1,c2,avg_cent_vec):
    if nullTest(c0-c1) and nullTest(c2-c1):
        A = (c0-c1)/numpy.linalg.norm(c0-c1)
        B = (c2-c1)/numpy.linalg.norm(c2-c1)
        rx = createAxisRotMatrix("x",A)
        B_prime = (B*rx).A1
        r2 = createPlaneRotMatrix('xy',B_prime)
        ROT = rx*r2
        support = (B_prime*r2).A1
        if support[1] < 0.0:
            flip = createRotMatrix('x',-1.0)
            ROT = ROT*flip

        return ((avg_cent_vec*ROT.T).A1)+c1
    else:
        return c1

def getCentroid(atomlist):
    n = 0.0
    t_x = 0.0
    t_y = 0.0
    t_z = 0.0
    c_a = None
    c_b = None
    #print atomlist
    """
    #This for-loop can be used, if C_alpha should be used as Sphere Center
    for atom in atomlist:
        (atomname,x,y,z) = atomlist[atom]
        if atomname == 'CA':
            return numpy.array([x,y,z])
    """
    #"""
    #If this part is commented out, the centroid is calculated by including the backbone
    for atom in atomlist:
        (atomname,x,y,z) = atomlist[atom]
        if len(atomname) > 1:
            t_x += x
            t_y += y
            t_z += z
            n += 1.0
            
    if n > 0.0:
        centroid = numpy.array([t_x/n,t_y/n,t_z/n])
        
    else:
    #"""
        for atom in atomlist:
            (atomname,x,y,z) = atomlist[atom]        
            t_x += x
            t_y += y
            t_z += z
            n += 1.0
            
        if n > 0.0:
            centroid = numpy.array([t_x/n,t_y/n,t_z/n])
        else:
            #Does this happen?
            print "This does happen"
            print atomlist
            centroid = None
    return centroid

def produceOutput(siss_map,coordinate_map,output_file):
    """
    Input:
    siss_map: {String:float} ; 
    coordinate_map: {String:[String,{String:(String,float,float,float)}]} ; Maps residue-id on residue name and atom map. atom map maps atom-id on atom name and atomic coordinates.
    """ 

    lines = []    

    for res in siss_map:
        res_name = coordinate_map[res][0]
        siss = siss_map[res]
        lines.append("%s %s\t%s" % (res,res_name,str(siss)))

    lines = sorted(lines,key=lambda x:int(x.split()[0]))

    lines = ["Residue\tSISS value"] + lines

    f = open(output_file,'w')
    f.write("\n".join(lines))
    f.close()

def calculateSiss(coordinate_map,centroid_map,dist_matrix,target_residues,c_alpha,manu_thresh=None,manu_angle_thresh=None,double_unknown_mode=False,manu_parameters=None,unknown_parameters=None):
    global threshs
    global angle_threshs
    global threshs_alpha
    global angle_threshs_alpha
    global unknown_thresh
    global unknown_angle_thresh
    siss_map = {}

    if unknown_parameters != None:
        (unknown_thresh,unknown_angle_thresh) = unknown_parameters
    if manu_parameters != None:
        [threshs,angle_threshs] = manu_parameters
        [threshs_alpha,angle_threshs_alpha] = manu_parameters
    if c_alpha:
        for res in target_residues:
            res_name = coordinate_map[res][0]

            if manu_thresh == None:
                thresh = threshs_alpha[res_name]
            else:
                thresh = manu_thresh
            if manu_angle_thresh == None:
                angle_thresh = angle_threshs_alpha[res_name]
            else:
                angle_thresh = manu_angle_thresh

            if thresh == None or double_unknown_mode:
                thresh = unknown_thresh
            if angle_thresh == None or double_unknown_mode:
                angle_thresh = unknown_angle_thresh

            (atomname,x,y,z) = coordinate_map[res][1].values()[0]
            c_alpha_1 = numpy.array([x,y,z])
            centroid_1 = centroid_map[res][0]

            
            #test: substract ssi of all residue spheres from the total sum
            #res_list = []

            siss = 0.0
            for res_2 in coordinate_map:
                #if res == res_2:
                #    siss += sphere_intersection(thresh,radii_map[res_name],0.0)
                #elif res != res_2:
                if res != res_2:
                    res_name_2 = coordinate_map[res_2][0]
                    if (res,res_2) in dist_matrix:
                        d = dist_matrix[(res,res_2)]
                    else:
                        raise NameError('Did not find residue pair in distance matrix for: %s and %s' % (res,res_2))

                    rad = radii_map[res_name_2]
                    if double_unknown_mode:
                        rad = unknown_rad
                    if d - rad <= thresh:
                        if angle_thresh > -1.0: 
                            centroid_2 = centroid_map[res_2][0]
                            
                            angle = getCosAngle(centroid_1-c_alpha_1,centroid_2-c_alpha_1)
                            #print centroid_1,centroid_2,c_alpha_1,angle
                            if angle == None:    
                                siss += sphere_intersection(thresh,rad,d)
                                #res_list.append(res_2)
                            elif angle_thresh <= angle:
                                siss += sphere_intersection(thresh,rad,d)
                                #res_list.append(res_2)
                        else:
                            siss += sphere_intersection(thresh,rad,d)
                            #res_list.append(res_2)

            """
            #test: substract ssi of all residue spheres from the total sum
            done = set([])
            for r1 in res_list:
                done.add(r1)
                for r2 in res_list:
                    if r2 in done:
                        continue
                    if (r1,r2) in dist_matrix:
                        d = dist_matrix[(r1,r2)]
                    else:
                        centroid_1 = centroid_map[r1][0]
                        centroid_2 = centroid_map[r2][0]
                        diff = centroid_2 - centroid_1
                        d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                    siss -= sphere_intersection(radii_map[coordinate_map[r1][0]],radii_map[coordinate_map[r2][0]],d)
            """

            siss = sissCorrectionVol(siss,thresh,angle_thresh)
            #siss_map[res] = siss + average_correction_alpha[res_name]
            siss_map[res] = siss
    else:
        for res in target_residues:
            res_name = coordinate_map[res][0]
            atomlist = coordinate_map[res][1]
            if manu_thresh == None:
                thresh = threshs[res_name]
            else:
                thresh = manu_thresh
            if manu_angle_thresh == None:
                angle_thresh = angle_threshs[res_name]
            else:
                angle_thresh = manu_angle_thresh

            if thresh == None or double_unknown_mode:
                thresh = unknown_thresh
            if angle_thresh == None or double_unknown_mode:
                angle_thresh = unknown_angle_thresh
            centroid_1 = centroid_map[res]
            #centroid_1 = c_alpha_1

            gly_c = [None]
            gly_n = [None]
            c_alpha_1 = [None]
            c_beta = [None]
            for atom in atomlist:
                (atomname,x,y,z) = atomlist[atom]
                if atomname == 'CA':
                    c_alpha_1 = numpy.array([x,y,z])
                if atomname == 'CB':
                    c_beta = numpy.array([x,y,z])
                if atomname == 'N':
                    gly_n = numpy.array([x,y,z])
                if atomname == 'C':
                    gly_c = numpy.array([x,y,z])

            siss = 0.0
            error_flag = False

            for res_2 in dist_matrix[res]:
                for atom in dist_matrix[res][res_2]:
                    (d,atomname,x,y,z) = dist_matrix[res][res_2][atom]
                    #print d,res,res_2
                    #if res == res_2:
                        #print 'this should happen alot'
                    #    siss += sphere_intersection(thresh,vdw_radius[atomname[0]],d)
                    if d - vdw_radius[atomname[0]] <= thresh:
                        if angle_thresh > -1.0:
                            coord = numpy.array([x,y,z])

                            if res_name != 'GLY':
                                if c_beta[0] == None or c_alpha_1[0] == None:
                                    angle = None
                                else:
                                    angle = getCosAngle(c_beta-c_alpha_1,coord-c_alpha_1)

                                #angle = None #Calculate everything with gly-vector

                                if angle == None:
                                    error_flag = True
                                    #print c_beta,res,res_2,atomlist
                                    #This happens for resiudes (beside Glycin), where only the C-Alpha atom is given. (Note: This is the full atom case)
                                    #The solution is to handle it as Glycin
                                    
                                    if gly_c[0] == None or gly_n[0] == None or c_alpha_1[0] == None:
                                        angle = 1.0
                                    else:
                                        gly_vec = gly_vector(gly_n,gly_c,c_alpha_1)
                                        angle = getCosAngle(gly_vec,coord-centroid_1)
                                
                            else:
                                if gly_c[0] == None or gly_n[0] == None or c_alpha_1[0] == None:
                                    error_flag = True
                                    #print res,res_2
                                    angle = 1.0
                                else:    
                                    gly_vec = gly_vector(gly_n,gly_c,c_alpha_1)
                                    angle = getCosAngle(gly_vec,coord-centroid_1)
                            
                            if angle_thresh <= angle or angle == None:
                                siss += sphere_intersection(thresh,vdw_radius[atomname[0]],d)
                            elif angle_thresh == -1.0:
                                print "WADDAFAQ"
                        else:
                            #If the full sphere is taken, there is no need for calculating the angle
                            siss += sphere_intersection(thresh,vdw_radius[atomname[0]],d)
                        
            #siss = sissCorrection(siss,res_name)
            #if not error_flag: #comment this out for final version
            siss = sissCorrectionVol(siss,thresh,angle_thresh)
            siss_map[res] = siss #+ average_correction[res_name]
            #else:
            #    print "Error_flage was hit"

    return siss_map

def calcAASiss(coordinate_map,target_residues,manu_thresh=None):
    siss_map = {}    
    
    for res in target_residues:
        siss = 0.0
        res_name = coordinate_map[res][0]
        atomlist = coordinate_map[res][1]
        
        if manu_thresh == None:
            thresh = threshs[res_name]
        else:
            thresh = manu_thresh
        
        for atom in atomlist:
            atom_siss = 0.0
            (atomname,x,y,z) = atomlist[atom]    
            for res_2 in coordinate_map:
                if res != res_2:
                    atomlist_2 = coordinate_map[res_2][1]
                    for atom_2 in atomlist_2:
                        (atomname_2,x_2,y_2,z_2) = atomlist_2[atom_2]
                        diff = numpy.array([x_2,y_2,z_2]) - numpy.array([x,y,z])
                        d = math.sqrt(diff[0]**2.0+diff[1]**2.0+diff[2]**2.0)
                        #If one atom is further than 10A+threshhold, then ignore the whole residue
                        if d - 15.0 > thresh:
                            break
                        if d - vdw_radius[atomname_2[0]] <= thresh:                    
                            atom_siss += sphere_intersection(thresh,vdw_radius[atomname_2[0]],d)
                        
        #siss = sissCorrection(siss,res_name)
            siss += sissCorrectionVol(atom_siss,thresh,-1.0)
        siss_map[res] = siss/float(len(atomlist))       

    return siss_map

def siss(input_file=None,output_file=None,target_residues=None,chain=None,c_alpha=False,coordinate_map_path=None,unknown_mode=False,double_unknown_mode=False):
    
    cwd = os.getcwd()
    if output_file == None:
        if input_file == None:
            output_file = "%s/siss_output.tsv" % cwd
        else:
            output_file = "%s_siss.tsv" % input_file

    if input_file == None:
        if coordinate_map_path == None:
            raise NameError("No input structure given")
        else:
            coordinate_map = parseCM(coordinate_map_path)
    else:
        #TODO remove centroid calculation
        coordinate_map,protein_centroid = parsePDB(input_file,chain,c_alpha)


    if target_residues == None:
        target_residues = coordinate_map.keys()

    centroid_map = calcCentroidMap(coordinate_map,target_residues,c_alpha,double_unknown_mode=double_unknown_mode)

    dist_matrix = calcDistMatrix(coordinate_map,centroid_map,target_residues,c_alpha)

    if not unkown_mode:
        siss_map = calculateSiss(coordinate_map,centroid_map,dist_matrix,target_residues,c_alpha)
    else:
        siss_map = calculateSiss(coordinate_map,centroid_map,dist_matrix,target_residues,c_alpha,manu_thresh=unknown_thresh,manu_angle_thresh=unknown_angle_thresh,double_unknown_mode=double_unknown_mode)
    
    produceOutput(siss_map,coordinate_map,output_file)
    
if __name__ == "__main__":
    
    argv = sys.argv[1:]
    helptext = """
\tSphere InterSection Sum

This tool calculates a measure for relative solvent accessible area for single residues or whole amino acid chains.
Usage:
siss.py -i /Path/To/Input/File [-o /Path/To/Output/File] [-c chain] [-r residues] [--ca] [--cm coordinate_map]

-i:\tPath to an input file in PDB (Protein Data Bank) file format. If --cm are given, this is not needed.

-o:\tPath to the output file produced as tab separated text file.
\tDefault: *InputFile*_siss.tsv

-c:\tChain identifier of the amino acid chain, which should be analysed denoted as the chain identifiers of the ATOM records in the PDB file.
\tDefault: The first chain found in the input file

-r:\tList of residue identifiers for all residues for which the SISS value should be computed. If not given, the SISS values for all residues are computed.
\tIf given with -i the residue identifiers denote as the residues identifiers of the ATOM records in the PDB file. If given with --cm the residues should denote as the residue map keys of --cm.
\tExamples: 234,78,368 | 17 | 34,35,36,37

--ca:\tC alpha version of SISS. Needs only coordinates of the C alpha atoms, but is more inaccurate as the full atom coordinates version.

--cm:\tPath to a text file containing an atomic cooridante map as a python map: {Residue-ID:(Residue-Name,{Atom-ID:(Atom-Name,x,y,z)})}.
"""
    if len(argv) == 0:
        print helptext
        sys.exit()

    input_file = None
    output_file = None
    target_residues = None
    chain = None
    c_alpha = False
    coordinate_map_path = None

    try:
        opts,args = getopt.getopt(argv,"hr:o:i:c:",['ca','cm='])
    except getopt.GetoptError:
        print "siss.py -h"
        sys.exit(2)
    for opt,arg in opts:
        if opt == '-h':
            print(helptext)
            sys.exit()
        elif opt == "-i":
            input_file = arg
        elif opt == '-o':
            output_file = arg
        elif opt == '-r':
            target_residues = arg.split(',')
        elif opt == '-c':
            chain = arg
        elif opt == '--ca':
            c_alpha = True
        elif opt == '--cm':
            coordinate_map_path = arg

    siss(input_file,output_file,target_residues,chain,c_alpha,coordinate_map_path)        
