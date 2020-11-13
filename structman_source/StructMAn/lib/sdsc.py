#sdsc: structman datastructures and classes
import rin
import re
import numpy as np

from sys import getsizeof, stderr
import itertools
from collections import deque
try:
    from reprlib import repr
except ImportError:
    pass


MobiDB_map = {'D_PA':'Polyampholite','D_WC':'Weak polyampholie','D_NPE':'D_NPE','D_PPE':'D_PPE'}

# PDB Three-letter-codes (spaces stripped)
# Generated from http://ligand-expo.rcsb.org/dictionaries/Components-smiles-oe.smi
# Using: grep -P "^\[?(Li|Be|Na|Mg|Al|K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|Fr|Ra|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs|Cn)[+-]?[0-9]*\]?\t" Components-smiles-oe.smi
# Does not include metalloids
metal_atoms = set(['0BE','3CO','3NI','4MO','4PU','4TI','6MO','AG','AL',
                    'AM','AU','AU3','BA','BS3','CA','CD','CE','CF',
                    'CO','CR','CS','CU','CU1','CU3','DY','ER3','EU','EU3',
                    'FE','FE2','GA','GD','GD3','HG','HO','HO3','IN','IR','IR3','K',
                    'LA','LI','LU','MG','MN','MN3','MO','NA','NI','OS','OS4',
                    'PB','PD','PR','PT','PT4','RB','RE','RH','RH3','RU','SM','SR',
                    'TA0','TB','TH','TL','U1','V','W','Y1','YB','YB2','YT3','ZCM',
                    'ZN','ZN2','ZR'])

# PDB Three-letter-codes (spaces stripped)
# Generated from http://ligand-expo.rcsb.org/dictionaries/Components-smiles-oe.smi
# Using: grep -P "^\[?(As|Si|F|I|Cl|Sb|Te|At|Br)[+-][0-9]*\]?\t" Components-smiles-oe.smi
# Includes (unused/obsolete) -group entries (ending in O) in addition to ions
ion_atoms = set(['BR', 'BRO', 'CL', 'CLO', 'F', 'FLO', 'IDO', 'IOD', 'SB'])

hydropathy = {'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'U':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5,'X':0.0,'B':0.0}
volume = {'I':166.7,'V':140.,'L':166.7,'F':189.9,'C':108.5,'U':108.5,'M':162.9,'A':88.6,'G':60.1,'T':116.1,'S':89.0,'W':227.8,'Y':193.6,'P':112.7,'H':153.2,'E':138.4,'Q':143.8,'D':111.1,'N':114.1,'K':168.6,'R':173.4,'X':150.0,'B':150.0}

aa_map_aliphatic = set(['I','L','V'])
aa_map_hydrophobic = set(['I','L','V','M','F','Y','W','H','K','T','C','U','A','G'])
aa_map_aromatic = set(['F','Y','W','H'])
aa_map_positive = set(['H','K','R'])
aa_map_polar = set(['Y','W','H','K','T','C','U','D','E','S','N','Q'])
aa_map_negative = set(['D','E'])
aa_map_charged = set(['H','K','R','D','E'])
aa_map_small = set(['P','G','C','U','A','S','N','D','T','V'])
aa_map_tiny = set(['G','C','U','A','S'])

blosum62 = {('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0, ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
            ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1, ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
            ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3, ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
            ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4, ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
            ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2, ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
            ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3, ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
            ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1, ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
            ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2, ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
            ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0, ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
            ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1, ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
            ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3, ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
            ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2, ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
            ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2, ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
            ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1, ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
            ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2, ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
            ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1, ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
            ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0, ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
            ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0, ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
            ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1, ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
            ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1, ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
            ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1, ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
            ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2, ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
            ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1, ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
            ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2, ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
            ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1, ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
            ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3, ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
            ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3, ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
            ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4, ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
            ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2, ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
            ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3, ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
            ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1, ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
            ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3, ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
            ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1, ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
            ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3, ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
            ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4}

def median(l):
    #print l
    n = len(l)
    l = sorted(l)
    if n % 2 == 0:
        med = (l[(n//2)-1]+l[n//2])/2.0
    else:
        med = l[(n-1)//2]
    return med

def processAlignmentData(alignment):
    lines = alignment.split("\n")
    nextline = False
    target_start = False
    template_start = False
    target_lines = []
    template_lines = []
    target_name = ""
    for line in lines:
        if len(line) > 0:
            if line[0] == ">":
                ids = line.split(";")
                if target_name == "":
                    target_name = ids[1]
                    target_name = target_name.replace(" ","")
                    target_name = target_name.replace("\n","")
                nextline = True
            elif nextline:
                if not target_start:
                    target_start = True
                else:
                    target_start = False
                    template_start = True
                    words = line.split(":")
                    startres = words[2]
                    endres = words[4]
                    chain = words[3]
                nextline = False
            elif line[0] == "\n":
                template_start = False
            elif target_start:
                target_lines.append(line)
            elif template_start:
                template_lines.append(line)

    target_seq = "".join(target_lines)
    target_seq = target_seq.replace("*","")
    template_seq = "".join(template_lines)
    template_seq = template_seq.replace("*","")
    return target_seq,template_seq

#Taken from https://code.activestate.com/recipes/577504/
def total_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}

    """
    dict_handler = lambda d: itertools.chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = getsizeof(o, default_size)

        if verbose:
            print(s, type(o), repr(o), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    return sizeof(o)

class Position:
    __slots__ = ['pos','wt_aa','mut_aas','pos_tags','mut_tags_map','stored','database_id','pdb_res_nr','checked','mappings',
                    'disorder_score','disorder_region']

    def __init__(self,pos = None,wt_aa='X',mut_aas = set(),tags = set(),pdb_res_nr = None,checked = False):
        self.pos = pos #the position in a protein sequence, int or None
        self.pdb_res_nr = pdb_res_nr #this is a string, due to insertion-codes
        self.wt_aa = wt_aa #wildtype amino acid type in one letter code
        self.mut_aas = mut_aas.copy() #set of mutant amino acid types in one letter code
        self.mut_tags_map = {}
        if len(mut_aas) == 0:
            self.pos_tags = tags.copy()
        else:
            self.pos_tags = set()
            for aa in mut_aas:
                 self.mut_tags_map[aa] = tags.copy()
        self.stored = False
        self.database_id = None
        self.checked = checked #Flag for sanity checks
        self.mappings = Mappings()
        self.disorder_score = None
        self.disorder_region = None

    def print_state(self):
        print(self.wt_aa,self.pos,self.mut_aas,self.mut_tags_map,self.pos_tags)
        return

    def check(self,wt_aa,overwrite=False):
        if self.wt_aa == 'X':
            self.wt_aa = wt_aa
            self.checked = True
            return True
        if self.wt_aa != wt_aa:
            if overwrite:
                self.wt_aa = wt_aa
                self.checked = True
                return True
            else:
                return False
        else:
            self.checked = True
            return True

    def fuse(self,position):
        warn = False
        if self.pos != position.pos:
            raise NameError('Cannot fuse positions with differing pos numbers')
            return
        if self.pdb_res_nr != position.pdb_res_nr:
            raise NameError('Cannot fuse positions with differing pdb_res_nr')
            return
        if self.wt_aa != position.wt_aa:
            print('Warning: fuse positions with different WT AAs:',self.pos,self.pdb_res_nr,self.wt_aa,position.wt_aa)
            warn = True
        self.mut_aas = self.mut_aas | position.mut_aas
        self.pos_tags = self.pos_tags | position.pos_tags
        for aa in position.mut_tags_map:
            if not aa in self.mut_tags_map:
                self.mut_tags_map[aa] = position.mut_tags_map[aa]
            else:
                self.mut_tags_map[aa] = self.mut_tags_map[aa] | position.mut_tags_map[aa]

        return warn

    def add_tags(self,tags):
        self.pos_tags = self.pos_tags | tags
        return

    def get_pos_tags(self):
        return self.pos_tags

    def getAACBase(self):
        return '%s%s' % (self.wt_aa,self.pos)

    def set_stored(self,value):
        self.stored = value
        return

    def get_stored(self):
        return self.stored

    def set_database_id(self,value):
        self.database_id = value
        return

    def get_database_id(self):
        return self.database_id

    def get_res_id(self):
        return self.pdb_res_nr

    def get_wt_aa(self):
        return self.wt_aa

    def get_mut_aas(self):
        return self.mut_aas

    def get_mut_tags(self,new_aa):
        return self.mut_tags_map[new_aa]

    def set_disorder_score(self,disorder_score):
        self.disorder_score = disorder_score
        return

    def get_disorder_score(self):
        return self.disorder_score

    def set_disorder_region(self,region):
        self.disorder_region = region
        return

    def get_disorder_region(self):
        return self.disorder_region

    def add_pos_res_mapping(self,pdb_id,chain,res_nr,mapping):
        self.mappings.add_mapping((pdb_id,chain,res_nr),mapping)

    def classify(self,config):
        self.mappings.weight_all(config,self.disorder_score,self.disorder_region)
        return

def is_mutant_ac(ac):
    if not ac.count('del') == 1 and not ac.count('ins') == 1:
        return False
    else:
        return True

class Indel:
    __slots__ = ['left_flank','left_flank_pos','right_flank','right_flank_pos','tags','stored','database_id','wt_prot','mut_prot','terminal']

    def __init__(self,left_flank = None,right_flank = None,tags = set()):
        self.terminal = False
        self.left_flank = left_flank
        if left_flank != None:
            self.left_flank_pos = int(left_flank[1:])
        else:
            self.left_flank_pos = 0
            self.terminal = True
        self.right_flank = right_flank
        if right_flank != None:
            self.right_flank_pos = int(right_flank[1:])
        else:
            self.right_flank_pos = None
            self.terminal = True
        self.tags = tags.copy()
        self.stored = False
        self.database_id = None
        self.wt_prot = None
        self.mut_prot = None

        return

    def set_proteins(self,wt_prot,mut_prot):
        self.wt_prot = wt_prot
        self.mut_prot = mut_prot
        return

    def set_database_id(self,database_id):
        self.database_id = database_id
        return

    def set_stored(self,boolean):
        self.stored = boolean
        return

class Deletion(Indel):

    def __init__(self,left_flank = None,right_flank = None,tags = set()):
        super().__init__(left_flank = left_flank,right_flank = right_flank,tags = tags)
        if self.left_flank_pos == 1:
            self.terminal = True
        return

    def create_protein_name(self,wt_name):
        return '%s_del_%s_%s' % (wt_name,self.left_flank[1:],self.right_flank[1:])

    def mutate_sequence(self,proteins):
        wt_seq = proteins[self.wt_prot].sequence
        if self.right_flank_pos == len(wt_seq):
            self.terminal = True
        mutated_sequence = '%s%s' % (wt_seq[:(self.left_flank_pos-1)],wt_seq[self.right_flank_pos:])
        proteins[self.mut_prot].sequence = mutated_sequence
        for (pos,aa) in enumerate(mutated_sequence):
            seq_pos = pos+1
            position = Position(pos=seq_pos,wt_aa=aa,checked=True)
            proteins[self.mut_prot].positions[seq_pos] = position
        return

    def get_notation(self):
        return '%s_%sdel' % (self.left_flank,self.right_flank)

    def get_scenario(self,n_wt,n_mut):
        if n_wt == 0 and n_mut == 0:
            return None
        if not self.terminal:
            if n_mut == 0:
                return ['FRA','InCl','CFRI']
            if n_wt == 0:
                return ['FRA']
            return ['FRA','InCl','CFRI','CFRA']
        else:
            if n_mut == 0:
                return ['FRC','InCl','CFRI']
            if n_wt == 0:
                return ['FRC']
            return ['FRC','InCl','CFRI','CFRA']

    def inclusion_check(self,wt_structure_annotations,mut_structure_annotations):
        del_list = []
        for struct_tuple in wt_structure_annotations:
            #The wildtype structure should contain both flanks
            if wt_structure_annotations[struct_tuple].is_covered(self.left_flank_pos) and wt_structure_annotations[struct_tuple].is_covered(self.right_flank_pos):
                continue #In this case do not delete annotation (suboptimal syntax to highlight this case)
            else:
                del_list.append(struct_tuple)
        for struct_tuple in del_list:
            del wt_structure_annotations[struct_tuple]

        #The flank positions are different for the mutant sequence!
        mut_left_flank_pos = self.left_flank_pos - 1
        mut_right_flank_pos = self.left_flank_pos

        del_list = []
        for struct_tuple in mut_structure_annotations:
            #The mutant structure must not include the WT flanks, but must include the mut flanks and they have to be adjacent, attention for terminal deletions here
            if self.terminal:
                if mut_right_flank_pos == 1: #Terminal deletion on the left side
                    if mut_structure_annotations[struct_tuple].is_covered(mut_right_flank_pos):
                        if not mut_structure_annotations[struct_tuple].is_terminal(mut_right_flank_pos):
                            del_list.append(struct_tuple)
                    else:
                        del_list.append(struct_tuple)
                else: #Terminal deletion on the right side
                    if mut_structure_annotations[struct_tuple].is_covered(mut_left_flank_pos):
                        if not mut_structure_annotations[struct_tuple].is_terminal(mut_left_flank_pos):
                            del_list.append(struct_tuple)
                    else:
                        del_list.append(struct_tuple)

            else:
                if mut_structure_annotations[struct_tuple].is_covered(mut_left_flank_pos) and mut_structure_annotations[struct_tuple].is_covered(mut_right_flank_pos):
                    if not mut_structure_annotations[struct_tuple].adjacency_check(mut_left_flank_pos,mut_right_flank_pos):
                        continue #In this case do not delete annotation
                    else:
                        del_list.append(struct_tuple)
                else:
                    del_list.append(struct_tuple)

        for struct_tuple in del_list:
            del mut_structure_annotations[struct_tuple]

        return wt_structure_annotations,mut_structure_annotations

class Insertion(Indel):
    __slots__ = ['inserted_sequence']

    def __init__(self,left_flank = None,right_flank = None,inserted_sequence = None,tags = set()):
        super().__init__(left_flank = left_flank,right_flank = right_flank,tags = tags)
        self.inserted_sequence = inserted_sequence
        if self.right_flank == None:
            self.right_flank_pos = 1
            self.left_flank_pos = 0
        if self.left_flank == None:
            self.left_flank_pos = self.right_flank_pos
        return

    def create_protein_name(self,wt_name,seq_name_length = 5):
        if self.left_flank != None:
            lfs = self.left_flank[1:]
        else:
            lfs = self.right_flank[1:]

        return '%s_ins_%s_%s' % (wt_name,lfs,self.inserted_sequence[:seq_name_length])

    def create_other_protein_name(self,wt_name,indel_protein_name):
        added_seq_len = len(indel_protein_name.split('_')[:-1])
        new_name = self.create_protein_name(self,wt_name,seq_name_length = added_seq_len + 1)
        if new_name == indel_protein_name:
            return None
        return new_name

    def mutate_sequence(self,proteins):
        wt_seq = proteins[self.wt_prot].sequence
        mutated_sequence = '%s%s%s' % (wt_seq[:(self.left_flank_pos)],self.inserted_sequence,wt_seq[(self.right_flank_pos - 1):])
        proteins[self.mut_prot].sequence = mutated_sequence
        for (pos,aa) in enumerate(mutated_sequence):
            seq_pos = pos+1
            position = Position(pos=seq_pos,wt_aa=aa,checked=True)
            proteins[self.mut_prot].positions[seq_pos] = position
        return

    def get_notation(self):
        if self.left_flank != None:
            lfs = self.left_flank
        else:
            lfs = ''
        if self.right_flank != None:
            if lfs != '':
                rfs = '_%s' % self.right_flank
            else:
                rfs = self.right_flank
        else:
            rfs = ''
        return '%s%sins%s' % (lfs,rfs,self.inserted_sequence)

    def get_scenario(self,n_wt,n_mut):
        if n_wt == 0 and n_mut == 0:
            return None
        if not self.terminal:
            if n_mut == 0:
                return ['FRA']
            if n_wt == 0:
                return ['FRA','InCl','CFRI']
            return ['FRA','InCl','CFRI','CFRA']
        else:
            if n_mut == 0:
                return ['FRC']
            if n_wt == 0:
                return ['FRC','InCl','CFRI']
            return ['FRC','InCl','CFRI','CFRA']

    def inclusion_check(self,wt_structure_annotations,mut_structure_annotations):
        del_list = []
        for struct_tuple in wt_structure_annotations:
            if self.terminal:
                if self.left_flank == None: #Terminal deletion on the left side
                    if wt_structure_annotations[struct_tuple].is_covered(self.right_flank_pos):
                        if not wt_structure_annotations[struct_tuple].is_terminal(self.right_flank_pos):
                            del_list.append(struct_tuple)
                    else:
                        del_list.append(struct_tuple)
                else: #Terminal insertion on the right side
                    if wt_structure_annotations[struct_tuple].is_covered(self.left_flank_pos):
                        if not wt_structure_annotations[struct_tuple].is_terminal(self.left_flank_pos):
                            del_list.append(struct_tuple)
                    else:
                        del_list.append(struct_tuple)
            else:
                #The wildtype structure should contain both flanks
                if wt_structure_annotations[struct_tuple].is_covered(self.left_flank_pos) and wt_structure_annotations[struct_tuple].is_covered(self.right_flank_pos):
                    #And they should be adjacent
                    if not wt_structure_annotations[struct_tuple].adjacency_check(self.left_flank_pos,self.right_flank_pos):
                        del_list.append(struct_tuple)
                else:
                    del_list.append(struct_tuple)

        for struct_tuple in del_list:
            del wt_structure_annotations[struct_tuple]

        #The flank positions are different for the mutant sequence! We take here the inner sides of the flanks.
        mut_left_flank_pos = self.left_flank_pos + 1
        mut_right_flank_pos = self.left_flank_pos + len(self.inserted_sequence) -1 #They can be identical if the only AA is inserted
        del_list = []

        for struct_tuple in mut_structure_annotations:
            #The mutant structure must include the Mut flanks, attention for terminal insertions here
            if self.terminal:
                if self.left_flank == None: #Terminal insertion on the left side
                    if not mut_structure_annotations[struct_tuple].is_covered(mut_right_flank_pos):
                        del_list.append(struct_tuple)

                else: #Terminal insertion on the right side
                    if not mut_structure_annotations[struct_tuple].is_covered(mut_left_flank_pos):
                        del_list.append(struct_tuple)

            else:
                if not (mut_structure_annotations[struct_tuple].is_covered(mut_left_flank_pos) and mut_structure_annotations[struct_tuple].is_covered(mut_right_flank_pos)):
                    del_list.append(struct_tuple)

        for struct_tuple in del_list:
            del mut_structure_annotations[struct_tuple]

        return wt_structure_annotations,mut_structure_annotations

class Substitution(Indel):
    __slots__ = ['inserted_sequence']

    def __init__(self,left_flank = None,right_flank = None,inserted_sequence = None,tags = set()):
        super().__init__(left_flank = left_flank,right_flank = right_flank,tags = tags)
        self.inserted_sequence = inserted_sequence
        if self.left_flank_pos == 1:
            self.terminal = True
        return

    def create_protein_name(self,wt_name,seq_name_length = 5):

        return '%s_delins_%s_%s_%s' % (wt_name,self.left_flank[1:],self.right_flank[1:],self.inserted_sequence[:seq_name_length])

    def create_other_protein_name(self,wt_name,indel_protein_name):
        added_seq_len = len(indel_protein_name.split('_')[:-1])
        new_name = self.create_protein_name(self,wt_name,seq_name_length = added_seq_len + 1)
        if new_name == indel_protein_name:
            return None
        return new_name

    def mutate_sequence(self,proteins):
        wt_seq = proteins[self.wt_prot].sequence
        if self.right_flank_pos == len(wt_seq):
            self.terminal = True
        mutated_sequence = '%s%s%s' % (wt_seq[:(self.left_flank_pos)],self.inserted_sequence,wt_seq[(self.right_flank_pos - 1):])
        proteins[self.mut_prot].sequence = mutated_sequence
        for (pos,aa) in enumerate(mutated_sequence):
            seq_pos = pos+1
            position = Position(pos=seq_pos,wt_aa=aa,checked=True)
            proteins[self.mut_prot].positions[seq_pos] = position
        return

    def get_notation(self):
        return '%s_%sdelins%s' % (self.left_flank,self.right_flank,self.inserted_sequence)

    def get_scenario(self,n_wt,n_mut):
        if n_wt == 0 and n_mut == 0:
            return None
        if not self.terminal:
            if n_mut == 0:
                return ['FRA','InCl','CFRI']
            if n_wt == 0:
                return ['FRA','InCl','CFRI']
            return ['FRA','InCl','CFRI','CFRA']
        else:
            if n_mut == 0:
                return ['FRC','InCl','CFRI']
            if n_wt == 0:
                return ['FRC','InCl','CFRI']
            return ['FRC','InCl','CFRI','CFRA']

    def inclusion_check(self,wt_structure_annotations,mut_structure_annotations):
        #This just checks if there is a region resolved in the structure of appropriate size, not if it is correct sequence. This should be checked elsewhere with the help of sequence identies
        del_list = []
        for struct_tuple in wt_structure_annotations:
            #The wildtype structure should contain both flanks (terminal or not), similar to deletions
            if not (wt_structure_annotations[struct_tuple].is_covered(self.left_flank_pos) and wt_structure_annotations[struct_tuple].is_covered(self.right_flank_pos)):
                del_list.append(struct_tuple)

        for struct_tuple in del_list:
            del wt_structure_annotations[struct_tuple]

        #The flank positions are different for the mutant sequence! We take here the inner sides of the (mutant) flanks.
        mut_left_flank_pos = self.left_flank_pos
        mut_right_flank_pos = self.left_flank_pos + len(self.inserted_sequence) - 1
        del_list = []

        for struct_tuple in wt_structure_annotations:
            #The mutant structure is here to be treated as the wt structure
            if not (mut_structure_annotations[struct_tuple].is_covered(mut_left_flank_pos) and mutt_structure_annotations[struct_tuple].is_covered(mut_right_flank_pos)):
                del_list.append(struct_tuple)

        for struct_tuple in del_list:
            del mut_structure_annotations[struct_tuple]

        return wt_structure_annotations,mut_structure_annotations

class Protein:
    __slots__ = ['u_ac','u_id','ref_ids','other_ids','pdb_id','positions','res_id_map','sequence','stored','completely_stored',
                    'go_terms','pathways','disorder_regions','disorder_tool','database_id','structure_annotations']

    def __init__(self,u_ac = None,u_id = None,ref_ids = set([]),pdb_id = None,positions = {},database_id = None, other_ids = []):
        self.u_ac = u_ac # UNIPROT accession number
        self.u_id = u_id # UNIPROT ID
        self.ref_ids = ref_ids.copy() #RefSeq IDs
        self.pdb_id = pdb_id
        self.positions = {}
        self.res_id_map = {}
        self.sequence = None
        self.stored = False #If the protein is stored in the database
        self.completely_stored = False #Gets true if all corresponding input positions are stored
        self.database_id = database_id #ID in the corresponding database table
        if database_id != None:
            self.stored = True
        self.go_terms = {}
        self.pathways = {}
        self.disorder_regions = None
        self.disorder_tool = None
        self.structure_annotations = {}

        if pdb_id == None:
            self.add_positions(positions)
        else:
            self.add_residues(positions)
        self.other_ids = {}
        for id_id,other_id in other_ids:
            self.other_ids[id_id] = other_id
        return

    def print_state(self):
        print('----State of %s----' % self.u_ac)
        print('Uniprot Id:',self.u_id)
        print('RefSeq Ids:',self.ref_ids)
        for pos in self.positions:
            self.positions[pos].print_state()

    def get_u_ac(self):
        if self.u_ac != None:
            return self.u_ac
        else:
            return self.pdb_id

    def add_positions(self,positions):
       for position in positions:
            warn = False
            if not position.pos in self.positions:
                self.positions[position.pos] = position
            else:
                warn = self.positions[position.pos].fuse(position)
                if warn:
                    print('Warning happened in fuse',self.u_ac,self.pdb_id)

    def add_residues(self,positions):
        for position in positions:
            warn = False
            if not position.pdb_res_nr in self.res_id_map:
                self.res_id_map[position.pdb_res_nr] = position
            else:
                warn = self.res_id_map[position.pdb_res_nr].fuse(position)
                if warn:
                    print('Warning happened in fuse',self.u_ac,self.pdb_id)

    def popNone(self):
        if self.pdb_id == None:
            if None in self.positions:
                tags = self.positions[None].pos_tags
                del self.positions[None]
                return True,tags
            else:
                return False,set()
        else:
            if None in self.res_id_map:
                tags = self.res_id_map[None].pos_tags
                del self.res_id_map[None]
                return True,tags
            else:
                return False,set()

    def getAACList(self):
        aaclist = {}
        for pos in self.positions:
            aac_base = self.positions[pos].getAACBase()
            aaclist[aac_base] = self.positions[pos].get_database_id()
        return aaclist

    def get_aac_base(self,pos):
        return self.positions[pos].getAACBase()

    def status(self):
        stat = '\n'.join([self.u_ac,self.pdb_id,str(len(self.positions)),self.sequence,self.stored,self.database_id,str(len(self.structure_annotations))])
        return stat

    def set_sequence(self,value):
        self.sequence = value
        return

    def get_sequence(self):
        return self.sequence

    def contains_position(self,pos):
        return pos in self.positions

    def set_stored(self,value):
        self.stored = value
        return

    def set_completely_stored(self):
        self.completely_stored = True
        return self.database_id

    def is_stored(self):
        return self.stored

    def add_ref_id(self,ref):
        self.ref_ids.add(ref)
        return

    def add_other_ids(self,id_id,other_id):
        self.other_ids[id_id] = other_id
        return

    def get_ref_ids(self):
        return self.ref_ids

    def get_u_id(self):
        return self.u_id

    def get_go_terms(self):
        return self.go_terms

    def get_pathways(self):
        return self.pathways

    def set_disorder_scores(self,disorder_scores):
        if disorder_scores == None:
            return
        for pos in disorder_scores:
            if not pos in self.positions:
                continue
            self.positions[pos].set_disorder_score(disorder_scores[pos])
        return

    def get_disorder_scores(self):
        disorder_scores = {}
        for pos in self.positions:
            disorder_scores[pos] = self.positions[pos].get_disorder_score()
        return disorder_scores

    def get_disorder_score(self,pos):
        return self.positions[pos].get_disorder_score()

    def set_disorder_regions(self,regions):
        self.disorder_regions = regions
        if regions == None:
            return
        for [a,b,region_type] in regions:
            for pos in range(int(a+1),int(b)):
                if not pos in self.positions:
                    continue
                self.positions[pos].set_disorder_region(region_type)
        return

    def get_disorder_regions(self):
        return self.disorder_regions

    def get_disorder_region(self,pos):
        return self.positions[pos].get_disorder_region()

    def set_disorder_tool(self,value):
        self.disorder_tool = value
        return

    def get_disorder_tool(self):
        return self.disorder_tool

    def is_position_stored(self,pos):
        if not self.contains_position(pos):
            return False
        else:
            return self.positions[pos].get_stored()

    def get_position(self,pos):
        if not pos in self.positions:
            return None
        return self.positions[pos]

    def get_position_ids(self):
        return self.positions.keys()

    def set_position_stored(self,pos,stored_value):
        self.positions[pos].set_stored(stored_value)
        return

    def set_position_database_id(self,pos,value):
        self.positions[pos].set_database_id(value)
        return

    def get_position_database_id(self,pos):
        return self.positions[pos].get_database_id()

    def get_pos_tags(self,pos):
        return self.positions[pos].get_pos_tags()

    def set_database_id(self,value):
        self.database_id = value
        return

    def get_database_id(self):
        return self.database_id

    def get_res_id(self,pos):
        return self.positions[pos].get_res_id()

    def get_wt_aa(self,pos):
        return self.positions[pos].get_wt_aa()

    def get_mut_aas(self,pos):
        return self.positions[pos].get_mut_aas()

    def get_mut_tags(self,pos,new_aa):
        if not pos in self.positions:
            return None
        return self.positions[pos].get_mut_tags(new_aa)

    def add_annotation(self,pdb_id,chain,anno_obj):
        self.structure_annotations[(pdb_id,chain)] = anno_obj
        return

    def remove_annotation(self,pdb_id,chain):
        if not (pdb_id,chain) in self.structure_annotations:
            return
        else:
            del self.structure_annotations[(pdb_id,chain)]
            return

    def remove_annotations(self):
        del self.structure_annotations
        return

    def get_annotation_list(self):
        return self.structure_annotations.keys()

    def get_structure_annotations(self):
        return self.structure_annotations

    def is_annotation_stored(self,pdb_id,chain):
        return self.structure_annotations[(pdb_id,chain)].get_stored()

    def set_alignment(self,pdb_id,chain,value):
        self.structure_annotations[(pdb_id,chain)].set_alignment(value)
        return

    def get_alignment(self,pdb_id,chain):
        return self.structure_annotations[(pdb_id,chain)].get_alignment()

    def set_coverage(self,pdb_id,chain,value):
        self.structure_annotations[(pdb_id,chain)].set_coverage(value)
        return

    def get_coverage(self,pdb_id,chain):
        return self.structure_annotations[(pdb_id,chain)].get_coverage()

    def set_sequence_id(self,pdb_id,chain,value):
        self.structure_annotations[(pdb_id,chain)].set_sequence_id(value)
        return

    def get_sequence_id(self,pdb_id,chain):
        return self.structure_annotations[(pdb_id,chain)].get_sequence_id()

    def set_annotation_db_id(self,pdb_id,chain,value):
        self.structure_annotations[(pdb_id,chain)].set_database_id(value)
        return

    def set_sub_infos(self,pdb_id,chain,value):
        self.structure_annotations[(pdb_id,chain)].set_sub_infos(value)
        return

    def get_sub_infos(self,pdb_id,chain):
        return self.structure_annotations[(pdb_id,chain)].get_sub_infos()

    def add_pos_res_mapping(self,pos,pdb_id,chain,res_nr,mapping):
        self.positions[pos].add_pos_res_mapping(pdb_id,chain,res_nr,mapping)
        return

    def classify(self,pos,config):
        self.positions[pos].classify(config)
        return

class Proteins:

    __slots__ = ['protein_map','stored_ids','stored_ids_mutant_excluded','completely_stored_ids','not_stored_ids','id_map','structures','complexes','indels','lite']

    def __init__(self,proteins,indels, lite = False):
        self.protein_map = proteins
        self.stored_ids = set()
        self.stored_ids_mutant_excluded = set()
        self.completely_stored_ids = set()
        self.not_stored_ids = set()
        self.id_map = {}
        self.structures = {}
        self.complexes = {}
        self.indels = {}
        for indel in indels:
            self.indels[indel.get_notation()] = indel
        self.lite = lite
        return

    def semi_deconstruct(self):
        del self.stored_ids
        del self.stored_ids_mutant_excluded
        del self.not_stored_ids
        del self.id_map
        del self.structures
        del self.complexes
        return

    def remove_protein_annotations(self,u_ac):
        self.protein_map[u_ac].remove_annotations()
        return

    def print_protein_state(self):
        for key in self.protein_map:
            self.protein_map[key].print_state()
        return

    def __getitem__(self,key):
        return self.protein_map[key]

    def __setitem__(self,key,prot):
        self.protein_map[key] = prot
        if prot.is_stored():
            self.stored_ids.add(prot.get_database_id())
            self.id_map[prot.get_database_id()] = key
        return

    def __delitem__(self,key):
        del self.protein_map[key]
        return

    def get_protein_map(self):
        return self.protein_map

    def get_protein_u_acs(self):
        return self.protein_map.keys()

    def get_position(self,u_ac,pos):
        return self.protein_map[u_ac].get_position(pos)

    def get_complexes(self):
        return self.complexes

    def get_complex_list(self):
        return self.complexes.keys()

    def add_complex(self,pdb_id,complex_obj):
        self.complexes[pdb_id] = complex_obj
        return

    def contains_complex(self,pdb_id):
        return pdb_id in self.complexes

    def is_complex_stored(self,pdb_id):
        return self.complexes[pdb_id].get_stored()

    def set_complex_db_id(self,pdb_id,value):
        self.complexes[pdb_id].set_database_id(value)
        return

    def get_complex_chains(self,pdb_id):
        return self.complexes[pdb_id].get_chains()

    def set_lig_profile(self,pdb_id,value):
        self.complexes[pdb_id].set_lig_profile(value)
        return

    def get_lig_profile(self,pdb_id):
        return self.complexes[pdb_id].get_lig_profile()

    def get_lig_str(self,pdb_id):
        return self.complexes[pdb_id].getLigProfileStr()

    def set_ion_profile(self,pdb_id,value):
        self.complexes[pdb_id].set_ion_profile(value)
        return

    def get_ion_profile(self,pdb_id):
        return self.complexes[pdb_id].get_ion_profile()

    def get_ion_str(self,pdb_id):
        return self.complexes[pdb_id].getIonProfileStr()

    def set_metal_profile(self,pdb_id,value):
        self.complexes[pdb_id].set_metal_profile(value)
        return

    def get_metal_profile(self,pdb_id):
        return self.complexes[pdb_id].get_metal_profile()

    def get_metal_str(self,pdb_id):
        return self.complexes[pdb_id].getMetalProfileStr()

    def set_chain_chain_profile(self,pdb_id,value):
        self.complexes[pdb_id].set_chain_chain_profile(value)
        return

    def get_chain_chain_profile(self,pdb_id):
        return self.complexes[pdb_id].get_chain_chain_profile()

    def get_chain_chain_str(self,pdb_id):
        return self.complexes[pdb_id].getChainChainProfileStr()

    def set_interaction_partners(self,pdb_id,value):
        self.complexes[pdb_id].set_interaction_partners(value)
        return

    def get_interaction_partners(self,pdb_id):
        return self.complexes[pdb_id].get_interaction_partners()

    def set_chain_type_map(self,pdb_id,value):
        self.complexes[pdb_id].set_chain_type_map(value)
        return

    def get_chains_str(self,pdb_id):
        return self.complexes[pdb_id].getChainStr()

    def get_homomers_str(self,pdb_id):
        return self.complexes[pdb_id].getHomomersStr()

    def get_complex_homomers(self,pdb_id,chain):
        return self.complexes[pdb_id].get_homomers(chain)

    def get_resolution(self,pdb_id):
        return self.complexes[pdb_id].get_resolution()

    def get_structures(self):
        return self.structures

    def add_structure(self,pdb_id,chain,struct_obj):
        self.structures[(pdb_id,chain)] = struct_obj
        return

    def contains_structure(self,pdb_id,chain):
        return (pdb_id,chain) in self.structures

    def set_oligo(self,pdb_id,chain,value):
        self.structures[(pdb_id,chain)].set_oligo(value)
        return

    def get_oligo(self,pdb_id,chain):
        return self.structures[(pdb_id,chain)].get_oligo()

    def set_structure_db_id(self,pdb_id,chain,value):
        self.structures[(pdb_id,chain)].set_database_id(value)
        return

    def get_structure_db_id(self,pdb_id,chain):
        return self.structures[(pdb_id,chain)].get_database_id()

    def set_structure_stored(self,pdb_id,chain,value):
        self.structures[(pdb_id,chain)].set_stored(value)
        return

    def is_structure_stored(self,pdb_id,chain):
        return self.structures[(pdb_id,chain)].get_stored()

    def add_residue(self,pdb_id,chain,res_nr,residue_obj):
        self.structures[(pdb_id,chain)].add_residue(res_nr,residue_obj)
        return

    def set_residue_db_id(self,pdb_id,chain,res_nr,value):
        self.structures[(pdb_id,chain)].set_residue_db_id(res_nr,value)
        return

    def contains_residue(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].contains_residue(res_nr)

    def get_residue_db_id(self,pdb_id,chain,res_nr):
        if not self.contains_residue(pdb_id,chain,res_nr):
            return None
        return self.structures[(pdb_id,chain)].get_residue_db_id(res_nr)

    def get_residue_aa(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_aa(res_nr)

    def get_residue_sld(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_sld(res_nr)

    def get_residue_scd(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_scd(res_nr)

    def get_residue_homomer_dists(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_homomer_dists(res_nr)

    def get_residue_centralities(self,pdb_id,chain,res_nr,get_whats_there=False):
        return self.structures[(pdb_id,chain)].get_residue_centralities(res_nr,get_whats_there=get_whats_there)

    def get_residue_modres(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_modres(res_nr)

    def get_residue_b_factor(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_b_factor(res_nr)

    def get_residue_rsa(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_rsa(res_nr)

    def get_residue_ssa(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_ssa(res_nr)

    def get_residue_phi(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_phi(res_nr)

    def get_residue_psi(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_psi(res_nr)

    def get_residue_link_information(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_link_information(res_nr)

    def get_residue_interaction_profile(self,pdb_id,chain,res_nr,get_whats_there = False):
        return self.structures[(pdb_id,chain)].get_residue_interaction_profile(res_nr,get_whats_there = get_whats_there)

    def get_residue_interaction_profile_str(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_interaction_profile_str(res_nr)

    def get_residue_milieu(self,pdb_id,chain,res_nr):
        return self.structures[(pdb_id,chain)].get_residue_milieu(res_nr)

    def add_residue_classification(self,pdb_id,chain,res_nr,Class,simpleClass):
        self.structures[(pdb_id,chain)].add_residue_classification(res_nr,Class,simpleClass)
        return

    def get_protein_annotation_list(self,u_ac):
        return self.protein_map[u_ac].get_annotation_list()

    def get_protein_structure_annotations(self,u_ac):
        return self.protein_map[u_ac].get_structure_annotations()

    def add_mapping_to_structure(self,pdb_id,chain,u_ac):
        self.structures[(pdb_id,chain)].add_mapping(u_ac)
        return

    def is_annotation_stored(self,pdb_id,chain,u_ac):
        return self.protein_map[u_ac].is_annotation_stored(pdb_id,chain)

    def set_alignment(self,u_ac,pdb_id,chain,alignment):
        self.protein_map[u_ac].set_alignment(pdb_id,chain,alignment)
        return

    def set_alignment_by_db_id(self,prot_id,pdb_id,chain,value):
        self.getByDbId(prot_id).set_alignment(pdb_id,chain,value)
        return

    def get_alignment(self,u_ac,pdb_id,chain):
        return self.protein_map[u_ac].get_alignment(pdb_id,chain)

    def set_coverage_by_db_id(self,prot_id,pdb_id,chain,value):
        self.getByDbId(prot_id).set_coverage(pdb_id,chain,value)
        return

    def set_coverage(self,u_ac,pdb_id,chain,value):
        self.protein_map[u_ac].set_coverage(pdb_id,chain,value)
        return

    def get_coverage(self,u_ac,pdb_id,chain):
        return self.protein_map[u_ac].get_coverage(pdb_id,chain)

    def set_sequence_id_by_db_id(self,prot_id,pdb_id,chain,value):
        self.getByDbId(prot_id).set_sequence_id(pdb_id,chain,value)
        return

    def set_sequence_id(self,u_ac,pdb_id,chain,value):
        self.protein_map[u_ac].set_sequence_id(pdb_id,chain,value)
        return

    def get_sequence_id(self,u_ac,pdb_id,chain):
        return self.protein_map[u_ac].get_sequence_id(pdb_id,chain)

    def set_annotation_db_id_by_db_id(self,prot_id,pdb_id,chain,value):
        self.getByDbId(prot_id).set_annotation_db_id(pdb_id,chain,value)
        return

    def set_sub_infos(self,u_ac,pdb_id,chain,value):
        self.protein_map[u_ac].set_sub_infos(pdb_id,chain,value)
        return

    def get_sub_infos(self,u_ac,pdb_id,chain):
        return self.protein_map[u_ac].get_sub_infos(pdb_id,chain)

    def set_not_stored_ids(self,value):
        self.not_stored_ids = value
        return

    def get_not_stored_ids(self):
        return self.not_stored_ids

    def get_not_stored_acs(self):
        if self.lite:
            return self.get_protein_u_acs()
        acs = []
        for db_id in self.not_stored_ids:
            acs.append(self.id_map[db_id])
        return acs

    def set_stored_ids(self,stored_ids,stored_ids_mutant_excluded):
        self.stored_ids = stored_ids
        self.stored_ids_mutant_excluded = stored_ids_mutant_excluded
        return

    def get_stored_ids(self,exclude_indel_mutants = False,exclude_completely_stored=False):
        if not exclude_indel_mutants:
            if exclude_completely_stored:
                return self.stored_ids - self.completely_stored_ids
            return self.stored_ids
        else:
            if exclude_completely_stored:
                return self.stored_ids_mutant_excluded - self.completely_stored_ids
            return self.stored_ids_mutant_excluded

    def contains(self,u_ac):
        return u_ac in self.protein_map

    def isEmpty(self):
        return len(self.protein_map) == 0

    def is_protein_stored(self,u_ac):
        return self.protein_map[u_ac].is_stored()

    def set_completely_stored(self,u_ac):
        prot_id = self.protein_map[u_ac].set_completely_stored()
        self.completely_stored_ids.add(prot_id)
        return

    def is_completely_stored(self,u_ac):
        return self.protein_map[u_ac].completely_stored

    def isStored(self,prot_id):
        return prot_id in self.stored_ids

    def generate_id_map(self):
        for u_ac in self.protein_map:
            self.id_map[self.protein_map[u_ac].get_database_id()] = u_ac
        return

    def getByDbId(self,database_id):
        return self.protein_map[self.id_map[database_id]]

    def set_protein_stored(self,u_ac,value):
        self.protein_map[u_ac].set_stored(value)
        return

    def set_protein_db_id(self,u_ac,value):
        self.protein_map[u_ac].set_database_id(value)
        return

    def get_protein_db_id(self,u_ac):
        return self.protein_map[u_ac].get_database_id()

    def set_protein_sequence(self,u_ac,value):
        self.protein_map[u_ac].set_sequence(value)
        return

    def get_sequence(self,u_ac):
        return self.protein_map[u_ac].get_sequence()

    def add_ref_id(self,u_ac,ref):
        self.protein_map[u_ac].add_ref_id(ref)
        return

    def get_ref_ids(self,u_ac):
        return self.protein_map[u_ac].get_ref_ids()

    def get_u_id(self,u_ac):
        return self.protein_map[u_ac].get_u_id()

    def get_go_terms(self,u_ac):
        return self.protein_map[u_ac].get_go_terms()

    def get_pathways(self,u_ac):
        return self.protein_map[u_ac].get_pathways()

    def set_disorder_scores(self,u_ac,value):
        if value == None:
            return
        self.protein_map[u_ac].set_disorder_scores(value)
        return

    def get_disorder_scores(self,u_ac):
        return self.protein_map[u_ac].get_disorder_scores()

    def get_disorder_score(self,u_ac,pos):
        return self.protein_map[u_ac].get_disorder_score(pos)

    def set_disorder_regions(self,u_ac,value):
        self.protein_map[u_ac].set_disorder_regions(value)
        return

    def get_disorder_regions(self,u_ac):
        return self.protein_map[u_ac].get_disorder_regions()

    def get_disorder_region(self,u_ac,pos):
        return self.protein_map[u_ac].get_disorder_region(pos)

    def set_disorder_tool(self,u_ac,value):
        self.protein_map[u_ac].set_disorder_tool(value)
        return

    def get_disorder_tool(self,u_ac):
        return self.protein_map[u_ac].get_disorder_tool()

    def position_in_protein_by_db_id(self,database_id,pos):
        return self.getByDbId(database_id).contains_position(pos)

    def set_position_stored(self,database_id,pos,stored_value):
        self.getByDbId(database_id).set_position_stored(pos,stored_value)
        return

    def set_position_database_id(self,database_id,pos,value):
        self.getByDbId(database_id).set_position_database_id(pos,value)
        return

    def get_position_database_id(self,u_ac,pos):
        return self.protein_map[u_ac].get_position_database_id(pos)

    def get_pos_tags(self,u_ac,pos):
        return self.protein_map[u_ac].get_pos_tags(pos)

    def get_wt_aa(self,u_ac,pos):
        return self.protein_map[u_ac].get_wt_aa(pos)

    def get_mut_aas(self,u_ac,pos):
        return self.protein_map[u_ac].get_mut_aas(pos)

    def get_mut_tags(self,u_ac,pos,new_aa):
        if new_aa == None:
            return None
        return self.protein_map[u_ac].get_mut_tags(pos,new_aa)

    def get_structure_list(self):
        return set(self.structures.keys())

    def getStoredStructureIds(self):
        stored_ids = {}
        for (pdb_id,chain) in self.structures:
            if self.structures[(pdb_id,chain)].get_stored():
                stored_ids[self.structures[(pdb_id,chain)].get_database_id()] = (pdb_id,chain)
        return stored_ids

    def get_protein_database_id(self,u_ac):
        return self.protein_map[u_ac].get_database_id()

    def get_position_ids(self,u_ac):
        return self.protein_map[u_ac].get_position_ids()

    def is_position_stored(self,u_ac,pos):
        return self.protein_map[u_ac].is_position_stored(pos)

    def get_aac_base(self,u_ac,pos):
        return self.protein_map[u_ac].get_aac_base(pos)

    def getAACList(self,u_ac):
        return self.protein_map[u_ac].getAACList()

    def get_res_id(self,u_ac,pos):
        return self.protein_map[u_ac].get_res_id(pos)

    def add_annotation(self,u_ac,pdb_id,chain,anno_obj):
        self.protein_map[u_ac].add_annotation(pdb_id,chain,anno_obj)
        return

    def remove_annotation(self,u_ac,pdb_id,chain):
        self.protein_map[u_ac].remove_annotation(pdb_id,chain)
        return

    def get_target_dict(self,pdb_id):
        target_dict = {}
        chains = self.get_complex_chains(pdb_id)
        for chain in chains:
            if chains[chain] != 'Protein':
                continue
            target_dict[chain] = {}
            if not (pdb_id,chain) in self.structures:
                continue
            mapped_proteins = self.structures[(pdb_id,chain)].get_mapped_proteins()
            for u_ac in mapped_proteins:
                sub_infos = self.get_sub_infos(u_ac,pdb_id,chain)
                for pos in sub_infos:
                    res_nr = sub_infos[pos][0]
                    if res_nr == None:
                        continue
                    if not res_nr in target_dict[chain]:
                        target_dict[chain][res_nr] = []
                    target_dict[chain][res_nr].append((u_ac,pos))
        return target_dict

    def add_pos_res_mapping(self,u_ac,pos,pdb_id,chain,res_nr,mapping):
        self.protein_map[u_ac].add_pos_res_mapping(pos,pdb_id,chain,res_nr,mapping)

    def classify(self,u_ac,pos,config):
        self.protein_map[u_ac].classify(pos,config)
        return

class Structure_annotation:
    __slots__ = ['u_ac','pdb_id','chain','alignment','coverage','sequence_identity','sub_infos','stored','database_id']

    def __init__(self,u_ac,pdb_id,chain,alignment = None,stored = False):
        self.u_ac = u_ac
        self.pdb_id = pdb_id
        self.chain = chain
        self.alignment = alignment
        self.coverage = None
        self.sequence_identity = None
        self.sub_infos = {} #{pos:(res_nr,res_aa,structure_sequence_number)}
        self.stored = stored
        self.database_id = None

    def set_alignment(self,value):
        self.alignment = value
        return

    def get_alignment(self):
        return self.alignment

    def set_coverage(self,value):
        self.coverage = value
        return

    def get_coverage(self):
        return self.coverage

    def set_sequence_id(self,value):
        self.sequence_identity = value
        return

    def get_sequence_id(self):
        return self.sequence_identity

    def set_database_id(self,value):
        self.database_id = value
        return

    def get_stored(self):
        return self.stored

    def set_sub_infos(self,value):
        self.sub_infos = value
        return

    def get_sub_infos(self):
        return self.sub_infos

    def is_covered(self,pos):
        if not pos in self.sub_infos:
            return False
        if self.sub_infos[pos][0] == None:
            return False
        return True

    def adjacency_check(self,left_pos,right_pos):
        return ((self.sub_infos[right_pos] - self.sub_infos[left_pos][2]) == 1)

    def is_terminal(self,pos):
        if not pos in self.sub_infos:
            return False
        structure_sequence_number = self.sub_infos[pos][2]
        if structure_sequence_number == 1:
            return True
        target_seq,template_seq = processAlignmentData(self.alignment)
        template_seq = template_seq.replace('-','')
        if structure_sequence_number == len(template_seq):
            return True
        return False

class Complex:
    __slots__ = ['pdb_id','resolution','interaction_partners','database_id','chains','stored','lig_profile','metal_profile','ion_profile','chain_chain_profile','homomers']

    def __init__(self,pdb_id,resolution = None,chains_str = None,lig_profile = None, lig_profile_str = None, metal_profile = None,
                    metal_profile_str = None, ion_profile = None, ion_profile_str = None, chain_chain_profile = None,
                    chain_chain_profile_str = None, stored = False, database_id = None, homomers = {}, homomers_str = None):
        self.pdb_id = pdb_id
        self.chains = {} #{One-letter-chain-id:chaintype} ; possible chaintypes: 'Protein','DNA','RNA,'Peptide'
        self.resolution = resolution
        self.database_id = database_id
        self.interaction_partners = []
        self.stored = stored
        self.lig_profile = lig_profile
        self.metal_profile = metal_profile
        self.ion_profile = ion_profile
        self.chain_chain_profile = chain_chain_profile
        self.homomers = homomers #Output of pdbParser.getInfo

        if homomers_str != None:
            self.parseHomomersStr(homomers_str)

        if chains_str != None:
            self.parseChainStr(chains_str)

        if lig_profile_str != None:
            self.lig_profile = self.parseProfileStr(lig_profile_str)

        if metal_profile_str != None:
            self.metal_profile = self.parseProfileStr(metal_profile_str)

        if ion_profile_str != None:
            self.ion_profile = self.parseProfileStr(ion_profile_str)

        if chain_chain_profile_str != None:
            self.chain_chain_profile = self.parseProfileStr(chain_chain_profile_str)

    def parseHomomersStr(self,homomers_str):
        homomers = {}
        for chains in homomers_str.split(','):
            for chain in chains:
                homomers[chain] = []
                for chain2 in chains:
                    homomers[chain].append(chain2)
        self.homomers = homomers

    def get_homomers(self,chain):
        if not chain in self.homomers:
            return chain
        return self.homomers[chain]

    def getHomomersStr(self):
        used_chains = set()
        str_parts = []
        for chain in self.homomers:
            if chain in used_chains:
                continue
            str_part = chain
            used_chains.add(chain)
            for chain2 in self.homomers[chain]:
                if chain2 == chain:
                    continue
                used_chains.add(chain2)
                str_part += chain2
            str_parts.append(str_part)
        return ','.join(str_parts)

    def getChainStr(self):
        chain_str_parts = []
        for chain in self.chains:
            chain_str_parts.append('%s:%s' % (chain,self.chains[chain]))
        return ','.join(chain_str_parts)

    def parseChainStr(self,chain_str):
        for chain_str_part in chain_str.split(','):
            if chain_str_part.count(':') < 1:
                continue
            chain,chaintype = chain_str_part.split(':')
            self.chains[chain] = chaintype
        return

    def parseProfileStr(self,profile_str):
        profile = {}
        if profile_str == '':
            return profile
        parts = profile_str.split(',')
        for part in parts:
            id_part,value_part = part.split(':')
            id_1,id_2 = id_part.split('_')
            deg,score = value_part.split('_')
            deg = int(deg)
            score = float(score)
            profile[(id_1,id_2)] = (deg,score)
        return profile

    def getProfileStr(self,profile):
        if profile == None:
            return ''
        profile_str_parts = []
        for id_1,id_2 in profile:
            deg,score = profile[(id_1,id_2)]
            profile_str_parts.append('%s_%s:%s_%s' % (id_1,id_2,str(deg),str(score)))
        return ','.join(profile_str_parts)

    def set_lig_profile(self,value):
        self.lig_profile = value
        return

    def get_lig_profile(self):
        return self.lig_profile

    def getLigProfileStr(self):
        return self.getProfileStr(self.lig_profile)

    def set_metal_profile(self,value):
        self.metal_profile = value
        return

    def get_metal_profile(self):
        return self.metal_profile

    def getMetalProfileStr(self):
        return self.getProfileStr(self.metal_profile)

    def set_ion_profile(self,value):
        self.ion_profile = value
        return

    def get_ion_profile(self):
        return self.ion_profile

    def getIonProfileStr(self):
        return self.getProfileStr(self.ion_profile)

    def set_chain_chain_profile(self,value):
        self.chain_chain_profile = value
        return

    def get_chain_chain_profile(self):
        return self.chain_chain_profile

    def getChainChainProfileStr(self):
        return self.getProfileStr(self.chain_chain_profile)

    def set_chain_type_map(self,value):
        self.chains = value
        return

    def set_interaction_partners(self,value):
        self.interaction_partners = value
        return

    def get_interaction_partners(self):
        return self.interaction_partners

    def get_stored(self):
        return self.stored

    def get_chains(self):
        return self.chains

    def get_resolution(self):
        return self.resolution

    def set_database_id(self,value):
        self.database_id = value
        return

class Structure:
    __slots__ = ['pdb_id','chain','oligo','database_id','stored','mapped_proteins','residues']

    def __init__(self,pdb_id,chain, oligo = set(),mapped_proteins = [], database_id = None):
        self.pdb_id = pdb_id
        self.chain = chain
        self.database_id = database_id
        if isinstance(oligo,str):
            self.oligo = set(oligo)
        else:
            self.oligo = oligo.copy()
        self.stored = (database_id != None)
        self.mapped_proteins = mapped_proteins
        self.residues = {}

    def status(self):
        stat = '\n'.join([self.pdb_id,self.chain,str(self.resolution),str(self.oligo),str(self.database_id),str(self.stored),str(len(self.residues))])
        return stat

    def add_mapping(self,u_ac):
        self.mapped_proteins.append(u_ac)
        return

    def get_mapped_proteins(self):
        return self.mapped_proteins

    def set_database_id(self,value):
        self.database_id = value
        return

    def get_database_id(self):
        return self.database_id

    def set_oligo(self,value):
        self.oligo = value
        return

    def get_oligo(self):
        return self.oligo

    def set_stored(self,value):
        self.stored = value
        return

    def get_stored(self):
        return self.stored

    def add_residue(self,res_nr,residue_obj):
        self.residues[res_nr] = residue_obj
        return

    def get_residue_list(self):
        return self.residues.keys()

    def set_residue_db_id(self,res_nr,value):
        self.residues[res_nr].set_database_id(value)
        return

    def contains_residue(self,res_nr):
        return res_nr in self.residues

    def get_residue_db_id(self,res_nr):
        return self.residues[res_nr].get_database_id()

    def get_residue_aa(self,res_nr):
        return self.residues[res_nr].get_aa()

    def get_residue_sld(self,res_nr):
        return self.residues[res_nr].get_ligand_distances()

    def get_residue_scd(self,res_nr):
        return self.residues[res_nr].get_chain_distances()

    def get_residue_homomer_dists(self,res_nr):
        return self.residues[res_nr].get_homomer_dists()

    def get_residue_centralities(self,res_nr,get_whats_there=False):
        return self.residues[res_nr].get_centralities(get_whats_there=get_whats_there)

    def get_residue_modres(self,res_nr):
        return self.residues[res_nr].get_modres()

    def get_residue_b_factor(self,res_nr):
        return self.residues[res_nr].get_b_factor()

    def get_residue_rsa(self,res_nr):
        return self.residues[res_nr].get_rsa()

    def get_residue_ssa(self,res_nr):
        return self.residues[res_nr].get_ssa()

    def get_residue_phi(self,res_nr):
        return self.residues[res_nr].get_phi()

    def get_residue_psi(self,res_nr):
        return self.residues[res_nr].get_psi()

    def get_residue_link_information(self,res_nr):
        return self.residues[res_nr].get_residue_link_information()

    def get_residue_interaction_profile(self,res_nr,get_whats_there = False):
        return self.residues[res_nr].get_interaction_profile(get_whats_there = get_whats_there)

    def get_residue_interaction_profile_str(self,res_nr):
        return self.residues[res_nr].get_interaction_profile_str()

    def get_residue_milieu(self,res_nr):
        return self.residues[res_nr].get_milieu()

    def add_residue_classification(self,res_nr,Class,simpleClass):
        self.residues[res_nr].set_classification(Class,simpleClass)
        return

class Mappings:
    __slots__ = ['qualities','covs','seq_ids','rsas','ssas','lig_dists','metal_dists','ion_dists','chain_dists','rna_dists','dna_dists','homo_dists',
                'profiles','weighted_profile','weighted_profile_str','centralities','weighted_centralities','weighted_centralities_str','b_factors',
                'weighted_b_factor','modres','weighted_modres','weighted_ssa',
                'phis','weighted_phi','psis','weighted_psi','intra_ssbonds','weighted_intra_ssbond','weighted_inter_ssbond', 'ssbond_lengths',
                'weighted_ssbond_length', 'intra_links', 'weighted_intra_link', 'weighted_inter_link', 'link_lengths', 'weighted_link_length',
                'cis_conformations', 'weighted_cis_conformation', 'cis_followers', 'weighted_cis_follower',
                'inter_chain_median_kds', 'weighted_inter_chain_median_kd', 'inter_chain_dist_weighted_kds', 'weighted_inter_chain_dist_weighted_kd',
                'inter_chain_median_rsas', 'weighted_inter_chain_median_rsa', 'inter_chain_dist_weighted_rsas',
                'weighted_inter_chain_dist_weighted_rsa', 'intra_chain_median_kds', 'weighted_intra_chain_median_kd',
                'intra_chain_dist_weighted_kds', 'weighted_intra_chain_dist_weighted_kd', 'intra_chain_median_rsas',
                'weighted_intra_chain_median_rsa', 'intra_chain_dist_weighted_rsas', 'weighted_intra_chain_dist_weighted_rsa',
                'weighted_lig_dist','lig_dist_conf','weighted_metal_dist','metal_dist_conf','weighted_ion_dist','ion_dist_conf',
                'weighted_chain_dist','chain_dist_conf','weighted_rna_dist','rna_dist_conf','weighted_dna_dist','dna_dist_conf',
                'weighted_homo_dist','homo_dist_conf','weighted_surface_value','weighted_location','max_cov','location_conf',
                'res_classes','aa_ids','max_seq_res','recommended_res','interaction_recommendations',
                'rin_class','rin_simple_class','Class','simple_class','classification_conf','max_cov_rsa','resolutions','res_aas']
    
    def __init__(self,raw_results = None):
        self.qualities = {}
        self.seq_ids = {}
        self.covs = {}
        self.rsas = {}
        self.ssas = {}
        self.lig_dists = {}
        self.metal_dists = {}
        self.ion_dists = {}
        self.chain_dists = {}
        self.rna_dists = {}
        self.dna_dists = {}
        self.homo_dists = {}
        self.profiles = {}
        self.centralities = {}
        self.phis = {}
        self.psis = {}
        self.intra_ssbonds = {}
        self.ssbond_lengths = {}
        self.intra_links = {}
        self.link_lengths = {}
        self.cis_conformations = {}
        self.cis_followers = {}
        self.inter_chain_median_kds = {}
        self.inter_chain_dist_weighted_kds = {}
        self.inter_chain_median_rsas = {}
        self.inter_chain_dist_weighted_rsas = {}
        self.intra_chain_median_kds = {}
        self.intra_chain_dist_weighted_kds = {}
        self.intra_chain_median_rsas = {}
        self.intra_chain_dist_weighted_rsas = {}
        self.b_factors = {}
        self.modres = {}
        self.res_classes = {}
        self.aa_ids = {}
        self.res_aas = {}
        self.resolutions = {}
        self.recommended_res = None
        self.max_seq_res = None
        self.weighted_surface_value = None
        self.weighted_location = None
        self.weighted_profile = None
        self.weighted_profile_str = None
        self.weighted_centralities = None
        self.weighted_centralities_str = None
        self.weighted_b_factor = None
        self.weighted_modres = None
        self.weighted_ssa = None
        self.weighted_phi = None
        self.weighted_psi = None
        self.weighted_intra_ssbond = None
        self.weighted_inter_ssbond = None
        self.weighted_ssbond_length = None
        self.weighted_intra_link = None
        self.weighted_inter_link = None
        self.weighted_link_length = None
        self.weighted_cis_conformation = None
        self.weighted_cis_follower = None
        self.weighted_inter_chain_median_kd = None
        self.weighted_inter_chain_dist_weighted_kd = None
        self.weighted_inter_chain_median_rsa = None
        self.weighted_inter_chain_dist_weighted_rsa = None
        self.weighted_intra_chain_median_kd = None
        self.weighted_intra_chain_dist_weighted_kd = None
        self.weighted_intra_chain_median_rsa = None
        self.weighted_intra_chain_dist_weighted_rsa = None
        self.weighted_lig_dist = None
        self.weighted_metal_dist = None
        self.weighted_ion_dist = None
        self.weighted_chain_dist = None
        self.weighted_rna_dist = None
        self.weighted_dna_dist = None
        self.weighted_homo_dist = None
        self.rin_class = None
        self.rin_simple_class = None
        self.Class = None
        self.simple_class = None
        self.interaction_recommendations = None
        self.lig_dist_conf = None
        self.metal_dist_conf = None
        self.ion_dist_conf = None
        self.chain_dist_conf = None
        self.rna_dist_conf = None
        self.dna_dist_conf = None
        self.homo_dist_conf = None
        self.location_conf = None
        self.classification_conf = None
        if raw_results != None:
            (self.recommended_res, self.max_seq_res, self.weighted_surface_value, self.weighted_location, self.weighted_profile_str,
            self.weighted_centralities_str, self.weighted_b_factor, self.weighted_modres, self.weighted_ssa, self.weighted_phi,
            self.weighted_psi, self.weighted_intra_ssbond, self.weighted_inter_ssbond, self.weighted_ssbond_length, self.weighted_intra_link,
            self.weighted_inter_link, self.weighted_link_length, self.weighted_cis_conformation, self.weighted_cis_follower,
            self.weighted_inter_chain_median_kd, self.weighted_inter_chain_dist_weighted_kd, self.weighted_inter_chain_median_rsa,
            self.weighted_inter_chain_dist_weighted_rsa, self.weighted_intra_chain_median_kd, self.weighted_intra_chain_dist_weighted_kd,
            self.weighted_intra_chain_median_rsa, self.weighted_intra_chain_dist_weighted_rsa, self.weighted_lig_dist,
            self.weighted_metal_dist, self.weighted_ion_dist, self.weighted_chain_dist, self.weighted_rna_dist, self.weighted_dna_dist,
            self.weighted_homo_dist, self.rin_class, self.rin_simple_class, self.Class, self.simple_class, self.interaction_recommendations,
            self.lig_dist_conf, self.metal_dist_conf, self.ion_dist_conf, self.chain_dist_conf, self.rna_dist_conf, self.dna_dist_conf,
            self.homo_dist_conf, self.location_conf, self.classification_conf) = raw_results
        return

    def add_mapping(self,mapping_id,mapping):
        (quality,seq_id,cov,rsa,ssa,lig_dist,metal_dist,ion_dist,chain_dist,rna_dist,dna_dist,homo_dist,profile,centralities,
        phi,psi,intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
        inter_chain_median_kd,inter_chain_dist_weighted_kd,
        inter_chain_median_rsa,inter_chain_dist_weighted_rsa,intra_chain_median_kd,
        intra_chain_dist_weighted_kd,intra_chain_median_rsa,intra_chain_dist_weighted_rsa,b_factor,modres,res_class,res_simple_class,
        identical_aa,resolution,res_aa) = mapping
        self.qualities[mapping_id] = quality
        self.seq_ids[mapping_id] = seq_id
        self.covs[mapping_id] = cov
        self.rsas[mapping_id] = rsa
        self.ssas[mapping_id] = ssa
        self.lig_dists[mapping_id] = lig_dist
        self.metal_dists[mapping_id] = metal_dist
        self.ion_dists[mapping_id] = ion_dist
        self.chain_dists[mapping_id] = chain_dist
        self.rna_dists[mapping_id] = rna_dist
        self.dna_dists[mapping_id] = dna_dist
        self.homo_dists[mapping_id] = homo_dist
        self.profiles[mapping_id] = profile
        self.centralities[mapping_id] = centralities
        self.phis[mapping_id] = phi
        self.psis[mapping_id] = psi
        self.intra_ssbonds[mapping_id] = intra_ssbond
        self.ssbond_lengths[mapping_id] = ssbond_length
        self.intra_links[mapping_id] = intra_link
        self.link_lengths[mapping_id] = link_length
        self.cis_conformations[mapping_id] = cis_conformation
        self.cis_followers[mapping_id] = cis_follower
        self.inter_chain_median_kds[mapping_id] = inter_chain_median_kd
        self.inter_chain_dist_weighted_kds[mapping_id] = inter_chain_dist_weighted_kd
        self.inter_chain_median_rsas[mapping_id] = inter_chain_median_rsa
        self.inter_chain_dist_weighted_rsas[mapping_id] = inter_chain_dist_weighted_rsa
        self.intra_chain_median_kds[mapping_id] = intra_chain_median_kd
        self.intra_chain_dist_weighted_kds[mapping_id] = intra_chain_dist_weighted_kd
        self.intra_chain_median_rsas[mapping_id] = intra_chain_median_rsa
        self.intra_chain_dist_weighted_rsas[mapping_id] = intra_chain_dist_weighted_rsa
        self.b_factors[mapping_id] = b_factor
        self.modres[mapping_id] = modres
        self.res_classes[mapping_id] = res_class,res_simple_class
        self.aa_ids[mapping_id] = identical_aa
        self.res_aas[mapping_id] = res_aa
        self.resolutions[mapping_id[0]] = resolution
        return

    def set_weighted_lig_dist(self,wd,conf):
        self.weighted_lig_dist = wd
        self.lig_dist_conf = conf
        return

    def set_weighted_metal_dist(self,wd,conf):
        self.weighted_metal_dist = wd
        self.metal_dist_conf = conf
        return

    def set_weighted_ion_dist(self,wd,conf):
        self.weighted_ion_dist = wd
        self.ion_dist_conf = conf
        return

    def set_weighted_chain_dist(self,wd,conf):
        self.weighted_chain_dist = wd
        self.chain_dist_conf = conf
        return

    def set_weighted_rna_dist(self,wd,conf):
        self.weighted_rna_dist = wd
        self.rna_dist_conf = conf
        return

    def set_weighted_dna_dist(self,wd,conf):
        self.weighted_dna_dist = wd
        self.dna_dist_conf = conf
        return

    def set_weighted_homo_dist(self,wd,conf):
        self.weighted_homo_dist = wd
        self.homo_dist_conf = conf
        return

    def set_weighted_phi(self,value):
        self.weighted_phi = value
        return

    def set_weighted_psi(self,value):
        self.weighted_psi = value
        return

    def set_weighted_intra_ssbond(self,value):
        self.weighted_intra_ssbond = value
        return

    def set_weighted_inter_ssbond(self,value):
        self.weighted_inter_ssbond = value
        return

    def set_weighted_ssbond_length(self,value):
        self.weighted_ssbond_length = value
        return

    def set_weighted_intra_link(self,value):
        self.weighted_intra_link = value
        return

    def set_weighted_inter_link(self,value):
        self.weighted_inter_link = value
        return

    def set_weighted_link_length(self,value):
        self.weighted_link_length = value
        return

    def set_weighted_cis_conformation(self,value):
        self.weighted_cis_conformation = value
        return

    def set_weighted_cis_follower(self,value):
        self.weighted_cis_follower = value
        return

    def set_weighted_inter_chain_median_kd(self,value):
        self.weighted_inter_chain_median_kd = value
        return

    def set_weighted_inter_chain_dist_weighted_kd(self,value):
        self.weighted_inter_chain_dist_weighted_kd = value
        return

    def set_weighted_inter_chain_median_rsa(self,value):
        self.weighted_inter_chain_median_rsa = value
        return

    def set_weighted_inter_chain_dist_weighted_rsa(self,value):
        self.weighted_inter_chain_dist_weighted_rsa = value
        return

    def set_weighted_intra_chain_median_kd(self,value):
        self.weighted_intra_chain_median_kd = value
        return

    def set_weighted_intra_chain_dist_weighted_kd(self,value):
        self.weighted_intra_chain_dist_weighted_kd = value
        return

    def set_weighted_intra_chain_median_rsa(self,value):
        self.weighted_intra_chain_median_rsa = value
        return

    def set_weighted_intra_chain_dist_weighted_rsa(self,value):
        self.weighted_intra_chain_dist_weighted_rsa = value
        return

    def set_weighted_modres(self,value):
        self.weighted_modres = value
        return

    def set_weighted_b_factor(self,value):
        self.weighted_b_factor = value
        return

    def get_raw_result(self):
        return (self.recommended_res, self.max_seq_res, self.weighted_surface_value, self.weighted_location, self.get_weighted_profile_str(),
                self.get_weighted_centralities_str(), self.weighted_b_factor, self.weighted_modres, self.weighted_ssa, self.weighted_phi,
                self.weighted_psi, self.weighted_intra_ssbond, self.weighted_inter_ssbond, self.weighted_ssbond_length, self.weighted_intra_link,
                self.weighted_inter_link, self.weighted_link_length, self.weighted_cis_conformation, self.weighted_cis_follower,
                self.weighted_inter_chain_median_kd, self.weighted_inter_chain_dist_weighted_kd, self.weighted_inter_chain_median_rsa,
                self.weighted_inter_chain_dist_weighted_rsa, self.weighted_intra_chain_median_kd, self.weighted_intra_chain_dist_weighted_kd,
                self.weighted_intra_chain_median_rsa, self.weighted_intra_chain_dist_weighted_rsa, self.weighted_lig_dist,
                self.weighted_metal_dist, self.weighted_ion_dist, self.weighted_chain_dist, self.weighted_rna_dist, self.weighted_dna_dist,
                self.weighted_homo_dist, self.rin_class, self.rin_simple_class, self.Class, self.simple_class, self.interaction_recommendations,
                self.lig_dist_conf, self.metal_dist_conf, self.ion_dist_conf, self.chain_dist_conf, self.rna_dist_conf, self.dna_dist_conf,
                self.homo_dist_conf, self.location_conf, self.classification_conf)

    def get_result_copy(self):
        result_obj = Mappings()
        result_obj.recommended_res = self.recommended_res
        result_obj.max_seq_res = self.max_seq_res
        result_obj.weighted_surface_value = self.weighted_surface_value
        result_obj.weighted_location = self.weighted_location
        result_obj.weighted_profile = self.weighted_profile
        result_obj.weighted_centralities = self.weighted_centralities
        result_obj.weighted_b_factor = self.weighted_b_factor
        result_obj.weighted_modres = self.weighted_modres
        result_obj.weighted_ssa = self.weighted_ssa
        result_obj.weighted_phi = self.weighted_phi
        result_obj.weighted_psi = self.weighted_psi
        result_obj.weighted_intra_ssbond = self.weighted_intra_ssbond
        result_obj.weighted_inter_ssbond = self.weighted_inter_ssbond
        result_obj.weighted_ssbond_length = self.weighted_ssbond_length
        result_obj.weighted_intra_link = self.weighted_intra_link
        result_obj.weighted_inter_link = self.weighted_inter_link
        result_obj.weighted_link_length = self.weighted_link_length
        result_obj.weighted_cis_conformation = self.weighted_cis_conformation
        result_obj.weighted_cis_follower = self.weighted_cis_follower
        result_obj.weighted_inter_chain_median_kd = self.weighted_inter_chain_median_kd
        result_obj.weighted_inter_chain_dist_weighted_kd = self.weighted_inter_chain_dist_weighted_kd
        result_obj.weighted_inter_chain_median_rsa = self.weighted_inter_chain_median_rsa
        result_obj.weighted_inter_chain_dist_weighted_rsa = self.weighted_inter_chain_dist_weighted_rsa
        result_obj.weighted_intra_chain_median_kd = self.weighted_intra_chain_median_kd
        result_obj.weighted_intra_chain_dist_weighted_kd = self.weighted_intra_chain_dist_weighted_kd
        result_obj.weighted_intra_chain_median_rsa = self.weighted_intra_chain_median_rsa
        result_obj.weighted_intra_chain_dist_weighted_rsa = self.weighted_intra_chain_dist_weighted_rsa
        result_obj.weighted_lig_dist = self.weighted_lig_dist
        result_obj.weighted_metal_dist = self.weighted_metal_dist
        result_obj.weighted_ion_dist = self.weighted_ion_dist
        result_obj.weighted_chain_dist = self.weighted_chain_dist
        result_obj.weighted_rna_dist = self.weighted_rna_dist
        result_obj.weighted_dna_dist = self.weighted_dna_dist
        result_obj.weighted_homo_dist = self.weighted_homo_dist
        result_obj.rin_class = self.rin_class
        result_obj.rin_simple_class = self.rin_simple_class
        result_obj.Class = self.Class
        result_obj.simple_class = self.simple_class
        result_obj.interaction_recommendations = self.interaction_recommendations
        result_obj.lig_dist_conf = self.lig_dist_conf
        result_obj.metal_dist_conf = self.metal_dist_conf
        result_obj.ion_dist_conf = self.ion_dist_conf
        result_obj.chain_dist_conf = self.chain_dist_conf
        result_obj.rna_dist_conf = self.rna_dist_conf
        result_obj.dna_dist_conf = self.dna_dist_conf
        result_obj.homo_dist_conf = self.homo_dist_conf
        result_obj.location_conf = self.location_conf
        result_obj.classification_conf = self.classification_conf
        return result_obj

    def weight_all(self,config,disorder_score,disorder_region):
        dist_maps = [(self.lig_dists,self.set_weighted_lig_dist), (self.metal_dists,self.set_weighted_metal_dist),
                        (self.ion_dists,self.set_weighted_ion_dist), (self.chain_dists,self.set_weighted_chain_dist),
                        (self.rna_dists,self.set_weighted_rna_dist), (self.dna_dists,self.set_weighted_dna_dist),
                        (self.homo_dists,self.set_weighted_homo_dist)]
        for dist_map,func in dist_maps:
            weighted_value,conf = self.weight(dist_map,config,distance_weighting = True)
            func(weighted_value,conf)

        value_maps = [(self.phis,self.set_weighted_phi), (self.psis,self.set_weighted_psi), (self.ssbond_lengths,self.set_weighted_ssbond_length),
                        (self.link_lengths,self.set_weighted_link_length), (self.inter_chain_median_kds,self.set_weighted_inter_chain_median_kd),
                        (self.inter_chain_dist_weighted_kds,self.set_weighted_inter_chain_dist_weighted_kd),
                        (self.inter_chain_median_rsas,self.set_weighted_inter_chain_median_rsa),
                        (self.inter_chain_dist_weighted_rsas,self.set_weighted_inter_chain_dist_weighted_rsa),
                        (self.intra_chain_median_kds,self.set_weighted_intra_chain_median_kd),
                        (self.intra_chain_dist_weighted_kds,self.set_weighted_intra_chain_dist_weighted_kd),
                        (self.intra_chain_median_rsas,self.set_weighted_intra_chain_median_rsa),
                        (self.intra_chain_dist_weighted_rsas,self.set_weighted_intra_chain_dist_weighted_rsa),
                        (self.b_factors,self.set_weighted_b_factor)]
        for value_map,func in value_maps:
            weighted_value = self.weight(value_map,config,calc_conf = False)
            func(weighted_value)

        prop_maps = [(self.cis_conformations,self.set_weighted_cis_conformation), (self.cis_followers,self.set_weighted_cis_follower),
                        (self.modres,self.set_weighted_modres)]
        for prop_map,func in prop_maps:
            weighted_prop = self.weight_prop(prop_map)
            func(weighted_prop)

        bool_prop_maps = [(self.intra_ssbonds,self.set_weighted_intra_ssbond,self.set_weighted_inter_ssbond),
                            (self.intra_links,self.set_weighted_intra_link,self.set_weighted_inter_link)]
        for prop_map,true_func,false_func in bool_prop_maps:
            weighted_true_prop,weighted_false_prop = self.weigthed_bool_prop(prop_map)
            true_func(weighted_true_prop)
            false_func(weighted_false_prop)

        self.weighted_ssa = self.weight_majority(self.ssas)

        self.weight_structure_location(config)
        self.weight_centralities()
        self.weight_profiles()
        self.rin_based_classify()
        self.distance_classify(config,disorder_score,disorder_region)
        self.set_recommended_residues()

        return

    def weight(self,value_map,config,distance_weighting = False,calc_conf = True):
        nom = 0.0
        denom = 0.0
        n = 0.0
        qs = []
        weighted_value = None
        conf = 0.0
        for mapping_id in value_map:
            value = value_map[mapping_id]
            if value == None:
                continue
            if distance_weighting and value > config.long_distance_threshold:
                continue
            qual = self.qualities[mapping_id]
            nom += value*qual
            denom += qual
            n += 1.0
            qs.append(qual)
        if denom > 0.0:
            weighted_value = nom/denom
            if calc_conf:
                conf = (1.0-1.0/(n+1.0))*(max(qs)+median(qs))/2
        if calc_conf:
            return weighted_value,conf
        else:
            return weighted_value

    def weight_prop(self,value_map):
        prop = 0.
        if len(value_map) == 0:
            return None
        for mapping_id in value_map:
            value = value_map[mapping_id]
            if value == None:
                continue
            prop += 1.
        weighted_prop = prop/float(len(value_map))
        return weighted_prop

    def weigthed_bool_prop(self,value_map):
        true_prop = 0.
        false_prop = 0.
        if len(value_map) == 0:
            return None,None
        for mapping_id in value_map:
            value = value_map[mapping_id]
            if value == None:
                continue
            if value:
                true_prop += 1.
            else:
                false_prop += 1.
        weighted_true_prop = true_prop/float(len(value_map))
        weighted_false_prop = false_prop/float(len(value_map))
        return weighted_true_prop,weighted_false_prop

    def weight_majority(self,value_map):
        voting = {}
        best_value = None
        for mapping_id in value_map:
            qual = self.qualities[mapping_id]
            value = value_map[mapping_id]
            if not value in voting:
                voting[value] = qual
            else:
                voting[value] += qual

        max_qual = 0.
        for value in voting:
            qual = voting[value]
            if qual > max_qual:
                max_qual = qual
                best_value = value

        return best_value

    def weight_structure_location(self,config):

        DSC = 0.0 #decision sum core
        DSS = 0.0 #decision sum surface
        structure_location_qualities = []
        max_cov = 0.
        max_cov_rsa = None

        for mapping_id in self.rsas:
            rsa = self.rsas[mapping_id]
            if rsa == None:
                continue
            qual = self.qualities[mapping_id]
            cov = self.covs[mapping_id]

            if cov > max_cov:
                max_cov = cov
                max_cov_rsa = rsa
            if rsa < config.surface_threshold:
                DSC += qual*(cov**5)
            else:
                DSS += qual*(cov**5)
            structure_location_qualities.append(qual)

        if len(structure_location_qualities) == 0:
            self.weighted_surface_value = None
            self.weighted_location = None
            self.max_cov = None
            self.location_conf = None
            return

        weighted_surface_value = DSS-2*DSC
        if weighted_surface_value > 0:
            weighted_location = "Surface"
        else:
            weighted_location = "Core"

        self.weighted_surface_value = weighted_surface_value
        self.weighted_location = weighted_location

        conf_sc = (1.0-1.0/(len(structure_location_qualities)+1.0))*(abs(weighted_surface_value)/(DSS+2*DSC))*median(structure_location_qualities)

        self.max_cov = max_cov
        self.location_conf = conf_sc
        return

    def weight_centralities(self):
        self.weighted_centralities = None

        total_qual = 0.0

        for mapping_id in self.centralities:
            centrality_scores = self.centralities[mapping_id]
            if centrality_scores == None:
                continue
            qual = self.qualities[mapping_id]
            if self.weighted_centralities == None:
                self.weighted_centralities = [0.]*len(centrality_scores.cent_list)

            for pos,cent_score in enumerate(centrality_scores.cent_list):
                if cent_score == None:
                    continue
                self.weighted_centralities[pos] += cent_score*qual

            total_qual += qual
        if self.weighted_centralities == None:
            return

        if total_qual > 0.0:
            for pos,cent_score in enumerate(self.weighted_centralities):
                self.weighted_centralities[pos] = self.weighted_centralities[pos]/total_qual
            self.weighted_centralities = rin.Centrality_scores(cent_list = self.weighted_centralities)

        return

    def get_weighted_centralities(self):
        if self.weighted_centralities != None:
            return self.weighted_centralities
        if self.weighted_centralities_str == None:
            return None
        return rin.Centrality_scores(code_str = self.weighted_centralities_str)

    def get_weighted_centralities_str(self):
        if self.weighted_centralities_str != None:
            return self.weighted_centralities_str
        if self.weighted_centralities == None:
            return None
        return self.weighted_centralities.str_encode()

    def get_weighted_profile(self):
        if self.weighted_profile != None:
            return self.weighted_profile
        if self.weighted_profile_str == None:
            return None
        return rin.Interaction_profile(profile_str= self.weighted_profile_str)

    def get_weighted_profile_str(self):
        if self.weighted_profile_str != None:
            return self.weighted_profile_str
        if self.weighted_profile == None:
            return None
        return self.weighted_profile.encode()

    def weight_profiles(self):
        weight_profile_tuples = []
        for mapping_id in self.profiles:
            weight_profile_tuples.append((self.qualities[mapping_id],self.profiles[mapping_id]))
        self.weighted_profile = rin.calculateAverageProfile(weight_profile_tuples)
        return

    def set_recommended_residues(self):
        max_seq_id = 0.
        max_identical_aa_identical_class = 0.
        max_identical_aa = 0.
        max_identical_class = 0.
        max_should_not_happen = 0.
        max_qual = 0.
        max_qual_res = None
        self.max_seq_res = None
        self.recommended_res = None
        self.interaction_recommendations = {}
        for mapping_id in self.seq_ids:
            seq_id = self.seq_ids[mapping_id]
            cov = self.covs[mapping_id]
            qual = self.qualities[mapping_id]
            if seq_id > max_seq_id:
                max_seq_id = seq_id
                self.max_seq_res = qual,(mapping_id,seq_id,cov)
            elif seq_id == max_seq_id:
                if mapping_id[0] > self.max_seq_res[1][0][0]:
                    max_seq_id = seq_id
                    self.max_seq_res = qual,(mapping_id,seq_id,cov)

            if qual > max_qual:
                max_qual = qual
                max_qual_res = mapping_id,seq_id,cov
            elif qual == max_qual:
                if mapping_id[0] > max_qual_res[0][0]:
                    max_qual = qual
                    max_qual_res = mapping_id,seq_id,cov

            residue_simple_class = self.res_classes[mapping_id][1]

            if self.rin_simple_class == residue_simple_class:
                if self.aa_ids[mapping_id]:
                    if qual > max_identical_aa_identical_class:
                        max_identical_aa_identical_class = qual
                        max_identical_aa_identical_class_res = mapping_id,seq_id,cov
                    elif qual == max_identical_aa_identical_class:
                        if mapping_id[0] > max_identical_aa_identical_class_res[0][0]:
                            max_identical_aa_identical_class = qual
                            max_identical_aa_identical_class_res = mapping_id,seq_id,cov
                else:
                    if qual > max_identical_class:
                        max_identical_class = qual
                        max_identical_class_res = mapping_id,seq_id,cov
                    elif qual == max_identical_class:
                        if mapping_id[0] > max_identical_class_res[0][0]:
                            max_identical_class = qual
                            max_identical_class_res = mapping_id,seq_id,cov
            else:
                if self.aa_ids[mapping_id]:
                    if qual > max_identical_aa:
                        max_identical_aa = qual
                        max_identical_aa_res = mapping_id,seq_id,cov
                    elif qual == max_identical_aa:
                        if mapping_id[0] > max_identical_aa_res[0][0]:
                            max_identical_aa = qual
                            max_identical_aa_res = mapping_id,seq_id,cov
                else:
                    if qual > max_should_not_happen:
                        max_should_not_happen = qual
                        max_should_not_happen_res = mapping_id,seq_id,cov
                    elif qual == max_should_not_happen:
                        if mapping_id[0] > max_should_not_happen_res[0][0]:
                            max_should_not_happen = qual
                            max_should_not_happen_res = mapping_id,seq_id,cov

            if residue_simple_class != None:
                if residue_simple_class.count('Interaction') > 0:
                    if residue_simple_class not in self.interaction_recommendations:
                        self.interaction_recommendations[residue_simple_class] = (qual,mapping_id,seq_id,cov)
                    elif qual > self.interaction_recommendations[residue_simple_class][0]:
                        self.interaction_recommendations[residue_simple_class] = (qual,mapping_id,seq_id,cov)

        if max_identical_aa_identical_class > 0.:
            self.recommended_res = max_identical_aa_identical_class,max_identical_aa_identical_class_res
        elif max_identical_class > 0.:
            self.recommended_res = max_identical_class,max_identical_class_res
        elif max_identical_aa > 0.:
            self.recommended_res = max_identical_aa,max_identical_aa_res
        else:
            #If this case is reached, then max_qual_res should be identical to max_should_not_happen_res
            #This would mean there is not mapped residue with identical simple class or identical aa, which should not happen (thus the name)
            self.recommended_res = max_qual,max_qual_res

        if self.recommended_res != None:
            if self.recommended_res[1] != None:
                qual,(recommended_res,seq_id,cov) = self.recommended_res
                pdb_id,chain,res_nr = recommended_res
                resolution = self.resolutions[pdb_id]
                res_aa = self.res_aas[(pdb_id,chain,res_nr)]
                self.recommended_res = '%s:%s %s:%s;%1.2f;%1.2f;%1.1f' % (pdb_id,chain,res_nr,res_aa,seq_id,cov,resolution)
            else:
                self.recommended_res = None

        if self.max_seq_res != None:
            if self.max_seq_res[1] != None:
                qual,(recommended_res,seq_id,cov) = self.max_seq_res
                pdb_id,chain,res_nr = recommended_res
                resolution = self.resolutions[pdb_id]
                res_aa = self.res_aas[(pdb_id,chain,res_nr)]
                self.max_seq_res = '%s:%s %s:%s;%1.2f;%1.2f;%1.1f' % (pdb_id,chain,res_nr,res_aa,seq_id,cov,resolution)
            else:
                self.max_seq_res = None

        return

    def get_recommended_res_str(self):
        return self.recommended_res

    def get_max_seq_structure_res_str(self):
        return self.max_seq_res

    def distance_classify(self,config,disorder_score,disorder_region):
        self.Class,self.classification_conf = self.getWeightedClass(config)
        if self.Class == None:
            self.Class,self.classification_conf = self.disorderCheck(config,disorder_score,disorder_region)
            self.simple_class = self.Class
            self.rin_class = self.Class
            self.rin_simple_class = self.Class
        else:
            self.simple_class = self.simplifyClass(self.Class,self.weighted_location)
        return

    def rin_based_classify(self):
        raw_rin_class,raw_rin_simple_class = self.weighted_profile.getClass()
        if raw_rin_class == 'No interaction':
            self.rin_class = self.weighted_location
        else:
            self.rin_class = raw_rin_class
        if raw_rin_simple_class == 'No interaction':
            self.rin_simple_class = self.weighted_location
        else:
            self.rin_simple_class = raw_rin_simple_class
        return

    def simplifyClass(self,c,sc):
        if c == "Surface" or c == "Core" or c == 'Disorder' or c == None:
            return c
        interactions = re.sub(r' Interaction$','',re.sub(r'^[^:]*: ','',c)).split(' and ')
        non_far_interactions = set([ x for x in interactions if ' far' not in x ])

        priorities = ['Metal', 'Ligand', 'DNA', 'RNA', 'Protein', 'Ion']

        if len(non_far_interactions) == 0:
            return sc

        for interaction in priorities:
            if interaction in non_far_interactions:
                return interaction + ' Interaction'

        print('Unknown Class: %s' % c)
        return sc

    def disorderCheck(self,config,disorder_score,region):
        if region == None or disorder_score == None:
            return None,None
        glob = True
        if region == 'disorder':
            glob = False
        elif region == 'globular':
            glob = True
        elif region in MobiDB_map:
            glob = False
        elif config.verbose:
            print('Unknown disorder region: ',region)
        if not glob:
            return 'Disorder',disorder_score
        else:
            return None,None

    def getWeightedClass(self,config):
        dt = config.short_distance_threshold

        contacts = []
        confs = []

        # Determine the near macromolecule with the highest conf.
        # If there is no near macromolecule, also consider far ones.
        # In case of multiple maximum confs use the first one added to macros.
        macros = []
        macro_confs = []

        if self.weighted_chain_dist and self.weighted_chain_dist < dt:
            macros.append('Protein')
            macro_confs.append(self.chain_dist_conf)

        if self.weighted_dna_dist and self.weighted_dna_dist < dt:
            macros.append('DNA')
            macro_confs.append(self.dna_dist_conf)

        if self.weighted_rna_dist and self.weighted_rna_dist < dt:
            macros.append('RNA')
            macro_confs.append(self.rna_dist_conf)

        contact_macro = None
        conf_macro = None

        if len(macros) > 0:
            best_macro_idx = np.min(np.argmax(macro_confs))
            contacts.append(macros[best_macro_idx])
            confs.append(macro_confs[best_macro_idx])

        # Additional interactions
        if self.weighted_lig_dist and self.weighted_lig_dist < dt:
            contacts.append('Ligand')
            confs.append(self.lig_dist_conf)

        if self.weighted_metal_dist and self.weighted_metal_dist < dt:
            contacts.append('Metal')
            confs.append(self.metal_dist_conf)

        if self.weighted_ion_dist and self.weighted_ion_dist < dt:
            contacts.append('Ion')
            confs.append(self.ion_dist_conf)

        if len(contacts) == 0:
            return self.weighted_location,self.location_conf

        clas = " and ".join(contacts)
        if len(contacts) == 4:
            clas = "Quadruple Interaction: " + clas
        elif len(contacts) == 3:
            clas = "Triple Interaction: " + clas
        elif len(contacts) == 2:
            clas = "Double Interaction: " + clas
        elif len(contacts) == 1:
            clas = clas + " Interaction"

        conf = sum(confs) / len(confs)

        return clas,conf

def get_shortest_distances(chains,lig_dists,chain_distances,homomer_distances):
    min_ld = None
    min_md = None
    min_id = None
    min_cd = None
    min_dd = None
    min_rd = None
    min_hd = None

    ldists = {}
    mdists = {}
    idists = {}
    if lig_dists != None:
        for lig_id in lig_dists:
            lig_name,res,chain = lig_id.split('_')
            (dist,atom_pair) = lig_dists[lig_id]
            if lig_name in metal_atoms:
                mdists[dist] = lig_name,res,chain
            elif lig_name in ion_atoms:
                idists[dist] = lig_name,res,chain
            else:
                ldists[dist] = lig_name,res,chain

    min_lig = None
    min_metal = None
    min_ion = None

    if len(ldists) > 0:
        min_ld = min(ldists.keys())
        min_lig = ldists[min_ld]

    if len(mdists) > 0:
        min_md = min(mdists.keys())
        min_metal = mdists[min_md]

    if len(idists) > 0:
        min_id = min(idists.keys())
        min_ion = idists[min_id]

    cdists = {}
    ddists = {}
    rdists = {}

    iacs = {}

    if chain_distances != None:
        for chain_id in chain_distances:
            (dist,atom_pair,min_resi) = chain_distances[chain_id]
            if dist == None:
                continue
            if not chain_id in chains:
                chaintype = 'Protein'
            else:
                chaintype = chains[chain_id]

            if chaintype == "Protein" or chaintype == 'Peptide':
                cdists[dist] = chain_id
            elif chaintype == "RNA":
                rdists[dist] = chain_id
            elif chaintype == "DNA":
                ddists[dist] = chain_id

    if len(cdists) > 0:
        min_cd = min(cdists.keys())
        iacs['Protein'] = cdists[min_cd]
    if len(rdists) > 0:
        min_rd = min(rdists.keys())
        iacs['RNA'] = rdists[min_rd]
    if len(ddists) > 0:
        min_dd = min(ddists.keys())
        iacs['DNA'] = ddists[min_dd]

    homo_dists = []
    if homomer_distances != None:
        for homo_chain in homomer_distances:
            dist = homomer_distances[homo_chain]
            homo_dists.append(dist)
    if len(homo_dists) > 0:
        min_hd = min(homo_dists)

    minimal_distances = []
    if min_cd != None:
        minimal_distances.append(min_cd)
    if min_dd != None:
        minimal_distances.append(min_dd)
    if min_rd != None:
        minimal_distances.append(min_rd)
    if min_ld != None:
        minimal_distances.append(min_ld)
    if min_md != None:
        minimal_distances.append(min_md)
    if min_id != None:
        minimal_distances.append(min_id)

    if len(minimal_distances) == 0:
        min_minimal_distances = 2.0
    else:
        min_minimal_distances = min(minimal_distances)

    if min_minimal_distances < 1.2:
        return None

    return min_hd,min_ld,min_md,min_id,min_cd,min_rd,min_dd,min_lig,min_metal,min_ion,iacs


class Residue:
    __slots__ = ['res_num','aa','lig_dist_str','lig_dists','chain_dist_str','chain_distances','RSA','SSA',
                'homo_dist_str','homomer_distances',
                'interaction_profile','interaction_profile_str','centrality_score_str','centralities','modres','b_factor','database_id','stored','phi','psi',
                'intra_ssbond','ssbond_length','intra_link','link_length','cis_conformation','cis_follower',
                'inter_chain_median_kd','inter_chain_dist_weighted_kd','inter_chain_median_rsa',
                'inter_chain_dist_weighted_rsa','intra_chain_median_kd','intra_chain_dist_weighted_kd',
                'intra_chain_median_rsa','intra_chain_dist_weighted_rsa','Class','simpleClass']

    def __init__(self,res_num,aa = 'X',lig_dist_str = None,lig_dists = None,chain_dist_str = None,chain_distances = None,RSA = None,
                    SSA = None,homo_dist_str = None,homomer_distances = None,interaction_profile = None,interaction_profile_str = None,
                    centrality_score_str = None, centralities = None,modres = None,b_factor = None,database_id = None,stored = False,phi = None,psi = None,
                    intra_ssbond = None,ssbond_length = None,intra_link = None,link_length = None,cis_conformation = None,cis_follower = None,
                    inter_chain_median_kd = None ,inter_chain_dist_weighted_kd = None,inter_chain_median_rsa = None,
                    inter_chain_dist_weighted_rsa = None ,intra_chain_median_kd = None, intra_chain_dist_weighted_kd = None,
                    intra_chain_median_rsa = None, intra_chain_dist_weighted_rsa = None):

        self.res_num = res_num
        self.aa = aa
        self.lig_dist_str = lig_dist_str
        self.lig_dists = lig_dists
        self.chain_dist_str = chain_dist_str
        self.chain_distances = chain_distances
        self.RSA = RSA #relative surface accessible area
        self.SSA = SSA #secondary structure assignment
        self.homomer_distances = homomer_distances
        self.homo_dist_str = homo_dist_str
        self.interaction_profile = interaction_profile
        self.interaction_profile_str = interaction_profile_str #interaction profile coded into a string, see rin.py for more information
        self.centralities = centralities
        self.centrality_score_str = centrality_score_str
        self.modres = modres
        self.b_factor = b_factor
        self.database_id = database_id
        self.stored = stored
        self.phi = phi
        self.psi = psi
        self.intra_ssbond = intra_ssbond #True if ssbond to another chain, False if ssbond to same chain, None if no ssbond
        self.ssbond_length = ssbond_length
        self.intra_link = intra_link
        self.link_length = link_length
        self.cis_conformation = cis_conformation
        self.cis_follower = cis_follower
        self.inter_chain_median_kd = inter_chain_median_kd
        self.inter_chain_dist_weighted_kd = inter_chain_dist_weighted_kd
        self.inter_chain_median_rsa = inter_chain_median_rsa
        self.inter_chain_dist_weighted_rsa = inter_chain_dist_weighted_rsa
        self.intra_chain_median_kd = intra_chain_median_kd
        self.intra_chain_dist_weighted_kd = intra_chain_dist_weighted_kd
        self.intra_chain_median_rsa = intra_chain_median_rsa
        self.intra_chain_dist_weighted_rsa = intra_chain_dist_weighted_rsa
        self.Class = None
        self.simpleClass = None

    def set_database_id(self,value):
        self.database_id = value

    def get_database_id(self):
        return self.database_id

    def get_aa(self):
        return self.aa

    def get_chain_distances(self):
        return self.chain_distances

    def get_chain_dist_str(self):
        if self.chain_dist_str == None:
            if self.chain_distances == None:
                return None
            self.convert_chain_dist_str()
        return self.chain_dist_str

    def convert_chain_dist_str(self):
        chain_dist_strings = []
        for chain_id in self.chain_distances:
            (dist,atom_pair,min_resi) = self.chain_distances[chain_id]
            if not dist == None:
                chain_dist_strings.append("%s.%s:%1.2f:(%s-%s)" % (chain_id,min_resi,dist,str(atom_pair[0]),str(atom_pair[1])))

        self.chain_dist_str = ",".join(chain_dist_strings)
        return

    def get_ligand_distances(self):
        return self.lig_dists

    def get_lig_dist_str(self):
        if self.lig_dist_str == None:
            if self.lig_dists == None:
                return None
            self.convert_lig_dist_str()
        return self.lig_dist_str

    def convert_lig_dist_str(self):
        lig_dist_strings = []

        for lig_id in self.lig_dists:
            (dist,atom_pair) = self.lig_dists[lig_id]
            if dist != None:
                lig_dist_strings.append("%s:%1.2f:(%s-%s)" % (lig_id,dist,str(atom_pair[0]),str(atom_pair[1])))

        self.lig_dist_str = ",".join(lig_dist_strings)
        return

    def get_shortest_distances(self,chains):
        min_ld = None
        min_md = None
        min_id = None
        min_cd = None
        min_dd = None
        min_rd = None
        min_hd = None

        ldists = {}
        mdists = {}
        idists = {}
        if self.lig_dists != None:
            for lig_id in self.lig_dists:
                lig_name,res,chain = lig_id.split('_')
                (dist,atom_pair) = self.lig_dists[lig_id]
                if lig_name in metal_atoms:
                    mdists[dist] = lig_name,res,chain
                elif lig_name in ion_atoms:
                    idists[dist] = lig_name,res,chain
                else:
                    ldists[dist] = lig_name,res,chain

        min_lig = None
        min_metal = None
        min_ion = None

        if len(ldists) > 0:
            min_ld = min(ldists.keys())
            min_lig = ldists[min_ld]

        if len(mdists) > 0:
            min_md = min(mdists.keys())
            min_metal = mdists[min_md]

        if len(idists) > 0:
            min_id = min(idists.keys())
            min_ion = idists[min_id]

        cdists = {}
        ddists = {}
        rdists = {}

        iacs = {}

        if self.chain_distances != None:
            for chain_id in self.chain_distances:
                (dist,atom_pair,min_resi) = self.chain_distances[chain_id]
                if dist == None:
                    continue
                if not chain_id in chains:
                    chaintype = 'Protein'
                else:
                    chaintype = chains[chain_id]

                if chaintype == "Protein" or chaintype == 'Peptide':
                    cdists[dist] = chain_id
                elif chaintype == "RNA":
                    rdists[dist] = chain_id
                elif chaintype == "DNA":
                    ddists[dist] = chain_id

        if len(cdists) > 0:
            min_cd = min(cdists.keys())
            iacs['Protein'] = cdists[min_cd]
        if len(rdists) > 0:
            min_rd = min(rdists.keys())
            iacs['RNA'] = rdists[min_rd]
        if len(ddists) > 0:
            min_dd = min(ddists.keys())
            iacs['DNA'] = ddists[min_dd]

        homo_dists = []
        if self.homomer_distances != None:
            for homo_chain in self.homomer_distances:
                dist = self.homomer_distances[homo_chain]
                homo_dists.append(dist)
        if len(homo_dists) > 0:
            min_hd = min(homo_dists)

        return min_hd,min_ld,min_md,min_id,min_cd,min_rd,min_dd,min_lig,min_metal,min_ion,iacs

    def get_homomer_dists(self):
        return self.homomer_distances

    def get_homo_dist_str(self):
        if self.homo_dist_str == None:
            if self.homomer_distances == None:
                return None
            self.convert_homo_dist_str()
        return self.homo_dist_str

    def convert_homo_dist_str(self):
        homo_strs = []
        for homo_chain in self.homomer_distances:
            min_d = self.homomer_distances[homo_chain]
            homo_strs.append('%s:%1.2f' % (homo_chain,min_d))
        self.homo_dist_str = ','.join(homo_strs)
        return

    def get_centralities(self,get_whats_there = False):
        if get_whats_there:
            if self.centralities != None:
                return self.centralities
            return self.centrality_score_str

        if self.centrality_score_str != None and self.centralities == None:
            self.centralities = rin.Centrality_scores(code_str = self.centrality_score_str)
        return self.centralities

    def get_centrality_str(self):
        if self.centrality_score_str == None:
            if self.centralities == None:
                return None
            else:
                self.convert_centrality_str()
        return self.centrality_score_str

    def convert_centrality_str(self):
        self.centrality_score_str = self.centralities.str_encode()
        return

    def get_modres(self):
        return self.modres

    def get_b_factor(self):
        return self.b_factor

    def get_rsa(self):
        return self.RSA

    def get_ssa(self):
        return self.SSA

    def get_phi(self):
        return self.phi

    def get_psi(self):
        return self.psi

    def get_angles(self):
        return self.phi,self.psi

    def get_residue_link_information(self):
        return (self.intra_ssbond, self.ssbond_length, self.intra_link, self.link_length, self.cis_conformation, self.cis_follower)

    def get_interaction_profile(self,get_whats_there = False):
        if get_whats_there:
            if self.interaction_profile != None:
                return self.interaction_profile
            return self.interaction_profile_str

        if self.interaction_profile_str != None and self.interaction_profile == None:
            self.interaction_profile = rin.Interaction_profile(profile_str = self.interaction_profile_str)
        return self.interaction_profile

    def get_interaction_profile_str(self):
        if self.interaction_profile_str == None:
            if self.interaction_profile == None:
                return None
            else:
                self.interaction_profile_str = self.interaction_profile.encode()
        return self.interaction_profile_str

    def get_milieu(self):
        return (self.inter_chain_median_kd,self.inter_chain_dist_weighted_kd,
                self.inter_chain_median_rsa,self.inter_chain_dist_weighted_rsa,self.intra_chain_median_kd,
                self.intra_chain_dist_weighted_kd,self.intra_chain_median_rsa,self.intra_chain_dist_weighted_rsa)

    def set_classification(self,Class,simpleClass):
        self.Class = Class
        self.simpleClass = simpleClass
        return
