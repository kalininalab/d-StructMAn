from structman.lib.sdsc.mappings import Mappings
from structman.lib.sdsc.mutations import SNV
from structman.lib.sdsc.utils import process_recommend_structure_str


class Position:
    __slots__ = ['pos', 'wt_aa', 'mut_aas', 'pos_tags', 'stored', 'database_id', 'pdb_res_nr', 'checked', 'mappings',
                 'disorder_score', 'disorder_region', 'recommended_structure']

    def __init__(self, pos=None, wt_aa='X', mut_aas=set(), tags=set(), mut_tags_map={}, pdb_res_nr=None,
                 checked=False, database_id=None, recommended_structure=None):
        self.pos = pos  # the position in a protein sequence, int or None
        self.pdb_res_nr = pdb_res_nr  # this is a string, due to insertion-codes
        self.wt_aa = wt_aa  # wildtype amino acid type in one letter code
        self.mut_aas = {}
        for new_aa in mut_aas:
            if new_aa in mut_tags_map:
                mut_tags = mut_tags_map[new_aa]
            else:
                mut_tags = set()
            snv = SNV(new_aa, tags=mut_tags)
            self.mut_aas[new_aa] = snv
        self.pos_tags = tags.copy()
        self.database_id = database_id
        self.stored = (database_id is not None)
        self.checked = checked  # Flag for sanity checks
        self.mappings = Mappings()
        self.disorder_score = None
        self.disorder_region = None
        self.recommended_structure = recommended_structure

    def print_state(self):
        print(self.wt_aa, self.pos, self.mut_aas, self.pos_tags)

    def check(self, wt_aa, overwrite=False):
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

    def fuse(self, position):
        warn = False
        if self.pos != position.pos:
            raise NameError('Cannot fuse positions with differing pos numbers')
            return
        if self.pdb_res_nr != position.pdb_res_nr:
            raise NameError('Cannot fuse positions with differing pdb_res_nr')
            return
        if self.wt_aa != position.wt_aa:
            #print('Warning: fuse positions with different WT AAs:',self.pos,self.pdb_res_nr,self.wt_aa,position.wt_aa)
            warn = True

        self.pos_tags = self.pos_tags | position.pos_tags
        for aa in position.mut_aas:
            if aa not in self.mut_aas:
                self.mut_aas[aa] = position.mut_aas[aa]
            else:
                self.mut_aas[aa].fuse(position.mut_aas[aa])

        return warn

    def add_tags(self, tags):
        self.pos_tags = self.pos_tags | tags

    def get_pos_tags(self):
        return self.pos_tags

    def getAACBase(self):
        return '%s%s' % (self.wt_aa, self.pos)

    def set_stored(self, value):
        self.stored = value

    def get_stored(self):
        return self.stored

    def set_database_id(self, value):
        self.database_id = value

    def get_database_id(self):
        return self.database_id

    def get_res_id(self):
        return self.pdb_res_nr

    def get_wt_aa(self):
        return self.wt_aa

    def get_mut_aas(self):
        return self.mut_aas

    def get_mut_tags(self, aa):
        return self.mut_aas[aa].tags

    def set_disorder_score(self, disorder_score):
        self.disorder_score = disorder_score

    def get_disorder_score(self):
        return self.disorder_score

    def set_disorder_region(self, region):
        self.disorder_region = region

    def get_disorder_region(self):
        return self.disorder_region

    def add_pos_res_mapping(self, pdb_id, chain, res_nr, mapping):
        self.mappings.add_mapping((pdb_id, chain, res_nr), mapping)

    def classify(self, config):
        self.mappings.weight_all(config, self.disorder_score, self.disorder_region)

    def get_classification(self):
        return self.mappings.rin_simple_class

    def get_recommended_structure(self):
        if self.recommended_structure is not None:
            return self.recommended_structure
        recommended_structure, seq_id, cov, resolution = process_recommend_structure_str(self.mappings.recommended_res)
        self.recommended_structure = recommended_structure
        return recommended_structure
