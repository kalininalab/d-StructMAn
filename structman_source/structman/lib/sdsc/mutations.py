from structman.lib import sdsc


class SNV:
    __slots__ = ['new_aa', 'database_id', 'stored', 'tags']

    def __init__(self, new_aa, tags=set(), database_id=None):
        self.new_aa = new_aa
        self.tags = tags.copy()
        self.database_id = database_id
        self.stored = (database_id is not None)

    def fuse(self, other_snv):
        self.tags = self.tags | other_snv.tags


class MultiMutation:
    __slots__ = ['wt_prot', 'mut_prot', 'snvs', 'indels', 'database_id', 'stored', 'tags']

    def __init__(self, wt_prot, mut_prot, mutation_list, tags=set()):
        from structman.lib.sdsc.indel import Indel

        self.wt_prot = wt_prot
        self.mut_prot = mut_prot
        self.snvs = {}
        self.indels = []
        self.tags = tags.copy()
        for mut in mutation_list:
            if not isinstance(mut, tuple):
                self.indels.append(mut)
            else:
                self.snvs[mut[0].pos] = (mut[0].mut_aas[mut[1]])
        self.database_id = None
        self.stored = False

    def get_snv_db_ids(self):
        db_ids = []
        for pos in self.snvs:
            snv = self.snvs[pos]
            db_ids.append(snv.database_id)
        return db_ids

    def get_indel_db_ids(self):
        db_ids = []
        for indel in self.indels:
            db_ids.append(indel.database_id)
        return db_ids

    def mutate(self, proteins, config):
        if proteins[self.mut_prot].sequence is None:
            if proteins[self.wt_prot].sequence is None:
                config.errorlog.add_error('Sequence is None: %s %s' % (self.wt_prot, self.mut_prot))
            mut_seq = proteins[self.wt_prot].sequence
            indel_position_order = {}
            for pos in self.snvs:
                wt_aa = proteins[self.wt_prot].positions[pos].wt_aa
                aa2 = self.snvs[pos].new_aa
                mut_seq = '%s%s%s' % (mut_seq[:pos], aa2, mut_seq[(pos + 1):])

            for indel in self.indels:
                indel_position_order[indel.left_flank_pos] = indel

            pos_sorted = sorted(indel_position_order.keys(), reverse=True)
            for pos in pos_sorted:
                mut_seq = indel_position_order[pos].return_mutate_sequence(mut_seq)
            proteins[self.mut_prot].sequence = mut_seq
        else:
            mut_seq = proteins[self.mut_prot].sequence

        for (pos, aa) in enumerate(mut_seq):
            seq_pos = pos + 1
            position = sdsc.Position(pos=seq_pos, wt_aa=aa, checked=True)
            proteins[self.mut_prot].positions[seq_pos] = position
