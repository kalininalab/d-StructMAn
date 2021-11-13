from structman.lib.sdsc.sdsc_utils import translate

from structman.lib.sdsc.mappings import Mappings
from structman.lib.sdsc import position as position_package

def majority_vote(secs):
    class_dict = {}
    for sec in secs:
        if sec not in class_dict:
            class_dict[sec] = 1
        else:
            class_dict[sec] += 1

    max_sec = None
    max_qual = 0
    for sec in class_dict:
        qual = class_dict[sec]
        if qual > max_qual:
            max_qual = qual
            max_sec = sec
    return max_sec

def aggregate_positions(position_objects, config):
    disorder_scores = []
    disorder_regions = []
    mappings = Mappings()
    if position_objects is None:
        return mappings, None, None
    for pos_obj in position_objects:
        result = pos_obj.mappings.get_result_for_indel_aggregation()
        mappings.add_result(pos_obj.pos, result, 1.0, 1.0, 1.0)
        disorder_scores.append(pos_obj.disorder_score)
        disorder_regions.append(pos_obj.disorder_region)

    if len(disorder_scores) > 0:
        disorder_score = sum(disorder_scores)/float(len(disorder_scores))
        disorder_region = majority_vote(disorder_regions)
    else:
        disorder_score = None
        disorder_region = None

    mappings.weight_all(config, disorder_score, disorder_region, for_indel_aggregation = True)
    return mappings, disorder_score, disorder_region

class Indel:
    __slots__ = ['left_flank', 'left_flank_pos', 'right_flank', 'right_flank_pos', 'tags', 'stored', 'database_id', 'size',
                 'wt_prot', 'mut_prot', 'terminal', 'terminal_end', 'delta_delta_classification', 'wt_aggregates', 'mut_aggregates',
                 'left_flank_wt_aggregates', 'left_flank_mut_aggregates', 'right_flank_wt_aggregates', 'right_flank_mut_aggregates'
                ]

    def __init__(self, left_flank=None, right_flank=None, tags=set()):
        self.terminal = False
        self.terminal_end = None
        self.left_flank = left_flank
        if left_flank is not None:
            self.left_flank_pos = int(left_flank[1:])
        else:
            self.left_flank_pos = 0
            self.terminal = True
            self.terminal_end = 'left'
        self.right_flank = right_flank
        if right_flank is not None:
            self.right_flank_pos = int(right_flank[1:])
        else:
            self.right_flank_pos = None
            self.terminal = True
            self.terminal_end = 'right'
        self.tags = tags.copy()
        self.stored = False
        self.database_id = None
        self.wt_prot = None
        self.mut_prot = None
        self.delta_delta_classification = None
        self.wt_aggregates = None
        self.mut_aggregates = None
        self.left_flank_wt_aggregates = None
        self.left_flank_mut_aggregates = None
        self.right_flank_wt_aggregates = None
        self.right_flank_mut_aggregates = None

    def set_proteins(self, wt_prot, mut_prot):
        self.wt_prot = wt_prot
        self.mut_prot = mut_prot

    def set_database_id(self, database_id):
        self.database_id = database_id

    def set_stored(self, boolean):
        self.stored = boolean

    def get_flank_size(self):
        flank_size = self.size // 2
        if flank_size % 2 != 0:
            flank_size += 1
        return flank_size

    def get_flanking_region_pos(self):
        if self.terminal_end == 'left':
            return None, self.right_flank_pos + 1
        elif self.terminal_end == 'right':
            return self.left_flank_pos - 1, None
        else:
            return self.left_flank_pos - 1, self.right_flank_pos + 1

    def aggregate(self, proteins, config):
        wt_pos_objs, mut_pos_objs = self.get_positions(proteins, config)

        wt_mappings, wt_disorder_score, wt_disorder_region = aggregate_positions(wt_pos_objs, config)
        self.wt_aggregates = wt_mappings.get_database_result_for_indel_aggregation(), wt_disorder_score, wt_disorder_region

        mut_mappings, mut_disorder_score, mut_disorder_region = aggregate_positions(mut_pos_objs, config)
        self.mut_aggregates = mut_mappings.get_database_result_for_indel_aggregation(), mut_disorder_score, mut_disorder_region

    def get_flank_positions(self, proteins):
        left_flank_wt_position_objs, left_flank_mut_position_objs = self.get_left_flank_positions(proteins)
        right_flank_wt_position_objs, right_flank_mut_position_objs = self.get_right_flank_positions(proteins)
        return left_flank_wt_position_objs, left_flank_mut_position_objs, right_flank_wt_position_objs, right_flank_mut_position_objs

    def flank_aggregates(self, proteins, config):
        left_flank_wt_position_objs, left_flank_mut_position_objs, right_flank_wt_position_objs, right_flank_mut_position_objs = self.get_flank_positions(proteins)

        left_flank_wt_mappings, left_flank_wt_disorder_score, left_flank_wt_disorder_region = aggregate_positions(left_flank_wt_position_objs, config)
        self.left_flank_wt_aggregates = left_flank_wt_mappings.get_database_result_for_indel_aggregation(), left_flank_wt_disorder_score, left_flank_wt_disorder_region

        left_flank_mut_mappings, left_flank_mut_disorder_score, left_flank_mut_disorder_region = aggregate_positions(left_flank_mut_position_objs, config)
        self.left_flank_mut_aggregates = left_flank_mut_mappings.get_database_result_for_indel_aggregation(), left_flank_mut_disorder_score, left_flank_mut_disorder_region

        right_flank_wt_mappings, right_flank_wt_disorder_score, right_flank_wt_disorder_region = aggregate_positions(right_flank_wt_position_objs, config)
        self.right_flank_wt_aggregates = right_flank_wt_mappings.get_database_result_for_indel_aggregation(), right_flank_wt_disorder_score, right_flank_wt_disorder_region

        right_flank_mut_mappings, right_flank_mut_disorder_score, right_flank_mut_disorder_region = aggregate_positions(right_flank_mut_position_objs, config)
        self.right_flank_mut_aggregates = right_flank_mut_mappings.get_database_result_for_indel_aggregation(), right_flank_mut_disorder_score, right_flank_mut_disorder_region

class Insertion(Indel):
    __slots__ = ['inserted_sequence']

    def __init__(self, left_flank=None, right_flank=None, inserted_sequence=None, tags=set(), raw_str=None):
        if raw_str is not None:
            flanks, inserted_sequence = raw_str.split('ins')
            if flanks.count('_') == 1:
                left_flank, right_flank = flanks.split('_')
            elif flanks[1:] == '1':
                left_flank = None
                right_flank = flanks
            else:
                left_flank = flanks
                right_flank = None
        super().__init__(left_flank=left_flank, right_flank=right_flank, tags=tags)
        self.inserted_sequence = inserted_sequence
        self.size = len(inserted_sequence)
        # if self.right_flank == None:
        #    self.right_flank_pos = 1
        #    self.left_flank_pos = 0
        if self.left_flank is None:
            self.left_flank_pos = self.right_flank_pos
            self.terminal = True
            if self.left_flank_pos == 1:
                self.terminal_end = 'left'
            else:
                self.terminal_end = 'right'

    def get_type(self):
        return 'Insertion', self.terminal, self.terminal_end

    def create_protein_name(self, wt_name, seq_name_length=5):
        if self.left_flank is not None:
            lfs = self.left_flank[1:]
        else:
            lfs = self.right_flank[1:]

        return '%s_ins_%s_%s' % (wt_name, lfs, self.inserted_sequence[:seq_name_length])

    def create_other_protein_name(self, wt_name, indel_protein_name):
        added_seq_len = len(indel_protein_name.split('_')[:-1])
        new_name = self.create_protein_name(wt_name, seq_name_length=(added_seq_len + 1))
        if new_name == indel_protein_name:
            return None
        return new_name

    def return_mutate_sequence(self, wt_seq):
        if self.right_flank_pos is not None:
            mutated_sequence = '%s%s%s' % (wt_seq[:(self.left_flank_pos)], self.inserted_sequence, wt_seq[(self.right_flank_pos - 1):])
        else:
            mutated_sequence = '%s%s' % (wt_seq[:(self.left_flank_pos)], self.inserted_sequence)
        return mutated_sequence

    def mutate_sequence(self, proteins):
        if proteins[self.wt_prot].nucleotide_sequence is not None:
            self.inserted_sequence = translate(self.inserted_sequence)
            if self.left_flank is not None:
                self.left_flank_pos = ((self.left_flank_pos - 1) // 3) + 1
                if self.right_flank is not None:
                    self.right_flank_pos = ((self.right_flank_pos - 1) // 3) + 1

        wt_seq = proteins[self.wt_prot].sequence

        mutated_sequence = self.return_mutate_sequence(wt_seq)
        proteins[self.mut_prot].sequence = mutated_sequence
        for (pos, aa) in enumerate(mutated_sequence):
            seq_pos = pos + 1
            position = position_package.Position(pos=seq_pos, wt_aa=aa, checked=True)
            proteins[self.mut_prot].positions[seq_pos] = position

    def get_positions(self, proteins, config):
        mut_pos_objects = []
        for pos in range(self.left_flank_pos + 1, self.left_flank_pos + len(self.inserted_sequence)):
            mut_pos_objects.append(proteins.protein_map[self.mut_prot].positions[pos])

        return ([], mut_pos_objects)

    def get_left_flank_positions(self, proteins):
        if self.terminal:
            if self.terminal_end == 'left':
                return None, None
        flank_size = self.get_flank_size()

        mut_pos_objects = []
        wt_pos_objects = []
        for pos in range((self.left_flank_pos - flank_size) + 1, self.left_flank_pos + 1):
            if pos <= 0:
                continue
            wt_pos_objects.append(proteins.protein_map[self.wt_prot].positions[pos])
            mut_pos_objects.append(proteins.protein_map[self.mut_prot].positions[pos])
        return wt_pos_objects, mut_pos_objects

    def get_right_flank_positions(self, proteins):
        if self.terminal:
            if self.terminal_end == 'right':
                return None, None
        flank_size = self.get_flank_size()

        mut_pos_objects = []
        wt_pos_objects = []
        for pos in range(self.right_flank_pos , self.right_flank_pos + flank_size):
            if pos not in proteins.protein_map[self.wt_prot].positions:
                continue
            wt_pos_objects.append(proteins.protein_map[self.wt_prot].positions[pos])

        for pos in range(self.right_flank_pos + self.size, self.right_flank_pos + self.size + flank_size):
            if pos not in proteins.protein_map[self.mut_prot].positions:
                continue
            mut_pos_objects.append(proteins.protein_map[self.mut_prot].positions[pos])
        return wt_pos_objects, mut_pos_objects

    def get_notation(self):
        if self.left_flank is not None:
            lfs = self.left_flank
        else:
            lfs = ''
        if self.right_flank is not None:
            if lfs != '':
                rfs = '_%s' % self.right_flank
            else:
                rfs = self.right_flank
        else:
            rfs = ''
        return '%s%sins%s' % (lfs, rfs, self.inserted_sequence)

    def get_scenario(self, n_wt, n_mut):
        if n_wt == 0 and n_mut == 0:
            return []
        if not self.terminal:
            if n_mut == 0:
                return ['FRA']
            if n_wt == 0:
                return ['FRA', 'InCl', 'CFRI']
            return ['FRA', 'InCl', 'CFRI', 'CFRA']
        else:
            if n_mut == 0:
                return ['FRC']
            if n_wt == 0:
                return ['FRC', 'InCl', 'CFRI']
            return ['FRC', 'InCl', 'CFRI', 'CFRA']

    def inclusion_check(self, wt_structure_annotations, mut_structure_annotations):
        del_list = []
        for struct_tuple in wt_structure_annotations:
            if self.terminal:
                if self.left_flank is None:  # Terminal deletion on the left side
                    if wt_structure_annotations[struct_tuple].is_covered(self.right_flank_pos):
                        if not wt_structure_annotations[struct_tuple].is_terminal(self.right_flank_pos):
                            del_list.append(struct_tuple)
                    else:
                        del_list.append(struct_tuple)
                else:  # Terminal insertion on the right side
                    if wt_structure_annotations[struct_tuple].is_covered(self.left_flank_pos):
                        if not wt_structure_annotations[struct_tuple].is_terminal(self.left_flank_pos):
                            del_list.append(struct_tuple)
                    else:
                        del_list.append(struct_tuple)
            else:
                # The wildtype structure should contain both flanks
                if wt_structure_annotations[struct_tuple].is_covered(self.left_flank_pos) and wt_structure_annotations[struct_tuple].is_covered(self.right_flank_pos):
                    # And they should be adjacent
                    if not wt_structure_annotations[struct_tuple].adjacency_check(self.left_flank_pos, self.right_flank_pos):
                        del_list.append(struct_tuple)
                else:
                    del_list.append(struct_tuple)

        for struct_tuple in del_list:
            del wt_structure_annotations[struct_tuple]

        # The flank positions are different for the mutant sequence! We take here the inner sides of the flanks.
        mut_left_flank_pos = self.left_flank_pos + 1
        mut_right_flank_pos = self.left_flank_pos + len(self.inserted_sequence) - 1  # They can be identical if the only AA is inserted
        del_list = []

        for struct_tuple in mut_structure_annotations:
            # The mutant structure must include the Mut flanks, attention for terminal insertions here
            if self.terminal:
                if self.left_flank is None:  # Terminal insertion on the left side
                    if not mut_structure_annotations[struct_tuple].is_covered(mut_right_flank_pos):
                        del_list.append(struct_tuple)

                else:  # Terminal insertion on the right side
                    if not mut_structure_annotations[struct_tuple].is_covered(mut_left_flank_pos):
                        del_list.append(struct_tuple)

            else:
                if not (mut_structure_annotations[struct_tuple].is_covered(mut_left_flank_pos) and mut_structure_annotations[struct_tuple].is_covered(mut_right_flank_pos)):
                    del_list.append(struct_tuple)

        for struct_tuple in del_list:
            del mut_structure_annotations[struct_tuple]

        return wt_structure_annotations, mut_structure_annotations


class Deletion(Indel):
    def __init__(self, left_flank=None, right_flank=None, tags=set(), raw_str=None):
        if raw_str is not None:
            flanks = raw_str.split('del')[0]
            if flanks.count('_') == 1:
                left_flank, right_flank = flanks.split('_')
            else:
                left_flank = flanks
                right_flank = left_flank
        super().__init__(left_flank=left_flank, right_flank=right_flank, tags=tags)
        if self.left_flank_pos == 1:
            self.terminal = True
            self.terminal_end = 'left'
        self.size = (self.right_flank_pos - self.left_flank_pos) + 1

    def get_type(self):
        return 'Deletion', self.terminal, self.terminal_end

    def create_protein_name(self, wt_name):
        return '%s_del_%s_%s' % (wt_name, self.left_flank[1:], self.right_flank[1:])

    def create_other_protein_name(self, wt_name, indel_protein_name):
        return None

    def return_mutate_sequence(self, wt_seq):
        mutated_sequence = '%s%s' % (wt_seq[:(self.left_flank_pos - 1)], wt_seq[self.right_flank_pos:])
        return mutated_sequence

    def mutate_sequence(self, proteins):
        if proteins[self.wt_prot].nucleotide_sequence is not None:
            self.left_flank_pos = ((self.left_flank_pos - 1) // 3) + 1
            self.right_flank_pos = ((self.right_flank_pos - 1) // 3) + 1
        wt_seq = proteins[self.wt_prot].sequence
        if self.right_flank_pos == len(wt_seq):
            self.terminal = True
            self.terminal_end = 'right'
        mutated_sequence = self.return_mutate_sequence(wt_seq)
        proteins[self.mut_prot].sequence = mutated_sequence
        for (pos, aa) in enumerate(mutated_sequence):
            seq_pos = pos + 1
            position = position_package.Position(pos=seq_pos, wt_aa=aa, checked=True)
            proteins[self.mut_prot].positions[seq_pos] = position

    def get_positions(self, proteins, config):
        pos_objects = []
        for pos in range(self.left_flank_pos, self.right_flank_pos + 1):
            try:
                pos_objects.append(proteins.protein_map[self.wt_prot].positions[pos])
            except:
                if config is not None:
                    if pos not in proteins.protein_map[self.wt_prot].positions:
                        config.errorlog.add_warning(f'Position {pos} not in {self.wt_prot}')
                continue
        return (pos_objects, [])

    def get_left_flank_positions(self, proteins):
        if self.terminal:
            if self.terminal_end == 'left':
                return None, None
        flank_size = self.get_flank_size()

        mut_pos_objects = []
        wt_pos_objects = []
        for pos in range((self.left_flank_pos - flank_size), self.left_flank_pos):
            if pos <= 0:
                continue
            wt_pos_objects.append(proteins.protein_map[self.wt_prot].positions[pos])
            mut_pos_objects.append(proteins.protein_map[self.mut_prot].positions[pos])
        return wt_pos_objects, mut_pos_objects

    def get_right_flank_positions(self, proteins):
        if self.terminal:
            if self.terminal_end == 'right':
                return None, None
        flank_size = self.get_flank_size()

        mut_pos_objects = []
        wt_pos_objects = []
        for pos in range(self.right_flank_pos + 1, self.right_flank_pos + flank_size + 1):
            if pos not in proteins.protein_map[self.wt_prot].positions:
                continue
            wt_pos_objects.append(proteins.protein_map[self.wt_prot].positions[pos])

        for pos in range(self.left_flank_pos + 1, self.left_flank_pos + 1 + flank_size):
            if pos not in proteins.protein_map[self.mut_prot].positions:
                continue
            mut_pos_objects.append(proteins.protein_map[self.mut_prot].positions[pos])
        return wt_pos_objects, mut_pos_objects

    def get_notation(self):
        return '%s_%sdel' % (self.left_flank, self.right_flank)

    def get_scenario(self, n_wt, n_mut):
        if n_wt == 0 and n_mut == 0:
            return []
        if not self.terminal:
            if n_mut == 0:
                return ['FRA', 'InCl', 'CFRI']
            if n_wt == 0:
                return ['FRA']
            return ['FRA', 'InCl', 'CFRI', 'CFRA']
        else:
            if n_mut == 0:
                return ['FRC', 'InCl', 'CFRI']
            if n_wt == 0:
                return ['FRC']
            return ['FRC', 'InCl', 'CFRI', 'CFRA']

    def inclusion_check(self, wt_structure_annotations, mut_structure_annotations):
        del_list = []
        mut_del_list = []
        for struct_tuple in wt_structure_annotations:
            # The wildtype structure should contain both flanks
            if wt_structure_annotations[struct_tuple].is_covered(self.left_flank_pos) and wt_structure_annotations[struct_tuple].is_covered(self.right_flank_pos):
                # In this case we can delete it from the mutation_annotations, if it is there
                if struct_tuple in mut_structure_annotations:
                    mut_del_list.append(struct_tuple)
            else:
                del_list.append(struct_tuple)
        for struct_tuple in del_list:
            del wt_structure_annotations[struct_tuple]

        for struct_tuple in mut_del_list:
            del mut_structure_annotations[struct_tuple]

        # The flank positions are different for the mutant sequence!
        mut_left_flank_pos = self.left_flank_pos - 1
        mut_right_flank_pos = self.left_flank_pos

        del_list = []
        for struct_tuple in mut_structure_annotations:
            # The mutant structure must not include the WT flanks, but must include the mut flanks and they have to be adjacent, attention for terminal deletions here
            if self.terminal:
                if mut_right_flank_pos == 1:  # Terminal deletion on the left side
                    if mut_structure_annotations[struct_tuple].is_covered(mut_right_flank_pos):
                        if not mut_structure_annotations[struct_tuple].is_terminal(mut_right_flank_pos):
                            del_list.append(struct_tuple)
                    else:
                        del_list.append(struct_tuple)
                else:  # Terminal deletion on the right side
                    if mut_structure_annotations[struct_tuple].is_covered(mut_left_flank_pos):
                        if not mut_structure_annotations[struct_tuple].is_terminal(mut_left_flank_pos):
                            del_list.append(struct_tuple)
                    else:
                        del_list.append(struct_tuple)

            else:
                if mut_structure_annotations[struct_tuple].is_covered(mut_left_flank_pos) and mut_structure_annotations[struct_tuple].is_covered(mut_right_flank_pos):
                    if mut_structure_annotations[struct_tuple].adjacency_check(mut_left_flank_pos, mut_right_flank_pos):
                        continue  # In this case do not delete annotation
                    else:
                        del_list.append(struct_tuple)
                else:
                    del_list.append(struct_tuple)

        for struct_tuple in del_list:
            del mut_structure_annotations[struct_tuple]

        return wt_structure_annotations, mut_structure_annotations


class Substitution(Indel):
    __slots__ = ['inserted_sequence']

    def __init__(self, left_flank=None, right_flank=None, inserted_sequence=None, tags=set(), raw_str=None):
        if raw_str is not None:
            flanks, inserted_sequence = raw_str.split('delins')
            if flanks.count('_') == 1:
                left_flank, right_flank = flanks.split('_')
            else:
                left_flank = flanks
                right_flank = left_flank
        super().__init__(left_flank=left_flank, right_flank=right_flank, tags=tags)
        self.inserted_sequence = inserted_sequence
        if self.left_flank_pos == 1:
            self.terminal = True
            self.terminal_end = 'left'

        self.size = max([(self.right_flank_pos - self.left_flank_pos) + 1, len(inserted_sequence)])

    def get_type(self):
        return 'Substitution', self.terminal, self.terminal_end

    def create_protein_name(self, wt_name, seq_name_length=5):
        return '%s_delins_%s_%s_%s' % (wt_name, self.left_flank[1:], self.right_flank[1:], self.inserted_sequence[:seq_name_length])

    def create_other_protein_name(self, wt_name, indel_protein_name):
        added_seq_len = len(indel_protein_name.split('_')[:-1])
        new_name = self.create_protein_name(wt_name, seq_name_length=(added_seq_len + 1))
        if new_name == indel_protein_name:
            return None
        return new_name

    def return_mutate_sequence(self, wt_seq):
        if self.right_flank_pos == len(wt_seq):
            self.terminal = True
            self.terminal_end = 'right'
        mutated_sequence = '%s%s%s' % (wt_seq[:(self.left_flank_pos - 1)], self.inserted_sequence, wt_seq[(self.right_flank_pos):])

        return mutated_sequence

    def mutate_sequence(self, proteins):
        if proteins[self.wt_prot].nucleotide_sequence is not None:
            self.inserted_sequence = translate(self.inserted_sequence)
            self.left_flank_pos = ((self.left_flank_pos - 1) // 3) + 1
            self.right_flank_pos = ((self.right_flank_pos - 1) // 3) + 1

        wt_seq = proteins[self.wt_prot].sequence
        mutated_sequence = self.return_mutate_sequence(wt_seq)

        proteins[self.mut_prot].sequence = mutated_sequence
        for (pos, aa) in enumerate(mutated_sequence):
            seq_pos = pos + 1
            position = position_package.Position(pos=seq_pos, wt_aa=aa, checked=True)
            proteins[self.mut_prot].positions[seq_pos] = position

    def get_positions(self, proteins, config):
        wt_pos_objects = []
        for pos in range(self.left_flank_pos, self.right_flank_pos + 1):
            wt_pos_objects.append(proteins.protein_map[self.wt_prot].positions[pos])

        mut_pos_objects = []
        for pos in range(self.left_flank_pos + 1, self.left_flank_pos + len(self.inserted_sequence)):
            mut_pos_objects.append(proteins.protein_map[self.mut_prot].positions[pos])

        return (wt_pos_objects, mut_pos_objects)

    def get_left_flank_positions(self, proteins):
        if self.terminal:
            if self.terminal_end == 'left':
                return None, None
        flank_size = self.get_flank_size()

        mut_pos_objects = []
        wt_pos_objects = []
        for pos in range((self.left_flank_pos - flank_size), self.left_flank_pos):
            if pos <= 0:
                continue
            wt_pos_objects.append(proteins.protein_map[self.wt_prot].positions[pos])
            mut_pos_objects.append(proteins.protein_map[self.mut_prot].positions[pos])
        return wt_pos_objects, mut_pos_objects

    def get_right_flank_positions(self, proteins):
        if self.terminal:
            if self.terminal_end == 'right':
                return None, None
        flank_size = self.get_flank_size()

        mut_pos_objects = []
        wt_pos_objects = []
        for pos in range(self.right_flank_pos + 1, self.right_flank_pos + flank_size + 1):
            if pos not in proteins.protein_map[self.wt_prot].positions:
                continue
            wt_pos_objects.append(proteins.protein_map[self.wt_prot].positions[pos])

        for pos in range(self.left_flank_pos +len(self.inserted_sequence) + 1, self.left_flank_pos + len(self.inserted_sequence) + 1 + flank_size):
            if pos not in proteins.protein_map[self.mut_prot].positions:
                continue
            mut_pos_objects.append(proteins.protein_map[self.mut_prot].positions[pos])
        return wt_pos_objects, mut_pos_objects

    def get_notation(self):
        return '%s_%sdelins%s' % (self.left_flank, self.right_flank, self.inserted_sequence)

    def get_scenario(self, n_wt, n_mut):
        if n_wt == 0 and n_mut == 0:
            return []
        if not self.terminal:
            if n_mut == 0:
                return ['FRA', 'InCl', 'CFRI']
            if n_wt == 0:
                return ['FRA', 'InCl', 'CFRI']
            return ['FRA', 'InCl', 'CFRI', 'CFRA']
        else:
            if n_mut == 0:
                return ['FRC', 'InCl', 'CFRI']
            if n_wt == 0:
                return ['FRC', 'InCl', 'CFRI']
            return ['FRC', 'InCl', 'CFRI', 'CFRA']

    def inclusion_check(self, wt_structure_annotations, mut_structure_annotations):
        # This just checks if there is a region resolved in the structure of appropriate size, not if it is correct sequence. This should be checked elsewhere with the help of sequence identies
        del_list = []
        for struct_tuple in wt_structure_annotations:
            # The wildtype structure should contain both flanks (terminal or not), similar to deletions
            if not (wt_structure_annotations[struct_tuple].is_covered(self.left_flank_pos) and wt_structure_annotations[struct_tuple].is_covered(self.right_flank_pos)):
                del_list.append(struct_tuple)

        for struct_tuple in del_list:
            del wt_structure_annotations[struct_tuple]

        # The flank positions are different for the mutant sequence! We take here the inner sides of the (mutant) flanks.
        mut_left_flank_pos = self.left_flank_pos
        mut_right_flank_pos = self.left_flank_pos + len(self.inserted_sequence) - 1
        del_list = []

        for struct_tuple in mut_structure_annotations:
            # The mutant structure is here to be treated as the wt structure
            if not (mut_structure_annotations[struct_tuple].is_covered(mut_left_flank_pos) and mut_structure_annotations[struct_tuple].is_covered(mut_right_flank_pos)):
                del_list.append(struct_tuple)

        for struct_tuple in del_list:
            del mut_structure_annotations[struct_tuple]

        return wt_structure_annotations, mut_structure_annotations
