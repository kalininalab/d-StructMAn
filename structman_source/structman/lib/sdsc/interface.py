from structman.lib.sdsc.sdsc_utils import doomsday_protocol

#Called by serializedPipeline

def calculate_aggregated_interface_map(config, protein_id, interface_map, backmaps, quality_measures, mapped_positions):
    position_interface_map = {}
    aggregated_interfaces = []
    if config.verbosity >= 5:
        print('Call of calculate_aggregated_interface_map:', protein_id)
    for structure_id, t_chain in interface_map:
        interfaces = interface_map[(structure_id, t_chain)]
        backmap = backmaps[(structure_id, t_chain)]

        if config.verbosity >= 6:
            print('\ncalculate_aggregated_interface_map:', protein_id, structure_id, t_chain, interfaces.keys())
        if interfaces is None:
            config.errorlog.add_error(f'Interfaces is None: {protein_id} {structure_id} {t_chain}')
        #    continue
        for (chain, i_chain) in interfaces:
            if chain != t_chain:
                continue

            #print('interface:', chain, i_chain)
            #print(position_interface_map)

            other_interfaces = {}
            back_mappable_positions = 0
            for res in interfaces[(chain, i_chain)].residues:
                if res not in backmap:
                    continue
                position = backmap[res]
                if position is None:
                    continue
                back_mappable_positions += 1
                if position in position_interface_map:
                    known_interfaces = position_interface_map[position]
                    for known_interface in known_interfaces:
                        if not known_interface in other_interfaces:
                            other_interfaces[known_interface] = 1
                        else:
                            other_interfaces[known_interface] += 1
                else:
                    position_interface_map[position] = set()
                    #potential_new_interface = len(aggregated_interfaces)

            if back_mappable_positions == 0:
                continue

            #If the currently checked interface overlaps more than half with a known interface, fuse them
            no_fusion_happened = True
            fused_with = []
            for known_interface_number in other_interfaces:
                overlap = other_interfaces[known_interface_number]

                #print(known_interface_number, chain, i_chain, overlap, back_mappable_positions, len(aggregated_interfaces[known_interface_number].positions), f'old recommended_complex: {aggregated_interfaces[known_interface_number].recommended_complex} {aggregated_interfaces[known_interface_number].chain} {aggregated_interfaces[known_interface_number].interacting_chain}')

                if overlap/back_mappable_positions > 0.5 or overlap/len(aggregated_interfaces[known_interface_number].positions) > 0.5:
                    no_fusion_happened = False
                    seq_id, cov = quality_measures[(protein_id, structure_id, chain)]
                    interface_annotation_score = calc_interface_annotation_score(seq_id, cov, len(interfaces[(chain, i_chain)].residues))

                    for res in interfaces[(chain, i_chain)].residues:
                        if res not in backmap:
                            #print(res, 'not in backmap',protein_id, structure_id, chain, i_chain)
                            continue
                        position = backmap[res]
                        if position is None:
                            #print('position is none',protein_id, structure_id, chain, i_chain)
                            continue
                        position_interface_map[position].add(known_interface_number)
                        if (not position in aggregated_interfaces[known_interface_number].positions) or interface_annotation_score > aggregated_interfaces[known_interface_number].interface_annotation_score:
                            aggregated_interfaces[known_interface_number].add_position(position, (structure_id, chain, res))
                        if res in interfaces[(chain, i_chain)].interactions:
                            i_res = interfaces[(chain, i_chain)].interactions[res]
                            if not (structure_id, i_chain) in mapped_positions:
                                #print(structure_id, i_chain, 'not in structures', protein_id, chain, res, i_res)
                                continue
                            if i_res not in mapped_positions[(structure_id, i_chain)]: #Happens for ligands part of chain in structure file
                                continue
                            i_res_mapped_pos = mapped_positions[(structure_id, i_chain)][i_res]
                            for i_protein in i_res_mapped_pos:
                                if position not in aggregated_interfaces[known_interface_number].pos_pos_interactions:
                                    pos_pos_i = Position_Position_Interaction(protein_id, position, i_protein, i_res_mapped_pos[i_protein])
                                    aggregated_interfaces[known_interface_number].pos_pos_interactions[position] = pos_pos_i
                                    aggregated_interfaces[known_interface_number].pos_pos_interactions[position].set_recommended_complex(structure_id, (chain, res, i_chain, i_res))
                                elif interface_annotation_score > aggregated_interfaces[known_interface_number].interface_annotation_score:
                                    aggregated_interfaces[known_interface_number].pos_pos_interactions[position].set_recommended_complex(structure_id, (chain, res, i_chain, i_res))

                    if config.verbosity >= 6:
                        print(f'Fusion happend: {protein_id} {known_interface_number} {structure_id} {chain} {i_chain}, old recommended_complex: {aggregated_interfaces[known_interface_number].recommended_complex} {aggregated_interfaces[known_interface_number].chain} {aggregated_interfaces[known_interface_number].interacting_chain}')

                    fused_with.append(known_interface_number)
                    if interface_annotation_score > aggregated_interfaces[known_interface_number].interface_annotation_score:
                        aggregated_interfaces[known_interface_number].set_recommended_complex(structure_id, chain, i_chain, interface_annotation_score)
                        if config.verbosity >= 5:
                            print(f'Overwrite of recommended_complex: {structure_id}, {chain}, {i_chain}')

            if len(fused_with) > 1:
                fused_with = sorted(fused_with)
                if config.verbosity >= 6:
                    print(f'Superfusion of {protein_id} {fused_with}')
                interface_number_1 = fused_with[0]
                to_be_removed = set()
                for interface_number_2 in fused_with[1:]:
                    aggregated_interfaces[interface_number_1] = aggregated_interfaces[interface_number_1].fuse(aggregated_interfaces[interface_number_2])
                    to_be_removed.add(interface_number_2)

                new_aggregated_interfaces = []
                new_position_interface_map = {}

                diff_count = 0
                interface_number_map = {}
                for interface_number, agg_int in enumerate(aggregated_interfaces):
                    if not interface_number in to_be_removed:
                        new_aggregated_interfaces.append(agg_int)
                        interface_number_map[interface_number] = interface_number - diff_count
                    else:
                        diff_count += 1
                        interface_number_map[interface_number] = interface_number - diff_count

                for position in position_interface_map:
                    new_position_interface_map[position] = set()
                    for known_interface in position_interface_map[position]:
                        new_position_interface_map[position].add(interface_number_map[known_interface])

                aggregated_interfaces = new_aggregated_interfaces
                position_interface_map = new_position_interface_map

            if no_fusion_happened:
                trigger = True
                for res in interfaces[(chain, i_chain)].residues:
                    if not res in backmap:
                        #print(res, '2 not in backmap',protein_id, structure_id, chain, i_chain)
                        continue
                    position = backmap[res]
                    if position is None:
                        #print('2 position is none',protein_id, structure_id, chain, i_chain)
                        continue

                    if trigger: #We need at least one position being backmapped
                        trigger = False
                        seq_id, cov = quality_measures[(protein_id, structure_id, chain)]
                        interface_annotation_score = calc_interface_annotation_score(seq_id, cov, len(interfaces[(chain, i_chain)].residues))
                        aggregated_interfaces.append(Aggregated_interface(protein_id, chain = chain, interacting_chain = i_chain, recommended_complex = structure_id, interface_annotation_score = interface_annotation_score))
                        if config.verbosity >= 6:
                            print(f'Creation of new Aggregated_interface: {protein_id} {len(aggregated_interfaces) - 1} {chain} {i_chain} {structure_id}')

                    aggregated_interfaces[-1].add_position(position, (structure_id, chain, res))
                    position_interface_map[position].add(len(aggregated_interfaces) - 1) 
                    if res in interfaces[(chain, i_chain)].interactions:
                        i_res = interfaces[(chain, i_chain)].interactions[res]
                        if not (structure_id, i_chain) in mapped_positions:
                            #print(structure_id, i_chain, 'not in structures 2', protein_id, chain, res, i_res)
                            continue
                        if config.verbosity >= 6:
                            #print(chain, i_chain, interfaces[(chain, i_chain)].interactions)
                            #print(interfaces[(chain, i_chain)].residues)
                            print(f'Mapped positions for: {structure_id} {i_chain}: {mapped_positions[(structure_id, i_chain)]}')
                        if i_res not in mapped_positions[(structure_id, i_chain)]: #Happens for ligands part of chain in structure file
                            continue
                        i_res_mapped_pos = mapped_positions[(structure_id, i_chain)][i_res]
                        for i_protein in i_res_mapped_pos:
                            pos_pos_i = Position_Position_Interaction(protein_id, position, i_protein, i_res_mapped_pos[i_protein])
                            aggregated_interfaces[-1].pos_pos_interactions[position] = pos_pos_i
                            aggregated_interfaces[-1].pos_pos_interactions[position].set_recommended_complex(structure_id, (chain, res, i_chain, i_res))

    return aggregated_interfaces

def calc_interface_annotation_score(seq_id, coverage, interface_coverage):
    return interface_coverage * seq_id * coverage

class Aggregated_interface:
    __slots__ = ['protein',  'recommended_complex', 'chain', 'interacting_chain', 'positions', 'interface_annotation_score', 'pos_pos_interactions', 'database_id']
    def __init__(self, protein_id, recommended_complex = None, chain = None, interacting_chain = None, positions = None, interface_annotation_score = None, database_id = None):
        self.protein = protein_id
        self.recommended_complex = recommended_complex
        self.chain = chain
        self.interacting_chain = interacting_chain
        self.interface_annotation_score = interface_annotation_score
        if positions is None:
            self.positions = {}
        else:
            self.positions = positions
        self.pos_pos_interactions = {}
        self.database_id = database_id

    def deconstruct(self):
        del self.positions
        for pos in self.pos_pos_interactions:
            self.pos_pos_interactions[pos].deconstruct()
        del self.pos_pos_interactions
        doomsday_protocol(self)

    def add_position(self, position, recommended_residue):
        self.positions[position] = recommended_residue

    def set_recommended_complex(self, recommended_complex, chain, interacting_chain, interface_annotation_score):
        self.recommended_complex = recommended_complex
        self.chain = chain
        self.interacting_chain = interacting_chain
        self.interface_annotation_score = interface_annotation_score

    def fuse(self, other_interface):
        if other_interface.interface_annotation_score > self.interface_annotation_score:
            return other_interface.fuse(self)
        for pos in other_interface.positions:
            self.positions[pos] = other_interface.positions[pos]
        for pos in other_interface.pos_pos_interactions:
            self.pos_pos_interactions[pos] = other_interface.pos_pos_interactions[pos]
        return self

class Position_Position_Interaction:
    __slots__ = ['protein_a', 'position_a', 'protein_b', 'position_b', 'recommended_complex', 'recommended_interaction', 'database_id']
    def __init__(self, protein_a, position_a, protein_b, position_b, database_id = None):
        self.protein_a = protein_a
        self.position_a = position_a
        self.protein_b = protein_b
        self.position_b = position_b
        self.database_id = database_id

    def deconstruct(self):
        del self.protein_a
        del self.position_a
        del self.protein_b
        del self.position_b
        doomsday_protocol(self)

    def set_recommended_complex(self, recommended_complex, recommended_interaction):
        self.recommended_complex = recommended_complex
        self.recommended_interaction = recommended_interaction



class Interface:
    __slots__ = ['chain', 'interacting_chain', 'residues', 'support_residues', 'interactions', 'stored']

    def __init__(self, chain, interacting_chain, stored = False):
        self.chain = chain
        self.interacting_chain = interacting_chain
        self.interactions = {}
        self.residues = set()
        self.support_residues = set()
        self.stored = stored

    def deconstruct(self):
        del self.residues
        del self.interactions
        del self.support_residues
        doomsday_protocol(self)

    def add_interaction(self, res_a, res_b):
        self.interactions[res_a] = res_b
        self.residues.add(res_a)

    def add_support(self, res):
        self.support_residues.add(res)
        self.residues.add(res)
