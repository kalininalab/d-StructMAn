import os
import sys
import time
from ast import literal_eval
import ray
import structman
from structman import settings
from structman.lib import database, rin, sdsc
from structman.lib.output.feature import init_feature_table
from structman.lib.output.output import OutputGenerator
from structman.lib.output.out_utils import makeViolins
from filelock import FileLock

def ray_init(config):
    if ray.is_initialized():
        if config.verbosity >= 2:
            print('Ray_init was called, but ray was already initialised')
        return
    os.environ["PYTHONPATH"] = f'{settings.ROOT_DIR}:{os.environ.get("PYTHONPATH", "")}'
    os.environ["PYTHONPATH"] = f'{settings.LIB_DIR}:{os.environ.get("PYTHONPATH", "")}'
    os.environ["PYTHONPATH"] = f'{settings.RINERATOR_DIR}:{os.environ.get("PYTHONPATH", "")}'
    os.environ["PYTHONPATH"] = f'{settings.OUTPUT_DIR}:{os.environ.get("PYTHONPATH", "")}'
    if config.iupred_path != '':
        os.environ["PYTHONPATH"] = f'{os.path.abspath(os.path.realpath(config.iupred_path))}:{os.environ.get("PYTHONPATH", "")}'

    logging_level = 20
    if config.verbosity <= 1:
        logging_level = 0
    ray.init(num_cpus=config.proc_n, include_dashboard=False, ignore_reinit_error=True, logging_level = logging_level)

def writeStatFile(out_file, mutation_dict, class_dict, tag_map, stat_dict=None):
    seq_id_threshold = 0.99
    startline = 'Tag\tTotal proteins\tTotal positions\tUnmapped proteins\tEntirely disordered proteins\tProteins mapped to at least one corresponding structure (seq-id > %s%%\tProteins mapped only to structure of homologs (seq-id <= %s%%)\tMapped positions\tMapped into at least one corresponding structure (seq-id > %s%%)\tMapped only in homologs (seq-id <= %s%%)\tUnmapped, Disorder\tUnmapped, Globular' % (str(seq_id_threshold), str(seq_id_threshold), str(seq_id_threshold), str(seq_id_threshold))
    outmap = {'All': [{}, 0, 0, 0, 0, 0, 0]}

    if stat_dict is not None:
        class_dict = stat_dict
        max_seq_key = 1
    else:
        max_seq_key = 21
    for m in tag_map:
        if tag_map[m] is not None and tag_map[m] != '':
            raw_tags = tag_map[m].split(',')
        else:
            raw_tags = []
        tags = set(['All'])
        for tag in raw_tags:
            if tag[0] == '#':
                tag = tag[1:].split(':')[0]
            tags.add(tag)

        for tag in tags:
            if tag not in outmap:
                outmap[tag] = [{}, 0, 0, 0, 0, 0, 0]
            g = mutation_dict[m][1]
            if g not in outmap[tag][0]:
                outmap[tag][0][g] = 1
            outmap[tag][1] += 1

            if m in class_dict:
                clas = class_dict[m][0]
                max_seq_id = class_dict[m][max_seq_key]
                if clas != 'Disorder' and clas is not None:
                    outmap[tag][2] += 1
                    if max_seq_id > seq_id_threshold:
                        outmap[tag][5] += 1
                        outmap[tag][0][g] = 2
                    else:
                        outmap[tag][6] += 1
                        if not outmap[tag][0][g] == 2:
                            outmap[tag][0][g] = 3
                elif clas == 'Disorder':
                    outmap[tag][3] += 1
                else:
                    outmap[tag][4] += 1
                    if not outmap[tag][0][g] > 1:
                        outmap[tag][0][g] = 0
            else:
                outmap[tag][4] += 1
                if not outmap[tag][0][g] > 1:
                    outmap[tag][0][g] = 0
    if None in outmap:
        del outmap[None]

    lines = [startline]
    for tag in outmap:
        tot_prot = len(outmap[tag][0])
        tot_pos = outmap[tag][1]
        mapped = outmap[tag][2]
        dis = outmap[tag][3]
        unmapped = outmap[tag][4]
        mapped_to_corr = outmap[tag][5]
        mapped_to_homolog = outmap[tag][6]

        prot_numbers = [0, 0, 0, 0]
        for g in outmap[tag][0]:
            prot_numbers[outmap[tag][0][g]] += 1

        if float(tot_pos) == 0.0:
            continue
        line = '%s\t%s\t%s\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)' % (
            tag,
            str(tot_prot),
            str(tot_pos),
            str(prot_numbers[0]), str(100. * float(prot_numbers[0]) / float(tot_prot)),
            str(prot_numbers[1]), str(100. * float(prot_numbers[1]) / float(tot_prot)),
            str(prot_numbers[2]), str(100. * float(prot_numbers[2]) / float(tot_prot)),
            str(prot_numbers[3]), str(100. * float(prot_numbers[3]) / float(tot_prot)),
            str(mapped), str(100. * float(mapped) / float(tot_pos)),
            str(mapped_to_corr), str(100. * float(mapped_to_corr) / float(tot_pos)),
            str(mapped_to_homolog), str(100. * float(mapped_to_homolog) / float(tot_pos)),
            str(dis), str(100. * float(dis) / float(tot_pos)),
            str(unmapped), str(100. * float(unmapped) / float(tot_pos)))
        lines.append(line)
    f = open(out_file, 'w')
    f.write("\n".join(lines))
    f.close()


def init_classification_table(class_file, obj_only = False):
    classification_output = OutputGenerator()
    headers = [
        'Input Protein ID', 'Primary Protein ID', 'Uniprot-Ac', 'Uniprot ID', 'Refseq', 'Residue ID', 'Amino Acid', 'Position', 'Tags',
        'Weighted Location', 'Weighted Mainchain Location', 'Weighted Sidechain Location',
        'Class', 'Simple Class', 'RIN Class', 'RIN Simple Class', 'Individual Interactions',
        'Confidence Value', 'Secondary Structure', 'Recommended Structure', 'Sequence-ID', 'Coverage', 'Resolution',
        'Max Seq Id Structure', 'Max Sequence-ID', 'Max Seq Id Coverage', 'Max Seq Id Resolution', 'Amount of mapped structures'
    ]


    classification_output.add_headers(headers)
    if obj_only:
        return classification_output

    f = open(class_file, 'a')
    f.write(classification_output.get_header())

    return classification_output, f


@ray.remote(max_calls = 1)
def unpack_rows(package, store_data, tag_map, snv_map, snv_tag_map, protein_dict):
    from structman.utils import unpack
    config, input_res_id_dict, classification_file_path, feature_file_path = store_data

    classification_output = init_classification_table(None, obj_only = True)
    feature_output = init_feature_table(None, obj_only = True)

    results = []
    c_lines = []
    f_lines = []
    stat_dict = {}

    max_lines = 10000

    for row in package:
        m = row[0]
        position_number = row[1]
        wt_aa = row[2]
        prot_db_id = row[3]
        iupred = row[4]
        disorder_state = row[5]
        if row[6] is not None:
            (recommended_structure_str, max_seq_structure_str) = unpack(row[6])
        else:
            recommended_structure_str = None
            max_seq_structure_str = None

        recommended_structure, seq_id, cov, resolution = sdsc.process_recommend_structure_str(recommended_structure_str)
        max_seq_structure, max_seq_seq_id, max_seq_cov, max_seq_resolution = sdsc.process_recommend_structure_str(max_seq_structure_str)

        if row[7] is None:
            (weighted_sc, mainchain_location, sidechain_location, rsa, mainchain_rsa, sidechain_rsa, Class, rin_class, simple_class, rin_simple_class,
                interaction_str, conf, mv_sec_ass, amount_of_structures, rin_profile, modres_score, b_factor, centrality_scores, weighted_phi,
                weighted_psi, intra_ssbond_prop, inter_ssbond_prop, intra_link_prop, inter_link_prop, cis_conformation_prop, cis_follower_prop,
                inter_chain_median_kd, inter_chain_dist_weighted_kd, inter_chain_median_rsa, inter_chain_dist_weighted_rsa,
                intra_chain_median_kd, intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                weighted_inter_chain_interactions_median, weighted_inter_chain_interactions_dist_weighted,
                weighted_intra_chain_interactions_median, weighted_intra_chain_interactions_dist_weighted) = [None]*38

            centrality_scores = rin.Centrality_scores()
            rin_profile = rin.Interaction_profile()

        else:
            try:
                (weighted_sc, mainchain_location, sidechain_location, rsa, mainchain_rsa, sidechain_rsa, Class, rin_class, simple_class, rin_simple_class,
                interaction_str, conf, mv_sec_ass, amount_of_structures, rin_profile_str, modres_score, b_factor, weighted_centrality_scores, weighted_phi,
                weighted_psi, intra_ssbond_prop, inter_ssbond_prop, intra_link_prop, inter_link_prop, cis_conformation_prop, cis_follower_prop,
                inter_chain_median_kd, inter_chain_dist_weighted_kd, inter_chain_median_rsa, inter_chain_dist_weighted_rsa,
                intra_chain_median_kd, intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                weighted_inter_chain_interactions_median, weighted_inter_chain_interactions_dist_weighted,
                weighted_intra_chain_interactions_median, weighted_intra_chain_interactions_dist_weighted) = unpack(row[7])

                centrality_scores = rin.Centrality_scores(code_str=weighted_centrality_scores)
                rin_profile = rin.Interaction_profile(profile_str=rin_profile_str)

            except:
                print('Error in retrieving position data, couldnt unpack:', m)
                continue

        input_res_id = input_res_id_dict[m]
        if input_res_id is None:
            input_res_id = ''

        if disorder_state is not None:
            if Class is None and (disorder_state == 'disorder' or disorder_state[0] == 'D'):
                Class = 'Disorder'
                simple_class = 'Disorder'
                rin_class = 'Disorder'
                rin_simple_class = 'Disorder'

        (prot_id, u_ac, refseq, u_id, error_code, error, input_id) = protein_dict[prot_db_id]

        if max_seq_seq_id == '-':
            stat_dict[m] = (Class, 0.)
        else:
            stat_dict[m] = (Class, float(max_seq_seq_id))

        if config.skipref:
            refseq = ''
        classification_output.add_value('Input Protein ID', input_id)
        classification_output.add_value('Primary Protein ID', prot_id)
        classification_output.add_value('Uniprot-Ac', u_ac)
        classification_output.add_value('Uniprot ID', u_id)
        classification_output.add_value('Refseq', refseq)
        classification_output.add_value('Residue ID', input_res_id)
        classification_output.add_value('Amino Acid', wt_aa)
        classification_output.add_value('Position', position_number)
        classification_output.add_value('Tags', tag_map[m])
        classification_output.add_value('Weighted Location', weighted_sc)
        classification_output.add_value('Weighted Mainchain Location', mainchain_location)
        classification_output.add_value('Weighted Sidechain Location', sidechain_location)
        classification_output.add_value('Class', Class)
        classification_output.add_value('Simple Class', simple_class)
        classification_output.add_value('RIN Class', rin_class)
        classification_output.add_value('RIN Simple Class', rin_simple_class)
        classification_output.add_value('Individual Interactions', interaction_str)
        classification_output.add_value('Confidence Value', conf)
        classification_output.add_value('Secondary Structure', mv_sec_ass)
        classification_output.add_value('Recommended Structure', recommended_structure)
        classification_output.add_value('Sequence-ID', seq_id)
        classification_output.add_value('Coverage', cov)
        classification_output.add_value('Resolution', resolution)
        classification_output.add_value('Max Seq Id Structure', max_seq_structure)
        classification_output.add_value('Max Sequence-ID', max_seq_seq_id)
        classification_output.add_value('Max Seq Id Coverage', max_seq_cov)
        classification_output.add_value('Max Seq Id Resolution', max_seq_resolution)
        classification_output.add_value('Amount of mapped structures', amount_of_structures)

        c_lines.append(classification_output.pop_line())

        if m not in snv_map:
            snv_map[m] = {0: wt_aa}
            snv_tag_map[0] = tag_map[m]
        for snv_database_id in snv_map[m]:
            new_aa = snv_map[m][snv_database_id]
            tags = snv_tag_map[snv_database_id]
            aac = "%s%s%s" % (wt_aa, str(position_number), new_aa)
            feature_output.add_value('Input Protein ID', input_id)
            feature_output.add_value('Primary Protein ID', prot_id)
            feature_output.add_value('Uniprot-Ac', u_ac)
            feature_output.add_value('WT Amino Acid', wt_aa)
            feature_output.add_value('Position', aac[1:-1])
            feature_output.add_value('Mut Amino Acid', new_aa)
            feature_output.add_value('AA change', '%s%s' % (wt_aa, new_aa))
            feature_output.add_value('Tags', tags)
            feature_output.add_value('Distance-based classification', Class)
            feature_output.add_value('Distance-based simple classification', simple_class)
            feature_output.add_value('RIN-based classification', rin_class)
            feature_output.add_value('RIN-based simple classification', rin_simple_class)
            feature_output.add_value('Classification confidence', conf)
            feature_output.add_value('Structure Location', weighted_sc)
            feature_output.add_value('Mainchain Location', mainchain_location)
            feature_output.add_value('Sidechain Location', sidechain_location)
            feature_output.add_value('RSA', rsa)
            feature_output.add_value('Mainchain RSA', mainchain_rsa)
            feature_output.add_value('Sidechain RSA', sidechain_rsa)
            feature_output.add_value('Amount of mapped structures', amount_of_structures)
            feature_output.add_value('Secondary structure assignment', mv_sec_ass)
            feature_output.add_value('IUPred value', iupred)
            feature_output.add_value('Region structure type', disorder_state)
            feature_output.add_value('Modres score', modres_score)
            feature_output.add_value('Phi', weighted_phi)
            feature_output.add_value('Psi', weighted_psi)

            KDmean = abs(sdsc.HYDROPATHY[wt_aa] - sdsc.HYDROPATHY[new_aa])
            feature_output.add_value('KD mean', KDmean)

            d_vol = abs(sdsc.VOLUME[wt_aa] - sdsc.VOLUME[new_aa])
            feature_output.add_value('Volume mean', d_vol)

            chemical_distance = database.getChemicalDistance(aac)
            feature_output.add_value('Chemical distance', chemical_distance)

            blosum_value = database.getBlosumValue(aac)
            feature_output.add_value('Blosum62', aac)

            aliphatic_change = int((wt_aa in sdsc.AA_MAP_ALIPHATIC) != (new_aa in sdsc.AA_MAP_ALIPHATIC))
            hydrophobic_change = int((wt_aa in sdsc.AA_MAP_HYDROPHOBIC) != (new_aa in sdsc.AA_MAP_HYDROPHOBIC))
            aromatic_change = int((wt_aa in sdsc.AA_MAP_AROMATIC) != (new_aa in sdsc.AA_MAP_AROMATIC))
            positive_change = int((wt_aa in sdsc.AA_MAP_POSITIVE) != (new_aa in sdsc.AA_MAP_POSITIVE))
            polar_change = int((wt_aa in sdsc.AA_MAP_POLAR) != (new_aa in sdsc.AA_MAP_POLAR))
            negative_change = int((wt_aa in sdsc.AA_MAP_NEGATIVE) != (new_aa in sdsc.AA_MAP_NEGATIVE))
            charged_change = int((wt_aa in sdsc.AA_MAP_CHARGED) != (new_aa in sdsc.AA_MAP_NEGATIVE))
            small_change = int((wt_aa in sdsc.AA_MAP_SMALL) != (new_aa in sdsc.AA_MAP_SMALL))
            tiny_change = int((wt_aa in sdsc.AA_MAP_TINY) != (new_aa in sdsc.AA_MAP_TINY))
            total_change = aliphatic_change + hydrophobic_change + aromatic_change + positive_change + polar_change + negative_change + charged_change + small_change + tiny_change
            feature_output.add_value('Aliphatic change', aliphatic_change)
            feature_output.add_value('Hydrophobic change', hydrophobic_change)
            feature_output.add_value('Aromatic change', aromatic_change)
            feature_output.add_value('Positive charged change', positive_change)
            feature_output.add_value('Polar change', polar_change)
            feature_output.add_value('Negative charge change', negative_change)
            feature_output.add_value('Charged change', charged_change)
            feature_output.add_value('Small change', small_change)
            feature_output.add_value('Tiny change', tiny_change)
            feature_output.add_value('Total change', total_change)
            feature_output.add_value('B Factor', b_factor)

            feature_output.add_value('AbsoluteCentrality', centrality_scores.AbsoluteCentrality)
            feature_output.add_value('LengthNormalizedCentrality', centrality_scores.LengthNormalizedCentrality)
            feature_output.add_value('MinMaxNormalizedCentrality', centrality_scores.MinMaxNormalizedCentrality)
            feature_output.add_value('AbsoluteCentralityWithNegative', centrality_scores.AbsoluteCentralityWithNegative)
            feature_output.add_value('LengthNormalizedCentralityWithNegative', centrality_scores.LengthNormalizedCentralityWithNegative)
            feature_output.add_value('MinMaxNormalizedCentralityWithNegative', centrality_scores.MinMaxNormalizedCentralityWithNegative)
            feature_output.add_value('AbsoluteComplexCentrality', centrality_scores.AbsoluteComplexCentrality)
            feature_output.add_value('LengthNormalizedComplexCentrality', centrality_scores.LengthNormalizedComplexCentrality)
            feature_output.add_value('MinMaxNormalizedComplexCentrality', centrality_scores.MinMaxNormalizedComplexCentrality)
            feature_output.add_value('AbsoluteComplexCentralityWithNegative', centrality_scores.AbsoluteComplexCentralityWithNegative)
            feature_output.add_value('LengthNormalizedComplexCentralityWithNegative', centrality_scores.LengthNormalizedComplexCentralityWithNegative)
            feature_output.add_value('MinMaxNormalizedComplexCentralityWithNegative', centrality_scores.MinMaxNormalizedComplexCentralityWithNegative)

            feature_output.add_value('Intra_SSBOND_Propensity', intra_ssbond_prop)
            feature_output.add_value('Inter_SSBOND_Propensity', inter_ssbond_prop)
            feature_output.add_value('Intra_Link_Propensity', intra_link_prop)
            feature_output.add_value('Inter_Link_Propensity', inter_link_prop)
            feature_output.add_value('CIS_Conformation_Propensity', cis_conformation_prop)
            feature_output.add_value('CIS_Follower_Propensity', cis_follower_prop)
            feature_output.add_value('Inter Chain Median KD', inter_chain_median_kd)
            feature_output.add_value('Inter Chain Distance Weighted KD', inter_chain_dist_weighted_kd)
            feature_output.add_value('Inter Chain Median RSA', inter_chain_median_rsa)
            feature_output.add_value('Inter Chain Distance Weighted RSA', inter_chain_dist_weighted_rsa)
            feature_output.add_value('Intra Chain Median KD', intra_chain_median_kd)
            feature_output.add_value('Intra Chain Distance Weighted KD', intra_chain_dist_weighted_kd)
            feature_output.add_value('Intra Chain Median RSA', intra_chain_median_rsa)
            feature_output.add_value('Intra Chain Distance Weighted RSA', intra_chain_dist_weighted_rsa)

            feature_output.add_value('Inter Chain Interactions Median', weighted_inter_chain_interactions_median)
            feature_output.add_value('Inter Chain Interactions Distance Weighted', weighted_inter_chain_interactions_dist_weighted)
            feature_output.add_value('Intra Chain Interactions Median', weighted_intra_chain_interactions_median)
            feature_output.add_value('Intra Chain Interactions Distance Weighted', weighted_intra_chain_interactions_dist_weighted)

            for chaintype in ['mc', 'sc']:
                for interaction_type in ['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide']:
                    feature_name = '%s %s score' % (chaintype, interaction_type)
                    value = rin_profile.getChainSpecificCombiScore(chaintype, interaction_type)
                    feature_output.add_value(feature_name, value)

                    feature_name = '%s %s degree' % (chaintype, interaction_type)
                    value = rin_profile.getChainSpecificCombiDegree(chaintype, interaction_type)
                    feature_output.add_value(feature_name, value)

                    feature_name = '%s %s H-bond score' % (chaintype, interaction_type)
                    value = rin_profile.getScore(chaintype, 'hbond', interaction_type)
                    feature_output.add_value(feature_name, value)
            f_lines.append(feature_output.pop_line())

            if len(f_lines) > max_lines:
                with FileLock(os.path.abspath('%s.lock' % feature_file_path)):
                    with open(feature_file_path, 'a') as feat_f:
                        feat_f.write(''.join(f_lines))
                f_lines = []

        if len(c_lines) > max_lines:
            with FileLock(os.path.abspath('%s.lock' % classification_file_path)):
                with open(classification_file_path, 'a') as f:
                    f.write(''.join(c_lines))
            c_lines = []

    with FileLock(os.path.abspath('%s.lock' % feature_file_path)):
        with open(feature_file_path, 'a') as feat_f:
            feat_f.write(''.join(f_lines))

    #print('Writing classification lines:', len(c_lines), classification_file_path)

    with FileLock(os.path.abspath('%s.lock' % classification_file_path)):
        with open(classification_file_path, 'a') as f:
            f.write(''.join(c_lines))
    return stat_dict

def classificationOutput(config, outfolder, session_name, session_id, ligand_filter=None):
    outfile = '%s/%s' % (outfolder, session_name)
    if config.verbosity >= 2:
        t0 = time.time()

    if ligand_filter is not None:
        f = open(ligand_filter, 'r')
        lines = f.readlines()
        f.close()
        ligand_filter = set()
        for line in lines:
            ligand_filter.add(line.replace('\n', '').replace(' ', ''))

    if config.verbosity >= 2:
        t1 = time.time()
        print("Time for classificationOutput part 1: ", t1 - t0)

    table = 'RS_Position_Session'
    rows = ['Position', 'Tag']
    eq_rows = {'Session': session_id}

    results = database.select(config, rows, table, equals_rows=eq_rows)

    if config.verbosity >= 2:
        t2 = time.time()
        print("Time for classificationOutput part 2: ", t2 - t1)

    tag_map = {}
    for row in results:
        mut_id = row[0]
        tags = row[1]
        tag_map[mut_id] = tags

    mutation_dict = database.getMutationDict(set(tag_map.keys()), config)

    table = 'RS_SNV_Session'
    columns = ['SNV', 'Tag']

    results = database.select(config, columns, table, equals_rows=eq_rows)

    snv_tag_map = {}
    snv_ids = []
    for row in results:
        snv_ids.append(row[0])
        snv_tag_map[row[0]] = row[1]

    table = 'SNV'
    columns = ['SNV_Id', 'Position', 'New_AA']

    results = database.binningSelect(snv_ids, columns, table, config)

    snv_map = {}
    for row in results:
        if not row[1] in snv_map:
            snv_map[row[1]] = {}
        snv_map[row[1]][row[0]] = row[2]

    prot_id_list = set()
    input_res_id_dict = {}
    for m in mutation_dict:
        prot_id_list.add(mutation_dict[m][1])
        input_res_id_dict[m] = mutation_dict[m][4]

    protein_dict = database.getProteinDict(prot_id_list, session_id, config)

    class_files = []

    class_file = "%s.classification.tsv" % (outfile)
    class_files.append(class_file)
    if os.path.isfile(class_file):
        os.remove(class_file)

    feature_file = "%s.features.tsv" % (outfile)
    if os.path.isfile(feature_file):
        os.remove(feature_file)

    stat_file = "%s.statistics.tsv" % (outfile)
    if os.path.isfile(stat_file):
        os.remove(stat_file)

    if config.verbosity >= 2:
        t3 = time.time()
        print("Time for classificationOutput part 3: ", t3 - t2)

    columns = ['Position_Id', 'Position_Number', 'Wildtype_Residue', 'Protein', 'IUPRED', 'IUPRED_Glob', 'Recommended_Structure_Data', 'Position_Data']

    table = 'Position'
    all_results = database.binningSelect(mutation_dict.keys(), columns, table, config)

    if config.verbosity >= 2:
        t4 = time.time()
        print("Time for classificationOutput part 4: ", t4 - t3)
        print('Total number of positions:', len(all_results))
        t5 = time.time()

    classification_output, f = init_classification_table(class_file)
    feature_output, feat_f = init_feature_table(feature_file)

    max_rows_at_a_time = 5000000

    use_ray = len(all_results) > config.proc_n

    if use_ray:
        #print('Store sizes:', sdsc.total_size(tag_map), sdsc.total_size(mutation_dict), sdsc.total_size(protein_dict))
        ray_init(config)
        data_store = ray.put((config, input_res_id_dict, class_file, feature_file))

    already_unpacked = False
    main_loop_counter = 0
    stat_dict = {}
    while((main_loop_counter)*max_rows_at_a_time <= len(all_results)):
        if config.verbosity >= 2:
                t4 = time.time()
        results = all_results[(main_loop_counter*max_rows_at_a_time):((main_loop_counter+1)*max_rows_at_a_time)]

        if config.verbosity >= 3:
            print('Processing row', main_loop_counter*max_rows_at_a_time, 'till', (main_loop_counter+1)*max_rows_at_a_time)

        main_loop_counter += 1

        if use_ray:
            f.close()
            feat_f.close()

            small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks = structman.utils.calculate_chunksizes(config.proc_n, len(results))
            if config.verbosity >= 3:
                print('calculate_chunksizes output:', small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks)

            unpackaging_remote_ids = []
            i = 0
            for i in range(n_of_big_chunks):

                left = i*big_chunksize
                right = (i+1)*big_chunksize

                sub_tag_map = {}
                sub_snv_map = {}
                sub_snv_tag_map = {}
                sub_protein_dict = {}
                for row in results[left:right]:
                    m = row[0]
                    position_number = row[1]
                    wt_aa = row[2]
                    prot_db_id = row[3]

                    sub_tag_map[m] = tag_map[m]

                    if m in snv_map:
                        sub_snv_map[m] = snv_map[m]

                        for snv_database_id in sub_snv_map[m]:
                            sub_snv_tag_map[snv_database_id] = snv_tag_map[snv_database_id]
                    if prot_db_id not in sub_protein_dict:
                        sub_protein_dict[prot_db_id] = protein_dict[prot_db_id]


                unpackaging_remote_ids.append(unpack_rows.remote(results[left:right], data_store, sub_tag_map, sub_snv_map, sub_snv_tag_map, sub_protein_dict))

                if config.verbosity >= 4:
                    print('Started unpack_rows, package number:', i)

            for j in range(n_of_small_chunks):

                left = ((i+1)*big_chunksize)+j*small_chunksize
                right = ((i+1)*big_chunksize)+(j+1)*small_chunksize

                sub_tag_map = {}
                sub_snv_map = {}
                sub_snv_tag_map = {}
                sub_protein_dict = {}
                for row in results[left:right]:
                    m = row[0]
                    position_number = row[1]
                    wt_aa = row[2]
                    prot_db_id = row[3]

                    sub_tag_map[m] = tag_map[m]

                    if m in snv_map:
                        sub_snv_map[m] = snv_map[m]

                        for snv_database_id in sub_snv_map[m]:
                            sub_snv_tag_map[snv_database_id] = snv_tag_map[snv_database_id]
                    if prot_db_id not in sub_protein_dict:
                        sub_protein_dict[prot_db_id] = protein_dict[prot_db_id]

                unpackaging_remote_ids.append(unpack_rows.remote(results[left:right], data_store, sub_tag_map, sub_snv_map, sub_snv_tag_map, sub_protein_dict))
                if config.verbosity >= 4:
                    print('Started unpack_rows, package number:', i+j)

            results = []
            while True:
                ready, not_ready = ray.wait(unpackaging_remote_ids)
                if len(ready) > 0:
                    for sub_stat_dict in ray.get(ready):
                        stat_dict.update(sub_stat_dict)
                unpackaging_remote_ids = not_ready
                if len(unpackaging_remote_ids) == 0:
                    break

            already_unpacked = True
        if config.verbosity >= 2:
            t5 = time.time()
            print("Time for classificationOutput part 5: ", t5 - t4, main_loop_counter)

        if not already_unpacked:
            from structman.utils import unpack
            for row in results:
                m = row[0]
                position_number = row[1]
                wt_aa = row[2]
                prot_db_id = row[3]
                iupred = row[4]
                disorder_state = row[5]

                if row[6] is not None:
                    (recommended_structure_str, max_seq_structure_str) = unpack(row[6])
                else:
                    recommended_structure_str = None
                    max_seq_structure_str = None
                row = list(row)
                if row[7] is None:
                    row[7] = [None]*38
                else:
                    try:
                        row[7] = unpack(row[7])
                    except:
                        print('Error in retrieving position data, couldnt unpack:', m)
                        continue

                (weighted_sc, mainchain_location, sidechain_location, rsa, mainchain_rsa, sidechain_rsa, Class, rin_class, simple_class, rin_simple_class,
                    interaction_str, conf, mv_sec_ass, amount_of_structures, rin_profile_str, modres_score, b_factor, weighted_centrality_scores, weighted_phi,
                    weighted_psi, intra_ssbond_prop, inter_ssbond_prop, intra_link_prop, inter_link_prop, cis_conformation_prop, cis_follower_prop,
                    inter_chain_median_kd, inter_chain_dist_weighted_kd, inter_chain_median_rsa, inter_chain_dist_weighted_rsa,
                    intra_chain_median_kd, intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                    weighted_inter_chain_interactions_median, weighted_inter_chain_interactions_dist_weighted,
                    weighted_intra_chain_interactions_median, weighted_intra_chain_interactions_dist_weighted) = row[7]

                centrality_scores = rin.Centrality_scores(code_str=weighted_centrality_scores)
                rin_profile = rin.Interaction_profile(profile_str=rin_profile_str)

                recommended_structure, seq_id, cov, resolution = sdsc.process_recommend_structure_str(recommended_structure_str)
                max_seq_structure, max_seq_seq_id, max_seq_cov, max_seq_resolution = sdsc.process_recommend_structure_str(max_seq_structure_str)

                input_res_id = mutation_dict[m][4]
                if input_res_id is None:
                    input_res_id = ''

                if disorder_state is not None:
                    if Class is None and (disorder_state == 'disorder' or disorder_state[0] == 'D'):
                        Class = 'Disorder'
                        simple_class = 'Disorder'
                        rin_class = 'Disorder'
                        rin_simple_class = 'Disorder'

                if max_seq_seq_id == '-':
                    stat_dict[m] = (Class, 0.)
                else:
                    stat_dict[m] = (Class, float(max_seq_seq_id))

                (prot_id, u_ac, refseq, u_id, error_code, error, input_id) = protein_dict[prot_db_id]

                if config.skipref:
                    refseq = ''
                classification_output.add_value('Input Protein ID', input_id)
                classification_output.add_value('Primary Protein ID', prot_id)
                classification_output.add_value('Uniprot-Ac', u_ac)
                classification_output.add_value('Uniprot ID', u_id)
                classification_output.add_value('Refseq', refseq)
                classification_output.add_value('Residue ID', input_res_id)
                classification_output.add_value('Amino Acid', wt_aa)
                classification_output.add_value('Position', position_number)
                classification_output.add_value('Tags', tag_map[m])
                classification_output.add_value('Weighted Location', weighted_sc)
                classification_output.add_value('Weighted Mainchain Location', mainchain_location)
                classification_output.add_value('Weighted Sidechain Location', sidechain_location)
                classification_output.add_value('Class', Class)
                classification_output.add_value('Simple Class', simple_class)
                classification_output.add_value('RIN Class', rin_class)
                classification_output.add_value('RIN Simple Class', rin_simple_class)
                classification_output.add_value('Individual Interactions', interaction_str)
                classification_output.add_value('Confidence Value', conf)
                classification_output.add_value('Secondary Structure', mv_sec_ass)
                classification_output.add_value('Recommended Structure', recommended_structure)
                classification_output.add_value('Sequence-ID', seq_id)
                classification_output.add_value('Coverage', cov)
                classification_output.add_value('Resolution', resolution)
                classification_output.add_value('Max Seq Id Structure', max_seq_structure)
                classification_output.add_value('Max Sequence-ID', max_seq_seq_id)
                classification_output.add_value('Max Seq Id Coverage', max_seq_cov)
                classification_output.add_value('Max Seq Id Resolution', max_seq_resolution)
                classification_output.add_value('Amount of mapped structures', amount_of_structures)

                f.write(classification_output.pop_line())

                if m not in snv_map:
                    snv_map[m] = {0: wt_aa}
                    snv_tag_map[0] = tag_map[m]

                for snv_database_id in snv_map[m]:
                    new_aa = snv_map[m][snv_database_id]
                    tags = snv_tag_map[snv_database_id]
                    aac = "%s%s%s" % (wt_aa, str(position_number), new_aa)
                    feature_output.add_value('Input Protein ID', input_id)
                    feature_output.add_value('Primary Protein ID', prot_id)
                    feature_output.add_value('Uniprot-Ac', u_ac)
                    feature_output.add_value('WT Amino Acid', wt_aa)
                    feature_output.add_value('Position', aac[1:-1])
                    feature_output.add_value('Mut Amino Acid', new_aa)
                    feature_output.add_value('AA change', '%s%s' % (wt_aa, new_aa))
                    feature_output.add_value('Tags', tags)
                    feature_output.add_value('Distance-based classification', Class)
                    feature_output.add_value('Distance-based simple classification', simple_class)
                    feature_output.add_value('RIN-based classification', rin_class)
                    feature_output.add_value('RIN-based simple classification', rin_simple_class)
                    feature_output.add_value('Classification confidence', conf)
                    feature_output.add_value('Structure Location', weighted_sc)
                    feature_output.add_value('Mainchain Location', mainchain_location)
                    feature_output.add_value('Sidechain Location', sidechain_location)
                    feature_output.add_value('RSA', rsa)
                    feature_output.add_value('Mainchain RSA', mainchain_rsa)
                    feature_output.add_value('Sidechain RSA', sidechain_rsa)
                    feature_output.add_value('Amount of mapped structures', amount_of_structures)
                    feature_output.add_value('Secondary structure assignment', mv_sec_ass)
                    feature_output.add_value('IUPred value', iupred)
                    feature_output.add_value('Region structure type', disorder_state)
                    feature_output.add_value('Modres score', modres_score)
                    feature_output.add_value('Phi', weighted_phi)
                    feature_output.add_value('Psi', weighted_psi)

                    KDmean = abs(sdsc.HYDROPATHY[wt_aa] - sdsc.HYDROPATHY[new_aa])
                    feature_output.add_value('KD mean', KDmean)

                    d_vol = abs(sdsc.VOLUME[wt_aa] - sdsc.VOLUME[new_aa])
                    feature_output.add_value('Volume mean', d_vol)

                    chemical_distance = database.getChemicalDistance(aac)
                    feature_output.add_value('Chemical distance', chemical_distance)

                    blosum_value = database.getBlosumValue(aac)
                    feature_output.add_value('Blosum62', aac)

                    aliphatic_change = int((wt_aa in sdsc.AA_MAP_ALIPHATIC) != (new_aa in sdsc.AA_MAP_ALIPHATIC))
                    hydrophobic_change = int((wt_aa in sdsc.AA_MAP_HYDROPHOBIC) != (new_aa in sdsc.AA_MAP_HYDROPHOBIC))
                    aromatic_change = int((wt_aa in sdsc.AA_MAP_AROMATIC) != (new_aa in sdsc.AA_MAP_AROMATIC))
                    positive_change = int((wt_aa in sdsc.AA_MAP_POSITIVE) != (new_aa in sdsc.AA_MAP_POSITIVE))
                    polar_change = int((wt_aa in sdsc.AA_MAP_POLAR) != (new_aa in sdsc.AA_MAP_POLAR))
                    negative_change = int((wt_aa in sdsc.AA_MAP_NEGATIVE) != (new_aa in sdsc.AA_MAP_NEGATIVE))
                    charged_change = int((wt_aa in sdsc.AA_MAP_CHARGED) != (new_aa in sdsc.AA_MAP_NEGATIVE))
                    small_change = int((wt_aa in sdsc.AA_MAP_SMALL) != (new_aa in sdsc.AA_MAP_SMALL))
                    tiny_change = int((wt_aa in sdsc.AA_MAP_TINY) != (new_aa in sdsc.AA_MAP_TINY))
                    total_change = aliphatic_change + hydrophobic_change + aromatic_change + positive_change + polar_change + negative_change + charged_change + small_change + tiny_change
                    feature_output.add_value('Aliphatic change', aliphatic_change)
                    feature_output.add_value('Hydrophobic change', hydrophobic_change)
                    feature_output.add_value('Aromatic change', aromatic_change)
                    feature_output.add_value('Positive charged change', positive_change)
                    feature_output.add_value('Polar change', polar_change)
                    feature_output.add_value('Negative charge change', negative_change)
                    feature_output.add_value('Charged change', charged_change)
                    feature_output.add_value('Small change', small_change)
                    feature_output.add_value('Tiny change', tiny_change)
                    feature_output.add_value('Total change', total_change)
                    feature_output.add_value('B Factor', b_factor)

                    feature_output.add_value('AbsoluteCentrality', centrality_scores.AbsoluteCentrality)
                    feature_output.add_value('LengthNormalizedCentrality', centrality_scores.LengthNormalizedCentrality)
                    feature_output.add_value('MinMaxNormalizedCentrality', centrality_scores.MinMaxNormalizedCentrality)
                    feature_output.add_value('AbsoluteCentralityWithNegative', centrality_scores.AbsoluteCentralityWithNegative)
                    feature_output.add_value('LengthNormalizedCentralityWithNegative', centrality_scores.LengthNormalizedCentralityWithNegative)
                    feature_output.add_value('MinMaxNormalizedCentralityWithNegative', centrality_scores.MinMaxNormalizedCentralityWithNegative)
                    feature_output.add_value('AbsoluteComplexCentrality', centrality_scores.AbsoluteComplexCentrality)
                    feature_output.add_value('LengthNormalizedComplexCentrality', centrality_scores.LengthNormalizedComplexCentrality)
                    feature_output.add_value('MinMaxNormalizedComplexCentrality', centrality_scores.MinMaxNormalizedComplexCentrality)
                    feature_output.add_value('AbsoluteComplexCentralityWithNegative', centrality_scores.AbsoluteComplexCentralityWithNegative)
                    feature_output.add_value('LengthNormalizedComplexCentralityWithNegative', centrality_scores.LengthNormalizedComplexCentralityWithNegative)
                    feature_output.add_value('MinMaxNormalizedComplexCentralityWithNegative', centrality_scores.MinMaxNormalizedComplexCentralityWithNegative)

                    feature_output.add_value('Intra_SSBOND_Propensity', intra_ssbond_prop)
                    feature_output.add_value('Inter_SSBOND_Propensity', inter_ssbond_prop)
                    feature_output.add_value('Intra_Link_Propensity', intra_link_prop)
                    feature_output.add_value('Inter_Link_Propensity', inter_link_prop)
                    feature_output.add_value('CIS_Conformation_Propensity', cis_conformation_prop)
                    feature_output.add_value('CIS_Follower_Propensity', cis_follower_prop)
                    feature_output.add_value('Inter Chain Median KD', inter_chain_median_kd)
                    feature_output.add_value('Inter Chain Distance Weighted KD', inter_chain_dist_weighted_kd)
                    feature_output.add_value('Inter Chain Median RSA', inter_chain_median_rsa)
                    feature_output.add_value('Inter Chain Distance Weighted RSA', inter_chain_dist_weighted_rsa)
                    feature_output.add_value('Intra Chain Median KD', intra_chain_median_kd)
                    feature_output.add_value('Intra Chain Distance Weighted KD', intra_chain_dist_weighted_kd)
                    feature_output.add_value('Intra Chain Median RSA', intra_chain_median_rsa)
                    feature_output.add_value('Intra Chain Distance Weighted RSA', intra_chain_dist_weighted_rsa)

                    feature_output.add_value('Inter Chain Interactions Median', weighted_inter_chain_interactions_median)
                    feature_output.add_value('Inter Chain Interactions Distance Weighted', weighted_inter_chain_interactions_dist_weighted)
                    feature_output.add_value('Intra Chain Interactions Median', weighted_intra_chain_interactions_median)
                    feature_output.add_value('Intra Chain Interactions Distance Weighted', weighted_intra_chain_interactions_dist_weighted)

                    for chaintype in ['mc', 'sc']:
                        for interaction_type in ['neighbor', 'short', 'long', 'ligand', 'ion', 'metal', 'Protein', 'DNA', 'RNA', 'Peptide']:
                            feature_name = '%s %s score' % (chaintype, interaction_type)
                            value = rin_profile.getChainSpecificCombiScore(chaintype, interaction_type)
                            feature_output.add_value(feature_name, value)

                            feature_name = '%s %s degree' % (chaintype, interaction_type)
                            value = rin_profile.getChainSpecificCombiDegree(chaintype, interaction_type)
                            feature_output.add_value(feature_name, value)

                            feature_name = '%s %s H-bond score' % (chaintype, interaction_type)
                            value = rin_profile.getScore(chaintype, 'hbond', interaction_type)
                            feature_output.add_value(feature_name, value)

                    feat_f.write(feature_output.pop_line())

            f.close()
            feat_f.close()

    if os.path.isfile('%s.lock' % feature_file):
        os.remove('%s.lock' % feature_file)

    if os.path.isfile('%s.lock' % class_file):
        os.remove('%s.lock' % class_file)

    if config.verbosity >= 2:
        t6 = time.time()
        print("Time for classificationOutput part 6: ", t6 - t5)

    writeStatFile(stat_file, mutation_dict, {}, tag_map, stat_dict=stat_dict)

    if config.verbosity >= 2:
        t7 = time.time()
        print("Time for classificationOutput part 7: ", t7 - t6)

    return class_files, []
