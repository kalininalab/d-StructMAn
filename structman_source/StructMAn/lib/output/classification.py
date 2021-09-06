import os
import time
from ast import literal_eval

from structman.lib import database, rin, sdsc
from structman.lib.output.feature import init_feature_table
from structman.lib.output.output import OutputGenerator
from structman.lib.output.utils import makeViolins


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


def init_classification_table(class_file):
    classification_output = OutputGenerator()
    headers = [
        'Input Protein ID', 'Primary Protein ID', 'Uniprot-Ac', 'Uniprot ID', 'Refseq', 'Residue ID', 'Amino Acid', 'Position', 'Tags',
        'Weighted Location', 'Weighted Mainchain Location', 'Weighted Sidechain Location',
        'Class', 'Simple Class', 'RIN Class', 'RIN Simple Class', 'Individual Interactions',
        'Confidence Value', 'Secondary Structure', 'Recommended Structure', 'Sequence-ID', 'Coverage', 'Resolution',
        'Max Seq Id Structure', 'Max Sequence-ID', 'Max Seq Id Coverage', 'Max Seq Id Resolution', 'Amount of mapped structures'
    ]
    classification_output.add_headers(headers)
    f = open(class_file, 'a')
    f.write(classification_output.get_header())

    return classification_output, f


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
    for m in mutation_dict:
        prot_id_list.add(mutation_dict[m][1])

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

    columns = ['Position_Id', 'Amino_Acid_Change', 'Protein', 'Location', 'Mainchain_Location', 'Sidechain_Location',
               'Weighted_Surface_Access', 'Weighted_Surface_Access_Main_Chain', 'Weighted_Surface_Access_Side_Chain',
               'Class', 'Simple_Class', 'RIN_Class', 'RIN_Simple_Class', 'Interactions', 'Confidence',
               'Secondary_Structure', 'Recommended_Structure', 'Max_Seq_Structure', 'Mapped_Structures', 'RIN_Profile', 'IUPRED', 'IUPRED_Glob',
               'Modres_Score',
               'B_Factor', 'Weighted_Centrality_Scores', 'Weighted_Phi', 'Weighted_Psi', 'Intra_SSBOND_Propensity',
               'Inter_SSBOND_Propensity', 'Intra_Link_Propensity', 'Inter_Link_Propensity', 'CIS_Conformation_Propensity', 'CIS_Follower_Propensity',
               'Weighted_Inter_Chain_Median_KD', 'Weighted_Inter_Chain_Dist_Weighted_KD', 'Weighted_Inter_Chain_Median_RSA',
               'Weighted_Inter_Chain_Dist_Weighted_RSA', 'Weighted_Intra_Chain_Median_KD', 'Weighted_Intra_Chain_Dist_Weighted_KD',
               'Weighted_Intra_Chain_Median_RSA', 'Weighted_Intra_Chain_Dist_Weighted_RSA',
               'Weighted_Inter_Chain_Interactions_Median', 'Weighted_Inter_Chain_Interactions_Dist_Weighted',
               'Weighted_Intra_Chain_Interactions_Median', 'Weighted_Intra_Chain_Interactions_Dist_Weighted'
               ]

    table = 'Position'
    results = database.binningSelect(mutation_dict.keys(), columns, table, config)

    if config.verbosity >= 2:
        t4 = time.time()
        print("Time for classificationOutput part 4: ", t4 - t3)

    classification_output, f = init_classification_table(class_file)
    feature_output, feat_f = init_feature_table(feature_file)

    stat_dict = {}
    for row in results:
        m = row[0]
        aac = row[1]
        prot_db_id = row[2]
        weighted_sc = row[3]
        mainchain_location = row[4]
        sidechain_location = row[5]
        rsa = row[6]
        mainchain_rsa = row[7]
        sidechain_rsa = row[8]
        Class = row[9]
        simple_class = row[10]
        rin_class = row[11]
        rin_simple_class = row[12]
        interaction_str = row[13]
        conf = row[14]
        mv_sec_ass = row[15]
        recommended_structure_str = row[16]
        max_seq_structure_str = row[17]
        amount_of_structures = row[18]
        rin_profile_str = row[19]
        iupred = row[20]
        disorder_state = row[21]
        modres_score = row[22]
        b_factor = row[23]
        weighted_centrality_scores = row[24]
        weighted_phi = row[25]
        weighted_psi = row[26]
        intra_ssbond_prop = row[27]
        inter_ssbond_prop = row[28]
        intra_link_prop = row[29]
        inter_link_prop = row[30]
        cis_conformation_prop = row[31]
        cis_follower_prop = row[32]
        inter_chain_median_kd = row[33]
        inter_chain_dist_weighted_kd = row[34]
        inter_chain_median_rsa = row[35]
        inter_chain_dist_weighted_rsa = row[36]
        intra_chain_median_kd = row[37]
        intra_chain_dist_weighted_kd = row[38]
        intra_chain_median_rsa = row[39]
        intra_chain_dist_weighted_rsa = row[40]
        weighted_inter_chain_interactions_median = row[41]
        weighted_inter_chain_interactions_dist_weighted = row[42]
        weighted_intra_chain_interactions_median = row[43]
        weighted_intra_chain_interactions_dist_weighted = row[44]

        input_res_id = mutation_dict[m][4]
        if input_res_id is None:
            input_res_id = ''

        recommended_structure, seq_id, cov, resolution = sdsc.process_recommend_structure_str(recommended_structure_str)

        if disorder_state is not None:
            if Class is None and (disorder_state == 'disorder' or disorder_state[0] == 'D'):
                Class = 'Disorder'
                simple_class = 'Disorder'
                rin_class = 'Disorder'
                rin_simple_class = 'Disorder'

        max_seq_structure, max_seq_seq_id, max_seq_cov, max_seq_resolution = sdsc.process_recommend_structure_str(max_seq_structure_str)
        if max_seq_seq_id == '-':
            stat_dict[m] = (Class, 0.)
        else:
            stat_dict[m] = (Class, float(max_seq_seq_id))

        (prot_id, u_ac, refseq, u_id, error_code, error, input_id) = protein_dict[prot_db_id]

        if config.skipref:
            refseq = ''
        aa1 = aac[0]
        classification_output.add_value('Input Protein ID', input_id)
        classification_output.add_value('Primary Protein ID', prot_id)
        classification_output.add_value('Uniprot-Ac', u_ac)
        classification_output.add_value('Uniprot ID', u_id)
        classification_output.add_value('Refseq', refseq)
        classification_output.add_value('Residue ID', input_res_id)
        classification_output.add_value('Amino Acid', aac[0])
        classification_output.add_value('Position', aac[1:])
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

        centrality_scores = rin.Centrality_scores(code_str=weighted_centrality_scores)

        rin_profile = rin.Interaction_profile(profile_str=rin_profile_str)
        if m not in snv_map:
            snv_map[m] = {0: aa1}
            snv_tag_map[0] = tag_map[m]
        for snv_database_id in snv_map[m]:
            new_aa = snv_map[m][snv_database_id]
            tags = snv_tag_map[snv_database_id]
            aac = "%s%s" % (aac, new_aa)
            feature_output.add_value('Input Protein ID', input_id)
            feature_output.add_value('Primary Protein ID', prot_id)
            feature_output.add_value('Uniprot-Ac', u_ac)
            feature_output.add_value('WT Amino Acid', aa1)
            feature_output.add_value('Position', aac[1:-1])
            feature_output.add_value('Mut Amino Acid', new_aa)
            feature_output.add_value('AA change', '%s%s' % (aa1, new_aa))
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

            KDmean = abs(sdsc.HYDROPATHY[aa1] - sdsc.HYDROPATHY[new_aa])
            feature_output.add_value('KD mean', KDmean)

            d_vol = abs(sdsc.VOLUME[aa1] - sdsc.VOLUME[new_aa])
            feature_output.add_value('Volume mean', d_vol)

            chemical_distance = database.getChemicalDistance(aac)
            feature_output.add_value('Chemical distance', chemical_distance)

            blosum_value = database.getBlosumValue(aac)
            feature_output.add_value('Blosum62', aac)

            aliphatic_change = int((aa1 in sdsc.AA_MAP_ALIPHATIC) != (new_aa in sdsc.AA_MAP_ALIPHATIC))
            hydrophobic_change = int((aa1 in sdsc.AA_MAP_HYDROPHOBIC) != (new_aa in sdsc.AA_MAP_HYDROPHOBIC))
            aromatic_change = int((aa1 in sdsc.AA_MAP_AROMATIC) != (new_aa in sdsc.AA_MAP_AROMATIC))
            positive_change = int((aa1 in sdsc.AA_MAP_POSITIVE) != (new_aa in sdsc.AA_MAP_POSITIVE))
            polar_change = int((aa1 in sdsc.AA_MAP_POLAR) != (new_aa in sdsc.AA_MAP_POLAR))
            negative_change = int((aa1 in sdsc.AA_MAP_NEGATIVE) != (new_aa in sdsc.AA_MAP_NEGATIVE))
            charged_change = int((aa1 in sdsc.AA_MAP_CHARGED) != (new_aa in sdsc.AA_MAP_NEGATIVE))
            small_change = int((aa1 in sdsc.AA_MAP_SMALL) != (new_aa in sdsc.AA_MAP_SMALL))
            tiny_change = int((aa1 in sdsc.AA_MAP_TINY) != (new_aa in sdsc.AA_MAP_TINY))
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

    if config.verbosity >= 2:
        t5 = time.time()
        print("Time for classificationOutput part 5: ", t5 - t4)

    writeStatFile(stat_file, mutation_dict, {}, tag_map, stat_dict=stat_dict)

    if config.verbosity >= 2:
        t6 = time.time()
        print("Time for classificationOutput part 6: ", t6 - t5)

    return class_files, []
