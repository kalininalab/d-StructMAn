import os
import time
from ast import literal_eval

from structman.lib.database import database
from structman.lib.output import out_generator

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    import warnings
    warnings.warn('Plotting not available, matplotlib not installed')


def makeViolins(violins, outfile, session_name, add=''):
    fs = 10  # fontsize
    plt.clf()

    for violin_tag in violins:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 12))
        classes = list(violins[violin_tag].keys())
        pos = list(range(len(classes)))
        data = [violins[violin_tag][cl] for cl in classes]
        # print classes
        axes.violinplot(data, pos, widths=0.5,
                        showmeans=True, showextrema=True, showmedians=True,
                        bw_method='silverman')

        axes.set_xticks(pos)
        axes.set_xticklabels(classes, rotation='vertical')
        # axes[n].set_xlabel(classes, rotation='vertical')

    # for ax in axes.flatten():
    #    ax.set_yticklabels([])

        fig.suptitle("Violin Plots %s - %s value" % (session_name, violin_tag))
        # fig.subplots_adjust(hspace=0.0)
        plt.tight_layout()
        # plt.show()
        violin_file = '%s.violin_plots_%s%s.png' % (outfile, violin_tag, add)
        plt.savefig(violin_file)


def plotScatter(scatter_plots, outfile, session_name):
    """
    x_values = []
    y_values = []

    for pdb in dssp_map:
        for (chain,res) in dssp_map[pdb]:
            (acc,aa) = dssp_map[pdb][(chain,res)]
            measure = measure_map[pdb][chain][res]
            x_values.append(acc)
            if not anti_corr:
                y_values.append(measure)
            else:
                y_values.append(-measure)

    fig, ax = plt.subplots()
    ax.scatter(x_values,y_values,alpha=0.01)
    ax.set_xlabel('RSA', fontsize=15)
    ax.set_ylabel(y_axis_name, fontsize=15)

    plt.savefig(outfile,dpi=300)
    """
    fs = 10  # fontsize
    plt.clf()
    scores = ['LI score', 'CI score', 'SI score', 'MI score', 'LoI score']
    degrees = ['LI degree', 'CI degree', 'SI degree', 'MI degree', 'LoI degree']
    labels = ['Ligand Interaction Score', 'Chain Interaction Score', 'Short Interaction Score', 'Medium Interaction Score', 'Long Interaction Score']
    for scatter_tag in scatter_plots:
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 12))

        n = 0
        for score in scores:
            b = n
            a = 0
            if n > 1:
                b = n - 2
                a = 1
            for Class in scatter_plots[scatter_tag][score]:
                if Class == 'Core':
                    axes[a, b].scatter(scatter_plots[scatter_tag][score][Class][0], scatter_plots[scatter_tag][score][Class][1], alpha=0.05, color='red')
                    axes[a, b].set_xlabel(labels[n], fontsize=fs)
                    axes[a, b].set_ylabel(scatter_tag, fontsize=fs)
                """else:
                    axes[a,b].scatter(scatter_plots[scatter_tag][score][Class][0],scatter_plots[scatter_tag][score][Class][1],alpha = 0.05,color='blue')
                    axes[a,b].set_xlabel(labels[n], fontsize=fs)
                    axes[a,b].set_ylabel(scatter_tag, fontsize=fs)"""
            n += 1

        fig.suptitle("Scatter Plots Interaction Scores %s - %s value" % (session_name, scatter_tag))
        # fig.subplots_adjust(hspace=0.0)
        plt.tight_layout()
        # plt.show()
        scatter_file = '%s.scatter_plots_Iscore_%s.png' % (outfile, scatter_tag)
        # print scatter_file
        plt.savefig(scatter_file)
        # print scatter_file
        plt.clf()
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 12))
        n = 0
        for score in degrees:

            b = n
            a = 0
            if n > 1:
                b = n - 2
                a = 1
            for Class in scatter_plots[scatter_tag][score]:
                if Class == 'Core':
                    axes[a, b].scatter(scatter_plots[scatter_tag][score][Class][0], scatter_plots[scatter_tag][score][Class][1], alpha=0.05, color='red')
                    axes[a, b].set_xlabel(labels[n], fontsize=fs)
                    axes[a, b].set_ylabel(scatter_tag, fontsize=fs)
                """else:
                    axes[a,b].scatter(scatter_plots[scatter_tag][score][Class][0],scatter_plots[scatter_tag][score][Class][1],alpha = 0.05,color='blue')
                    axes[a,b].set_xlabel(labels[n], fontsize=fs)
                    axes[a,b].set_ylabel(scatter_tag, fontsize=fs)"""
            n += 1

        fig.suptitle("Scatter Plots Interaction Degrees %s - %s value" % (session_name, scatter_tag))
        # fig.subplots_adjust(hspace=0.0)
        plt.tight_layout()
        # plt.show()
        scatter_file = '%s.scatter_plots_Idegrees_%s.png' % (outfile, scatter_tag)
        # print scatter_file
        plt.savefig(scatter_file)
        # print scatter_file


def create_ppi_network(session, config, outfile):

    verbose = config.verbose
    proteins = database.proteinsFromDb(session, config)

    # Step 1: build a residue -> position and a position -> residue map
    r_m_map = {}

    u_acs = proteins.get_protein_ids()

    if verbose:
        print('Amount of proteins:', len(u_acs))

    n_struct_anno = 0
    n_positions = 0
    n_sub_infos = 0

    for u_ac in u_acs:

        positions = proteins.get_position_ids(u_ac)

        annotation_list = proteins.get_protein_annotation_list(u_ac)

        n_struct_anno += len(annotation_list)
        n_positions += len(positions)

        for (pdb_id, chain) in annotation_list:

            sub_infos = proteins.get_sub_infos(u_ac, pdb_id, chain)
            n_sub_infos += len(sub_infos)
            for pos in positions:
                if pos not in sub_infos:
                    continue
                sub_info = sub_infos[pos]
                res_nr = sub_info[0]

                if not proteins.contains_residue(pdb_id, chain, res_nr):
                    continue
                r_id = proteins.get_residue_db_id(pdb_id, chain, res_nr)
                res_aa = proteins.get_residue_aa(pdb_id, chain, res_nr)
                if r_id not in r_m_map:
                    r_m_map[r_id] = []
                r_m_map[r_id].append((u_ac, pos))

    if verbose:
        print('Size of r_m_map:', len(r_m_map))
        print('Amount of structure annotations:', n_struct_anno)
        print('Amount of positions:', n_positions)
        print('Amount of sub_infos:', n_sub_infos)

    # Step 2: build a residue,chain A - residue,chain B datastructure
    struct_res_res_dict = {}
    for (pdb_id, chain) in proteins.structures:
        struct_res_res_dict[(pdb_id, chain)] = {}
        for res_nr in proteins.structures[(pdb_id, chain)].residues:
            scd = proteins.get_residue_scd(pdb_id, chain, res_nr)
            r_id = proteins.get_residue_db_id(pdb_id, chain, res_nr)

            scds = scd.split(",")

            for sd in scds:
                sdi = sd.split(":")
                if len(sdi) > 1:
                    ichain = sdi[0][0]
                    res_nr_2 = sdi[0].split('.')[1]
                    mc_d = float(sdi[1])
                    if mc_d < 5.0:
                        if not proteins.contains_structure(pdb_id, ichain):
                            continue
                        if not (r_id, res_nr) in struct_res_res_dict[(pdb_id, chain)]:
                            struct_res_res_dict[(pdb_id, chain)][(r_id, res_nr)] = set([])
                        struct_res_res_dict[(pdb_id, chain)][(r_id, res_nr)].add((ichain, res_nr_2))

    if verbose:
        print('Size of struct_res_res_dict:', len(struct_res_res_dict))

    if os.path.exists(outfile):
        os.remove(outfile)

    f = open(outfile, 'a')
    f.write('u_ac_1\taacbase_1\tpdb_id\tchain_1\tres_nr_1\tu_ac_2\taacbase_2\tchain_2\tres_nr_2\n')

    # Step 3 produce the table

    for pdb_id, chain_1 in struct_res_res_dict:
        for r_1, res_nr_1 in struct_res_res_dict[(pdb_id, chain_1)]:
            if r_1 not in r_m_map:
                continue
            for (chain_2, res_nr_2) in struct_res_res_dict[(pdb_id, chain_1)][(r_1, res_nr_1)]:

                r_2 = proteins.get_residue_db_id(pdb_id, chain_2, res_nr_2)

                if r_2 not in r_m_map:
                    continue

                for (u_ac_1, aacbase_1) in r_m_map[r_1]:
                    for (u_ac_2, aacbase_2) in r_m_map[r_2]:

                        f.write('\t'.join([u_ac_1, aacbase_1, pdb_id, chain_1, res_nr_1, u_ac_2, aacbase_2, pdb_id, chain_2, res_nr_2]) + '\n')
    f.close()


def classDistributionFromFile(annotationfile, outfolder, session_name, config, by_conf=False, rin_classes=False):
    #"Uniprot-Ac\tUniprot Id\tRefseq\tPDB-ID (Input)\tResidue-Id\tAmino Acid\tPosition\tSpecies\tTag\tWeighted Surface/Core\tClass\tSimple Class\tConfidence Value\tSecondary Structure\tRecommended Structure\tSequence-ID\tCoverage\tResolution\tMax Seq Id Structure\tMax Sequence-ID\tMax Seq Id Coverage\tMax Seq Id Resolution\tAmount of mapped structures"
    outfile = '%s/%s' % (outfolder, session_name)

    t0 = time.time()

    f = open(annotationfile, 'r')
    lines = f.readlines()
    f.close()

    tag_map = {}
    simple_tag_map = {}
    tag_sizes = {}

    class_map = {}
    simple_class_map = {}
    size = 0

    simple_high_confidence_map = {}
    hc_size = 0
    hc_threshs = [0.0, 0.2, 0.4, 0.6, 0.8]

    violins = {}
    comp_violins = {}

    tag_pos = None
    class_pos = None
    simple_class_pos = None
    confidence_pos = None
    rin_class_pos = None
    rin_simple_class_pos = None

    for pos, column_name in enumerate(lines[0].split('\t')):
        if column_name == 'Tags':
            tag_pos = pos
        if column_name == 'Class':
            class_pos = pos
        if column_name == 'Simple Class':
            simple_class_pos = pos
        if column_name == 'Confidence Value':
            confidence_pos = pos
        if column_name == 'RIN Class':
            rin_class_pos = pos
        if column_name == 'RIN Simple Class':
            rin_simple_class_pos = pos

    for line in lines[1:]:
        words = line.replace('\n', '').split('\t')

        tag = words[tag_pos]
        if not rin_classes:
            classification = words[class_pos]
            simple_classification = words[simple_class_pos]
        else:
            classification = words[rin_class_pos]
            simple_classification = words[rin_simple_class_pos]
        if words[confidence_pos] != out_generator.OutputGenerator.null_symbol:
            confidence = float(words[confidence_pos])
        else:
            confidence = 0.

        size += 1

        if classification not in class_map:
            class_map[classification] = 1
        else:
            class_map[classification] += 1

        if simple_classification not in simple_class_map:
            simple_class_map[simple_classification] = 1
        else:
            simple_class_map[simple_classification] += 1

        try:
            tag_dict = literal_eval(tag)
        except:
            tag_dict = {}
        for aa2 in tag_dict:
            tags = tag_dict[aa2]
            for tag in tags.split(','):
                if tag == '':
                    continue
                if tag[0] == '#':
                    if tag.count(':') == 1:
                        violin_tag, violin_value = tag[1:].split(':')
                    else:
                        violin_tag, violin_value = tag[1:].split('=')
                    violin_value = float(violin_value)

                    if violin_tag not in violins:
                        violins[violin_tag] = {}

                    if simple_classification not in violins[violin_tag]:
                        violins[violin_tag][simple_classification] = []
                    violins[violin_tag][simple_classification].append(violin_value)

                    if violin_tag not in comp_violins:
                        comp_violins[violin_tag] = {}

                    if classification not in comp_violins[violin_tag]:
                        comp_violins[violin_tag][classification] = []
                    comp_violins[violin_tag][classification].append(violin_value)

                    continue

                if tag not in tag_map:
                    tag_map[tag] = {}
                    simple_tag_map[tag] = {}
                    tag_sizes[tag] = 0
                if classification not in tag_map[tag]:
                    tag_map[tag][classification] = 1
                else:
                    tag_map[tag][classification] += 1
                if simple_classification not in simple_tag_map[tag]:
                    simple_tag_map[tag][simple_classification] = 1
                else:
                    simple_tag_map[tag][simple_classification] += 1
                tag_sizes[tag] += 1

        if by_conf:  # TODO
            for hc_thresh in hc_threshs:
                if confidence > hc_thresh:
                    if simple_classification not in simple_high_confidence_map:
                        simple_high_confidence_map[simple_classification] = 1
                    else:
                        simple_high_confidence_map[simple_classification] += 1
                    hc_size += 1

    t1 = time.time()
    if config.verbosity >= 2:
        print(('Time for classDistribution Part1: %s' % str(t1 - t0)))

    makeViolins(violins, outfile, session_name)
    makeViolins(comp_violins, outfile, session_name, add='_complex_classes')

    t2 = time.time()
    if config.verbosity >= 2:
        print(('Time for classDistribution Part2: %s' % str(t2 - t1)))

    classes = list(class_map.keys())
    outlines = ['Tag\t%s' % '\t'.join(classes)]

    words = ['total']
    for classification in classes:
        r = float(class_map[classification]) / float(size)
        words.append(str(r))
    outlines.append('\t'.join(words))

    if len(tag_map) > 1:
        for tag in tag_map:
            words = [tag]
            for classification in classes:
                if classification not in tag_map[tag]:
                    r = 0.0
                else:
                    r = float(tag_map[tag][classification]) / float(tag_sizes[tag])
                words.append(str(r))
            outlines.append('\t'.join(words))

    simple_classes = list(simple_class_map.keys())
    simple_outlines = ['Tag\t%s' % '\t'.join(simple_classes)]

    words = ['total']
    for classification in simple_classes:
        r = float(simple_class_map[classification]) / float(size)
        words.append(str(r))
    simple_outlines.append('\t'.join(words))

    if len(simple_tag_map) > 1:
        for tag in simple_tag_map:
            words = [tag]
            for classification in simple_classes:
                if classification not in simple_tag_map[tag]:
                    r = 0.0
                else:
                    r = float(simple_tag_map[tag][classification]) / float(tag_sizes[tag])
                words.append(str(r))
            simple_outlines.append('\t'.join(words))

    if by_conf:
        simple_high_confidence_outlines = ['\t'.join(list(simple_high_confidence_map.keys()))]
        words = []
        for simple_classification in simple_high_confidence_map:
            r = float(simple_high_confidence_map[simple_classification]) / float(hc_size)
            words.append(str(r))
        simple_high_confidence_outlines.append('\t'.join(words))

    if rin_classes:
        file_name_tag = 'rin_'
    else:
        file_name_tag = ''

    f = open('%s.%sclass_distribution.tsv' % (outfile, file_name_tag), 'w')
    f.write('\n'.join(outlines))
    f.close()

    f = open('%s.%ssimple_class_distribution.tsv' % (outfile, file_name_tag), 'w')
    f.write('\n'.join(simple_outlines))
    f.close()

    if by_conf:
        f = open('%s.simple_high_confidence_class_distribution.tsv' % outfile, 'w')
        f.write('\n'.join(simple_high_confidence_outlines))
        f.close()

    t3 = time.time()
    if config.verbosity >= 2:
        print(('Time for c lassDistribution Part3: %s' % str(t3 - t2)))


def InteractionScoreAveragesFromFile(InteractionProfilesfile, outfile, session_name, by_tag=False):
    outfile = "%s/%s" % (outfile, session_name)
    f = open(InteractionProfilesfile, 'r')
    lines = f.readlines()
    f.close()

    min_degree = None
    max_degree = None
    min_score = None
    max_score = None

    if len(lines) == 1:
        return

    for line in lines[1:]:
        row = line.replace('\n', '').split('\t')
        Ligand_Interaction_Degree = float(row[4])
        Ligand_Interaction_Score = float(row[5])
        Chain_Interaction_Degree = float(row[6])
        Chain_Interaction_Score = float(row[7])
        Short_Interaction_Degree = float(row[8])
        Short_Interaction_Score = float(row[9])
        Medium_Interaction_Degree = float(row[10])
        Medium_Interaction_Score = float(row[11])
        Long_Interaction_Degree = float(row[12])
        Long_Interaction_Score = float(row[13])

        degrees = [Ligand_Interaction_Degree, Chain_Interaction_Degree, Short_Interaction_Degree, Medium_Interaction_Degree, Long_Interaction_Degree]
        scores = [Ligand_Interaction_Score, Chain_Interaction_Score, Short_Interaction_Score, Medium_Interaction_Score, Long_Interaction_Score]
        for degree in degrees:
            if min_degree is None or degree < min_degree:
                min_degree = degree
            if max_degree is None or degree > max_degree:
                max_degree = degree

        for score in scores:
            if min_score is None or score < min_score:
                min_score = score
            if max_score is None or score > max_score:
                max_score = score

    bins = 50
    degree_bin_size = (max_degree - min_degree) / float(bins)
    score_bin_size = (max_score - min_score) / float(bins)

    degree_tag_map = {}
    score_tag_map = {}
    degree_tag_sizes = {}
    score_tag_sizes = {}

    scatter_plots = {}

    for line in lines[1:]:
        row = line.replace('\n', '').split('\t')
        Ligand_Interaction_Degree = float(row[4])
        Ligand_Interaction_Score = float(row[5])
        Chain_Interaction_Degree = float(row[6])
        Chain_Interaction_Score = float(row[7])
        Short_Interaction_Degree = float(row[8])
        Short_Interaction_Score = float(row[9])
        Medium_Interaction_Degree = float(row[10])
        Medium_Interaction_Score = float(row[11])
        Long_Interaction_Degree = float(row[12])
        Long_Interaction_Score = float(row[13])
        Class = row[14]
        comp_class = row[15]

        degrees = {'LI degree': Ligand_Interaction_Degree, 'CI degree': Chain_Interaction_Degree, 'SI degree': Short_Interaction_Degree, 'MI degree': Medium_Interaction_Degree, 'LoI degree': Long_Interaction_Degree}
        scores = {'LI score': Ligand_Interaction_Score, 'CI score': Chain_Interaction_Score, 'SI score': Short_Interaction_Score, 'MI score': Medium_Interaction_Score, 'LoI score': Long_Interaction_Score}

        if by_tag:
            tags = row[3]
            tags = '%s,All' % tags
            for tag in tags.split(','):
                if tag[0] == '#':
                    scatter_tag, scatter_value = tag[1:].split(':')
                    scatter_value = float(scatter_value)

                    if scatter_tag not in scatter_plots:
                        scatter_plots[scatter_tag] = {}

                        for score in scores:
                            scatter_plots[scatter_tag][score] = {}

                        for degree in degrees:
                            scatter_plots[scatter_tag][degree] = {}

                    for score in scores:
                        if Class not in scatter_plots[scatter_tag][score]:
                            scatter_plots[scatter_tag][score][Class] = [[], []]
                    for score in scores:
                        sc = scores[score]
                        scatter_plots[scatter_tag][score][Class][0].append(sc)
                        scatter_plots[scatter_tag][score][Class][1].append(scatter_value)

                    for degree in degrees:
                        if Class not in scatter_plots[scatter_tag][degree]:
                            scatter_plots[scatter_tag][degree][Class] = [[], []]
                    for degree in degrees:
                        deg = degrees[degree]
                        scatter_plots[scatter_tag][degree][Class][0].append(deg)
                        scatter_plots[scatter_tag][degree][Class][1].append(scatter_value)

                    continue
                for degree_name in degrees:
                    combined_tag = "%s %s" % (tag, degree_name)
                    if combined_tag not in degree_tag_map:
                        degree_tag_map[combined_tag] = [0] * bins
                    degree = degrees[degree_name]
                    for i in range(1, bins + 1):
                        if degree < min_degree + i * degree_bin_size:
                            degree_tag_map[combined_tag][i - 1] += 1
                            break
                    if combined_tag not in degree_tag_sizes:
                        degree_tag_sizes[combined_tag] = 1
                    else:
                        degree_tag_sizes[combined_tag] += 1
                for score_name in scores:
                    combined_tag = "%s %s" % (tag, score_name)
                    if combined_tag not in score_tag_map:
                        score_tag_map[combined_tag] = [0] * bins
                    score = scores[score_name]
                    for i in range(1, bins + 1):
                        if score < min_score + i * score_bin_size:
                            score_tag_map[combined_tag][i - 1] += 1
                            break
                    if combined_tag not in score_tag_sizes:
                        score_tag_sizes[combined_tag] = 1
                    else:
                        score_tag_sizes[combined_tag] += 1
        else:
            pass  # TODO

    plotScatter(scatter_plots, outfile, session_name)

    # print tag_sizes
    startline_words = ['']
    for i in range(0, bins):
        startline_words.append('>%s' % str(min_degree + i * degree_bin_size))

    outlines = ['\t'.join(startline_words)]
    for tag in degree_tag_map:
        words = [tag]
        for amount in degree_tag_map[tag]:
            r = amount / float(degree_tag_sizes[tag])
            words.append(str(r))
        outlines.append('\t'.join(words))

    startline_words = ['']
    for i in range(0, bins):
        startline_words.append('>%s' % str(min_score + i * score_bin_size))

    outlines.append('\t'.join(startline_words))
    for tag in score_tag_map:
        words = [tag]
        for amount in score_tag_map[tag]:
            r = amount / float(score_tag_sizes[tag])
            words.append(str(r))
        outlines.append('\t'.join(words))

    f = open('%s.interaction_profile_means.tsv' % outfile, 'w')
    f.write('\n'.join(outlines))
    f.close()
