import pymysql as MySQLdb
import database
import sdsc
import rin
import os
import babel
import time
from ast import literal_eval
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except:
    pass
import numpy as np

def create_ppi_network(session,config,outfile):

    verbose = config.verbose
    proteins = database.proteinsFromDb(session,config)

    #Step 1: build a residue -> position and a position -> residue map
    r_m_map = {}

    u_acs = proteins.get_protein_u_acs()

    if verbose:
        print('Amount of proteins:',len(u_acs))

    n_struct_anno = 0
    n_positions = 0
    n_sub_infos = 0

    for u_ac in u_acs:

        positions = proteins.get_position_ids(u_ac)

        annotation_list = proteins.get_protein_annotation_list(u_ac)

        n_struct_anno += len(annotation_list)
        n_positions += len(positions)

        for (pdb_id,chain) in annotation_list:

            sub_infos = proteins.get_sub_infos(u_ac,pdb_id,chain)
            n_sub_infos += len(sub_infos)
            for pos in positions:
                aacbase = proteins.get_aac_base(u_ac,pos)
                if not aacbase in sub_infos:
                    continue
                sub_info = sub_infos[aacbase]
                res_nr = sub_info[0]

                if not proteins.contains_residue(pdb_id,chain,res_nr):
                    continue
                r_id = proteins.get_residue_db_id(pdb_id,chain,res_nr)
                res_aa = proteins.get_residue_aa(pdb_id,chain,res_nr)
                if not r_id in r_m_map:
                    r_m_map[r_id] = []
                r_m_map[r_id].append((u_ac,aacbase))

    if verbose:
        print('Size of r_m_map:',len(r_m_map))
        print('Amount of structure annotations:',n_struct_anno)
        print('Amount of positions:',n_positions)
        print('Amount of sub_infos:',n_sub_infos)

    #Step 2: build a residue,chain A - residue,chain B datastructure
    struct_res_res_dict = {}
    for (pdb_id,chain) in proteins.structures:
        struct_res_res_dict[(pdb_id,chain)] = {}
        for res_nr in proteins.structures[(pdb_id,chain)].residues:
            scd = proteins.get_residue_scd(pdb_id,chain,res_nr)
            r_id = proteins.get_residue_db_id(pdb_id,chain,res_nr)

            scds = scd.split(",")

            for sd in scds:
                sdi = sd.split(":")
                if len(sdi) > 1:
                    ichain = sdi[0][0]
                    res_nr_2 = sdi[0].split('.')[1]
                    mc_d = float(sdi[1])
                    if mc_d < 5.0:
                        if not proteins.contains_structure(pdb_id,ichain):
                            continue
                        if not (r_id,res_nr) in struct_res_res_dict[(pdb_id,chain)]:
                            struct_res_res_dict[(pdb_id,chain)][(r_id,res_nr)] = set([])
                        struct_res_res_dict[(pdb_id,chain)][(r_id,res_nr)].add((ichain,res_nr_2))

    if verbose:
        print('Size of struct_res_res_dict:',len(struct_res_res_dict))

    if os.path.exists(outfile):
        os.remove(outfile)

    f = open(outfile,'a')
    f.write('u_ac_1\taacbase_1\tpdb_id\tchain_1\tres_nr_1\tu_ac_2\taacbase_2\tchain_2\tres_nr_2\n')

    #Step 3 produce the table

    for pdb_id,chain_1 in struct_res_res_dict:
        for r_1,res_nr_1 in struct_res_res_dict[(pdb_id,chain_1)]:
            if not r_1 in r_m_map:
                continue
            for (chain_2,res_nr_2) in struct_res_res_dict[(pdb_id,chain_1)][(r_1,res_nr_1)]:

                r_2 = proteins.get_residue_db_id(pdb_id,chain_2,res_nr_2)
                
                if not r_2 in r_m_map:
                    continue
                
                for (u_ac_1,aacbase_1) in r_m_map[r_1]:
                    for (u_ac_2,aacbase_2) in r_m_map[r_2]:

                        f.write('\t'.join([u_ac_1,aacbase_1,pdb_id,chain_1,res_nr_1,u_ac_2,aacbase_2,pdb_id,chain_2,res_nr_2]) + '\n')
    f.close()

def makeViolins(violins,outfile,session_name,add=''):
    fs = 10  # fontsize
    plt.clf()

    for violin_tag in violins:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 12))
        classes = list(violins[violin_tag].keys())
        pos = list(range(len(classes)))
        data = [violins[violin_tag][cl] for cl in classes]
        #print classes
        axes.violinplot(data, pos, widths=0.5,
                              showmeans=True, showextrema=True, showmedians=True,
                              bw_method='silverman')

        axes.set_xticks(pos)
        axes.set_xticklabels(classes,rotation='vertical')
        #axes[n].set_xlabel(classes, rotation='vertical')


    #for ax in axes.flatten():
    #    ax.set_yticklabels([])

        fig.suptitle("Violin Plots %s - %s value" % (session_name,violin_tag))
        #fig.subplots_adjust(hspace=0.0)
        plt.tight_layout()
        #plt.show()
        violin_file = '%s.violin_plots_%s%s.png' % (outfile,violin_tag,add)
        plt.savefig(violin_file)


def plotScatter(scatter_plots,outfile,session_name):

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
    scores = ['LI score','CI score','SI score','MI score','LoI score']
    degrees = ['LI degree','CI degree','SI degree','MI degree','LoI degree']
    labels = ['Ligand Interaction Score','Chain Interaction Score','Short Interaction Score','Medium Interaction Score','Long Interaction Score']
    for scatter_tag in scatter_plots:
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(12, 12))

        n = 0
        for score in scores:
            b = n
            a = 0
            if n > 1:
                b = n -2
                a = 1
            for Class in scatter_plots[scatter_tag][score]:
                if Class == 'Core':
                    axes[a,b].scatter(scatter_plots[scatter_tag][score][Class][0],scatter_plots[scatter_tag][score][Class][1],alpha = 0.05,color='red')
                    axes[a,b].set_xlabel(labels[n], fontsize=fs)
                    axes[a,b].set_ylabel(scatter_tag, fontsize=fs)
                """else:
                    axes[a,b].scatter(scatter_plots[scatter_tag][score][Class][0],scatter_plots[scatter_tag][score][Class][1],alpha = 0.05,color='blue')
                    axes[a,b].set_xlabel(labels[n], fontsize=fs)
                    axes[a,b].set_ylabel(scatter_tag, fontsize=fs)"""
            n += 1


        fig.suptitle("Scatter Plots Interaction Scores %s - %s value" % (session_name,scatter_tag))
        #fig.subplots_adjust(hspace=0.0)
        plt.tight_layout()
        #plt.show()
        scatter_file = '%s.scatter_plots_Iscore_%s.png' % (outfile,scatter_tag)
        #print scatter_file
        plt.savefig(scatter_file)
        #print scatter_file
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
                    axes[a,b].scatter(scatter_plots[scatter_tag][score][Class][0],scatter_plots[scatter_tag][score][Class][1],alpha = 0.05,color='red')
                    axes[a,b].set_xlabel(labels[n], fontsize=fs)
                    axes[a,b].set_ylabel(scatter_tag, fontsize=fs)
                """else:
                    axes[a,b].scatter(scatter_plots[scatter_tag][score][Class][0],scatter_plots[scatter_tag][score][Class][1],alpha = 0.05,color='blue')
                    axes[a,b].set_xlabel(labels[n], fontsize=fs)
                    axes[a,b].set_ylabel(scatter_tag, fontsize=fs)"""
            n += 1


        fig.suptitle("Scatter Plots Interaction Degrees %s - %s value" % (session_name,scatter_tag))
        #fig.subplots_adjust(hspace=0.0)
        plt.tight_layout()
        #plt.show()
        scatter_file = '%s.scatter_plots_Idegrees_%s.png' % (outfile,scatter_tag)
        #print scatter_file
        plt.savefig(scatter_file)
        #print scatter_file


def classDistributionFromFile(annotationfile,outfolder,session_name,config,by_conf=False,rin_classes=False):
    #"Uniprot-Ac\tUniprot Id\tRefseq\tPDB-ID (Input)\tResidue-Id\tAmino Acid\tPosition\tSpecies\tTag\tWeighted Surface/Core\tClass\tSimple Class\tConfidence Value\tSecondary Structure\tRecommended Structure\tSequence-ID\tCoverage\tResolution\tMax Seq Id Structure\tMax Sequence-ID\tMax Seq Id Coverage\tMax Seq Id Resolution\tAmount of mapped structures"
    outfile = '%s/%s' % (outfolder,session_name)

    t0 = time.time()

    f = open(annotationfile,'r')
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
    hc_threshs = [0.0,0.2,0.4,0.6,0.8]

    violins = {}
    comp_violins = {}

    tag_pos = None
    class_pos = None
    simple_class_pos = None
    confidence_pos = None
    rin_class_pos = None
    rin_simple_class_pos = None

    for pos,column_name in enumerate(lines[0].split('\t')):
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
        words = line.replace('\n','').split('\t')

        tag = words[tag_pos]
        if not rin_classes:
            classification = words[class_pos]
            simple_classification = words[simple_class_pos]
        else:
            classification = words[rin_class_pos]
            simple_classification = words[rin_simple_class_pos]
        if words[confidence_pos] != Output_generator.null_symbol:
            confidence = float(words[confidence_pos])
        else:
            confidence = 0.

        size += 1

        if not classification in class_map:
            class_map[classification] = 1
        else:
            class_map[classification] += 1

        if not simple_classification in simple_class_map:
            simple_class_map[simple_classification] = 1
        else:
            simple_class_map[simple_classification] += 1


        tag_dict = literal_eval(tag)
        for aa2 in tag_dict:
            tags = tag_dict[aa2]
            for tag in tags.split(','):
                if tag == '':
                    continue
                if tag[0] == '#':
                    if tag.count(':') == 1:
                        violin_tag,violin_value = tag[1:].split(':')
                    else:
                        violin_tag,violin_value = tag[1:].split('=')
                    violin_value = float(violin_value)

                    if not violin_tag in violins:
                        violins[violin_tag] = {}

                    if not simple_classification in violins[violin_tag]:
                        violins[violin_tag][simple_classification] = []
                    violins[violin_tag][simple_classification].append(violin_value)

                    if not violin_tag in comp_violins:
                        comp_violins[violin_tag] = {}

                    if not classification in comp_violins[violin_tag]:
                        comp_violins[violin_tag][classification] = []
                    comp_violins[violin_tag][classification].append(violin_value)

                    continue

                if not tag in tag_map:
                    tag_map[tag] = {}
                    simple_tag_map[tag] = {}
                    tag_sizes[tag] = 0
                if not classification in tag_map[tag]:
                    tag_map[tag][classification] = 1
                else:
                    tag_map[tag][classification] += 1
                if not simple_classification in simple_tag_map[tag]:
                    simple_tag_map[tag][simple_classification] = 1
                else:
                    simple_tag_map[tag][simple_classification] += 1
                tag_sizes[tag] += 1

        if by_conf: #TODO
            for hc_thresh in hc_threshs:
                if confidence > hc_thresh:
                    if not simple_classification in simple_high_confidence_map:
                        simple_high_confidence_map[simple_classification] = 1
                    else:
                        simple_high_confidence_map[simple_classification] += 1
                    hc_size += 1

    t1 = time.time()
    if config.verbosity >= 2:
        print(('Time for classDistribution Part1: %s' % str(t1-t0)))

    makeViolins(violins,outfile,session_name)
    makeViolins(comp_violins,outfile,session_name,add='_complex_classes')

    t2 = time.time()
    if config.verbosity >= 2:
        print(('Time for classDistribution Part2: %s' % str(t2-t1)))

    classes = list(class_map.keys())
    outlines = ['Tag\t%s' % '\t'.join(classes)]

    words = ['total']
    for classification in classes:
        r = float(class_map[classification])/float(size)
        words.append(str(r))
    outlines.append('\t'.join(words))

    if len(tag_map) > 1:
        for tag in tag_map:
            words = [tag]
            for classification in classes:
                if not classification in tag_map[tag]:
                    r = 0.0
                else:
                    r = float(tag_map[tag][classification])/float(tag_sizes[tag])
                words.append(str(r))
            outlines.append('\t'.join(words))

    simple_classes = list(simple_class_map.keys())
    simple_outlines = ['Tag\t%s' % '\t'.join(simple_classes)]

    words = ['total']
    for classification in simple_classes:
        r = float(simple_class_map[classification])/float(size)
        words.append(str(r))
    simple_outlines.append('\t'.join(words))

    if len(simple_tag_map) > 1:
        for tag in simple_tag_map:
            words = [tag]
            for classification in simple_classes:
                if not classification in simple_tag_map[tag]:
                    r = 0.0
                else:
                    r = float(simple_tag_map[tag][classification])/float(tag_sizes[tag])
                words.append(str(r))
            simple_outlines.append('\t'.join(words))

    if by_conf:
        simple_high_confidence_outlines = ['\t'.join(list(simple_high_confidence_map.keys()))]
        words = []
        for simple_classification in simple_high_confidence_map:
            r = float(simple_high_confidence_map[simple_classification])/float(hc_size)
            words.append(str(r))
        simple_high_confidence_outlines.append('\t'.join(words))

    if rin_classes:
        file_name_tag = 'rin_'
    else:
        file_name_tag = ''

    f = open('%s.%sclass_distribution.tsv' % (outfile,file_name_tag),'w')
    f.write('\n'.join(outlines))
    f.close()

    f = open('%s.%ssimple_class_distribution.tsv' % (outfile,file_name_tag),'w')
    f.write('\n'.join(simple_outlines))
    f.close()

    if by_conf:
        f = open('%s.simple_high_confidence_class_distribution.tsv' % outfile,'w')
        f.write('\n'.join(simple_high_confidence_outlines))
        f.close()

    t3 = time.time()
    if config.verbosity >= 2:
       print(('Time for c lassDistribution Part3: %s' % str(t3-t2)))

    return

def InteractionScoreAveragesFromFile(InteractionProfilesfile,outfile,session_name,by_tag=False):
    outfile = "%s/%s" % (outfile,session_name)
    f = open(InteractionProfilesfile,'r')
    lines = f.readlines()
    f.close()


    min_degree = None
    max_degree = None
    min_score = None
    max_score = None

    if len(lines) == 1:
        return

    for line in lines[1:]:
        row = line.replace('\n','').split('\t')
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


        degrees = [Ligand_Interaction_Degree,Chain_Interaction_Degree,Short_Interaction_Degree,Medium_Interaction_Degree,Long_Interaction_Degree]
        scores = [Ligand_Interaction_Score,Chain_Interaction_Score,Short_Interaction_Score,Medium_Interaction_Score,Long_Interaction_Score]
        for degree in degrees:
            if min_degree == None or degree < min_degree:
                min_degree = degree
            if max_degree == None or degree > max_degree:
                max_degree = degree

        for score in scores:
            if min_score == None or score < min_score:
                min_score = score
            if max_score == None or score > max_score:
                max_score = score


    bins = 50
    degree_bin_size = (max_degree-min_degree)/float(bins)
    score_bin_size = (max_score-min_score)/float(bins)


    degree_tag_map = {}
    score_tag_map = {}
    degree_tag_sizes = {}
    score_tag_sizes = {}

    scatter_plots = {}

    for line in lines[1:]:
        row = line.replace('\n','').split('\t')
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

        degrees = {'LI degree':Ligand_Interaction_Degree,'CI degree':Chain_Interaction_Degree,'SI degree':Short_Interaction_Degree,'MI degree':Medium_Interaction_Degree,'LoI degree':Long_Interaction_Degree}
        scores = {'LI score':Ligand_Interaction_Score,'CI score':Chain_Interaction_Score,'SI score':Short_Interaction_Score,'MI score':Medium_Interaction_Score,'LoI score':Long_Interaction_Score}

        if by_tag:
            tags = row[3]
            tags = '%s,All' % tags
            for tag in tags.split(','):
                if tag[0] == '#':
                    scatter_tag,scatter_value = tag[1:].split(':')
                    scatter_value = float(scatter_value)

                    if not scatter_tag in scatter_plots:
                        scatter_plots[scatter_tag] = {}

                        for score in scores:
                            scatter_plots[scatter_tag][score] = {}

                        for degree in degrees:
                            scatter_plots[scatter_tag][degree] = {}

                    for score in scores:
                        if not Class in scatter_plots[scatter_tag][score]:
                            scatter_plots[scatter_tag][score][Class] = [[],[]]
                    for score in scores:
                        sc = scores[score]
                        scatter_plots[scatter_tag][score][Class][0].append(sc)
                        scatter_plots[scatter_tag][score][Class][1].append(scatter_value)

                    for degree in degrees:
                        if not Class in scatter_plots[scatter_tag][degree]:
                            scatter_plots[scatter_tag][degree][Class] = [[],[]]
                    for degree in degrees:
                        deg = degrees[degree]
                        scatter_plots[scatter_tag][degree][Class][0].append(deg)
                        scatter_plots[scatter_tag][degree][Class][1].append(scatter_value)

                    continue
                for degree_name in degrees:
                    combined_tag = "%s %s" % (tag,degree_name)
                    if not combined_tag in degree_tag_map:
                        degree_tag_map[combined_tag] = [0]*bins
                    degree = degrees[degree_name]
                    for i in range(1,bins+1):
                        if degree < min_degree+i*degree_bin_size:
                            degree_tag_map[combined_tag][i-1] += 1
                            break
                    if not combined_tag in degree_tag_sizes:
                        degree_tag_sizes[combined_tag] = 1
                    else:
                        degree_tag_sizes[combined_tag] += 1
                for score_name in scores:
                    combined_tag = "%s %s" % (tag,score_name)
                    if not combined_tag in score_tag_map:
                        score_tag_map[combined_tag] = [0]*bins
                    score = scores[score_name]
                    for i in range(1,bins+1):
                        if score < min_score+i*score_bin_size:
                            score_tag_map[combined_tag][i-1] += 1
                            break
                    if not combined_tag in score_tag_sizes:
                        score_tag_sizes[combined_tag] = 1
                    else:
                        score_tag_sizes[combined_tag] += 1
        else:
            pass #TODO

    plotScatter(scatter_plots,outfile,session_name)

    #print tag_sizes
    startline_words = ['']
    for i in range(0,bins):
        startline_words.append('>%s'%str(min_degree+i*degree_bin_size))

    outlines = ['\t'.join(startline_words)]
    for tag in degree_tag_map:
        words = [tag]
        for amount in degree_tag_map[tag]:
            r = amount/float(degree_tag_sizes[tag])
            words.append(str(r))
        outlines.append('\t'.join(words))

    startline_words = ['']
    for i in range(0,bins):
        startline_words.append('>%s'%str(min_score+i*score_bin_size))

    outlines.append('\t'.join(startline_words))
    for tag in score_tag_map:
        words = [tag]
        for amount in score_tag_map[tag]:
            r = amount/float(score_tag_sizes[tag])
            words.append(str(r))
        outlines.append('\t'.join(words))

    f = open('%s.interaction_profile_means.tsv' % outfile,'w')
    f.write('\n'.join(outlines))
    f.close()


def makeDiffDict(f):
    f = open(f, "r")
    lines = f.readlines()[1:]
    f.close()
    diff_dict = {}
    for line in lines:
        words = line.split("\t")
        diff_dict[words[1]] = [words[0],words[5]]
    return diff_dict

def makeGeneDict(f):
    f = open(f, "r")
    lines = f.readlines()[1:]
    f.close()
    diff_dict = {}
    for line in lines:
        words = line.split("\t")
        diff_dict[words[0]] = words[2]
    return diff_dict

def godiffAna(fileA,fileB,papermode=False):
    go_dict_a = makeDiffDict(fileA)
    go_dict_b = makeDiffDict(fileB)
    result_list = []
    for Id in go_dict_a:
        if Id in go_dict_b:
            diff = float(go_dict_a[Id][1]) - float(go_dict_b[Id][1])
            result_list.append([Id,go_dict_a[Id][0],diff])
    result_list = sorted(result_list,key = lambda x: x[2],reverse = True)

    lines = ["GO-Term\tGO-ID\tScore-Difference"]
    for result in result_list:
        #print result
        if not papermode:
            lines.append("%s\t%s\t%s" % (result[1],result[0],str(result[2])))
        elif result[1][:2] == ' P':
            lines.append("%s%s\t%s\t%s" % (result[1][3].upper(),result[1][4:],result[0],str(result[2])))
    page = "\n".join(lines)

    name_a = (fileA.rsplit("/",1)[1]).rsplit(".",1)[0]
    name_b = (fileB.rsplit("/",1)[1]).rsplit(".",1)[0]
    outfile = "%s/GO_Diff_%s_%s.tsv" % (fileA.rsplit("/",1)[0],name_a,name_b)

    f = open(outfile, "wb")
    f.write(page)
    f.close()

def pathdiffAna(fileA,fileB,papermode=False):
    path_dict_a = makeDiffDict(fileA)
    path_dict_b = makeDiffDict(fileB)
    result_list = []
    for Id in path_dict_a:
        if Id in path_dict_b:
            diff = float(path_dict_a[Id][1]) - float(path_dict_b[Id][1])
            result_list.append([Id,path_dict_a[Id][0],diff])
    result_list = sorted(result_list,key = lambda x: x[2],reverse = True)

    lines = ["Pathway\tReactome-ID\tScore-Difference"]
    for result in result_list:
        lines.append("%s\t%s\t%s" % (result[1],result[0],str(result[2])))
    page = "\n".join(lines)

    name_a = (fileA.rsplit("/",1)[1]).rsplit(".",1)[0]
    name_b = (fileB.rsplit("/",1)[1]).rsplit(".",1)[0]
    outfile = "%s/Path_Diff_%s_%s.tsv" % (fileA.rsplit("/",1)[0],name_a,name_b)

    f = open(outfile, "wb")
    f.write(page)
    f.close()

def genediffAna(fileA,fileB):
    gene_dict_a = makeGeneDict(fileA)
    gene_dict_b = makeGeneDict(fileB)
    result_list = []
    for Id in gene_dict_a:
        if Id in gene_dict_b:
            diff = float(gene_dict_a[Id]) - float(gene_dict_b[Id])
            result_list.append([Id,diff])
    result_list = sorted(result_list,key = lambda x: x[1],reverse = True)

    lines = ["Gene\tScore-Difference"]
    for result in result_list:
        lines.append("%s\t%s" % (result[0],str(result[1])))
    page = "\n".join(lines)

    name_a = (fileA.rsplit("/",1)[1]).rsplit(".",1)[0]
    name_b = (fileB.rsplit("/",1)[1]).rsplit(".",1)[0]
    outfile = "%s/Gene_Diff_%s_%s.tsv" % (fileA.rsplit("/",1)[0],name_a,name_b)

    f = open(outfile, "wb")
    f.write(page)
    f.close()

#structure of anno_anno_map: {(Mutation_id_1,Mutation_id_2):(distance,atom,atom2,template_id_1,template_id_2,chain1,chain2)} // Mutation_id_1 < Mutation_id_2
#structure of m_aac_map: {Mutation_id:(gene_id,AAC_base)}
#structure of gn_map: {Gene_id:Uniprot_Ac}
def annoAnnoNetwork(anno_anno_map,m_aac_map,gn_map,outfile,distance_threshold = None):
    lines = []
    for (m_id1,m_id2) in anno_anno_map:
        (min_d,atom,atom2,t_id1,t_id2,chain1,chain2) = anno_anno_map[(m_id1,m_id2)]
        if distance_threshold != None:
            if min_d > distance_threshold:
                continue
        (g_id1,aac1) = m_aac_map[m_id1]
        (g_id2,aac2) = m_aac_map[m_id2]
        u_ac1 = gn_map[g_id1]
        u_ac2 = gn_map[g_id2]

        lines.append("%s:%s %s %s:%s" % (u_ac1,aac1,str(min_d),u_ac2,aac2))

    f = open(outfile,'w')
    f.write('\n'.join(lines))
    f.close()

class Output_generator:
    null_symbol = '-'

    def __init__(self,):
        self.columns = []
        self.current_line = []
        self.header_map = {}
        return

    def add_headers(self,headers):
        n = len(self.columns)
        for pos,header in enumerate(headers):
            self.header_map[header] = pos + n
            self.columns.append(header)
        return

    def get_header(self):
        header = '\t'.join(self.columns) + '\n'
        return header

    def add_value(self,header,value):
        if value == None:
            value = Output_generator.null_symbol
        if not isinstance(value,str):
            value = str(value)

        pos = self.header_map[header]
        if len(self.current_line) <= pos:
            self.current_line += [Output_generator.null_symbol]*(1+(pos-len(self.current_line)))
        self.current_line[pos] = value
        return

    def pop_line(self):
        if len(self.current_line) < len(self.columns):
            self.current_line += [Output_generator.null_symbol]*(len(self.columns)-len(self.current_line))
        line = '\t'.join(self.current_line) + '\n'
        self.current_line = []
        return line

def process_recommend_structure_str(recommended_structure_str):
    if recommended_structure_str != None and recommended_structure_str != '-':
        words = recommended_structure_str.split(';')
        if len(words) < 4:
            resolution = '-'
            cov = '-'
            seq_id = '-'
            recommended_structure = '-'
        else:
            recommended_structure,seq_id,cov,resolution = words
    else:
        resolution = '-'
        cov = '-'
        seq_id = '-'
        recommended_structure = '-'
    return recommended_structure,seq_id,cov,resolution

def appendOutput(proteins,outfolder,session_name,out_objects = None):

    if out_objects == None:
        outfile = '%s/%s' % (outfolder,session_name)
        class_file = "%s.classification.tsv" % (outfile)

        if os.path.isfile(class_file):
            os.remove(class_file)

        feature_file = "%s.features.tsv" % (outfile)
        if os.path.isfile(feature_file):
            os.remove(feature_file)

        classification_output,f = init_classification_table(class_file)
        feature_output,feat_f = init_feature_table(feature_file)
    else:
        outfile = '%s/%s' % (outfolder,session_name)
        class_file = "%s.classification.tsv" % (outfile)
        feature_file = "%s.features.tsv" % (outfile)
        f = open(class_file,'a')
        feat_f = open(feature_file,'a')
        classification_output,feature_output = out_objects

    u_acs = proteins.get_protein_u_acs()
    
    for u_ac in u_acs:
        u_id = proteins.get_u_id(u_ac)
        refseq = proteins.get_ref_ids(u_ac)

        if len(u_ac) == 6 and u_ac[4] == ':':
            input_pdb_id = u_ac
            prot_name = input_pdb_id
        else:
            prot_name = u_ac

        for pos in proteins.get_position_ids(u_ac):

            m = (u_ac,pos)

            tags = proteins.get_pos_tags(u_ac,pos)

            position = proteins.get_position(u_ac,pos)

            aa1 = position.wt_aa
            
            if position == None:
                continue

            mappings = position.mappings

            Class = mappings.Class
            conf = mappings.classification_conf
            weighted_sc = mappings.weighted_location
            recommended_structure_str = mappings.get_recommended_res_str()
            recommended_structure,seq_id,cov,resolution =process_recommend_structure_str(recommended_structure_str)
            max_seq_structure_str = mappings.get_max_seq_structure_res_str()
            max_seq_structure,max_seq_seq_id,max_seq_cov,max_seq_resolution = process_recommend_structure_str(max_seq_structure_str)

            amount_of_structures = len(mappings.qualities)
            mv_sec_ass = mappings.weighted_ssa
            simple_class = mappings.simple_class
            interaction_str = str(mappings.interaction_recommendations)

            input_res_id = proteins.get_res_id(u_ac,pos)
            input_pdb_id = ''

            if input_res_id != None:
                input_pdb_id = u_ac

            classification_output.add_value('Uniprot-Ac',u_ac)
            classification_output.add_value('Uniprot Id',u_id)
            classification_output.add_value('Refseq',refseq)
            classification_output.add_value('PDB-ID (Input)',input_pdb_id)
            classification_output.add_value('Residue-Id',input_res_id)
            classification_output.add_value('Amino Acid',aa1)
            classification_output.add_value('Position',pos)
            classification_output.add_value('Tags',tags)
            classification_output.add_value('Weighted Surface/Core',weighted_sc)
            classification_output.add_value('Class',Class)
            classification_output.add_value('Simple Class',simple_class)
            classification_output.add_value('RIN Class',mappings.rin_class)
            classification_output.add_value('RIN Simple Class',mappings.rin_simple_class)
            classification_output.add_value('Individual Interactions',interaction_str)
            classification_output.add_value('Confidence Value',conf)
            classification_output.add_value('Secondary Structure',mv_sec_ass)
            classification_output.add_value('Recommended Structure',recommended_structure)
            classification_output.add_value('Sequence-ID',seq_id)
            classification_output.add_value('Coverage',cov)
            classification_output.add_value('Resolution',resolution)
            classification_output.add_value('Max Seq Id Structure',max_seq_structure)
            classification_output.add_value('Max Sequence-ID',max_seq_seq_id)
            classification_output.add_value('Max Seq Id Coverage',max_seq_cov)
            classification_output.add_value('Max Seq Id Resolution',max_seq_resolution)
            classification_output.add_value('Amount of mapped structures',amount_of_structures)

            f.write(classification_output.pop_line())

            iupred = position.get_disorder_score()
            disorder_state = position.get_disorder_region()
            centrality_scores = mappings.get_weighted_centralities()
            rin_profile = mappings.get_weighted_profile()

            mut_aas = position.get_mut_aas()
            mut_aas.add(aa1)
            for new_aa in mut_aas:
                try:
                     mut_tags = position.get_mut_tags(new_aa)
                except:
                     mut_tags = tags
                aac = "%s%s%s" % (aa1,str(pos),new_aa)
                feature_output.add_value('Protein (Uniprot-Ac or PDB-Id:Chain-Id)',prot_name)
                feature_output.add_value('WT Amino Acid',aa1)
                feature_output.add_value('Position',pos)
                feature_output.add_value('Mut Amino Acid',new_aa)
                feature_output.add_value('AA change','%s%s' % (aa1,new_aa))
                feature_output.add_value('Tags',mut_tags)
                feature_output.add_value('Distance-based classification',Class)
                feature_output.add_value('Distance-based simple classification',simple_class)
                feature_output.add_value('RIN-based classification',mappings.rin_class)
                feature_output.add_value('RIN-based simple classification',mappings.rin_simple_class)
                feature_output.add_value('Classification confidence',conf)
                feature_output.add_value('Structure location',weighted_sc)
                feature_output.add_value('Amount of mapped structures',amount_of_structures)
                feature_output.add_value('Secondary structure assignment',mv_sec_ass)
                feature_output.add_value('IUPred value',iupred)
                feature_output.add_value('Region structure type',disorder_state)
                feature_output.add_value('Modres score',mappings.weighted_modres)
                feature_output.add_value('Phi',mappings.weighted_phi)
                feature_output.add_value('Psi',mappings.weighted_psi)

                KDmean = abs(sdsc.hydropathy[aa1]-sdsc.hydropathy[new_aa])
                feature_output.add_value('KD mean',KDmean)

                d_vol = abs(sdsc.volume[aa1]-sdsc.volume[new_aa])
                feature_output.add_value('Volume mean',d_vol)

                chemical_distance = database.getChemicalDistance(aac)
                feature_output.add_value('Chemical distance',chemical_distance)

                blosum_value = database.getBlosumValue(aac)
                feature_output.add_value('Blosum62',aac)

                aliphatic_change = int((aa1 in sdsc.aa_map_aliphatic) != (new_aa in sdsc.aa_map_aliphatic))
                hydrophobic_change = int((aa1 in sdsc.aa_map_hydrophobic) != (new_aa in sdsc.aa_map_hydrophobic))
                aromatic_change = int((aa1 in sdsc.aa_map_aromatic) != (new_aa in sdsc.aa_map_aromatic))
                positive_change = int((aa1 in sdsc.aa_map_positive) != (new_aa in sdsc.aa_map_positive))
                polar_change = int((aa1 in sdsc.aa_map_polar) != (new_aa in sdsc.aa_map_polar))
                negative_change = int((aa1 in sdsc.aa_map_negative) != (new_aa in sdsc.aa_map_negative))
                charged_change = int((aa1 in sdsc.aa_map_charged) != (new_aa in sdsc.aa_map_charged))
                small_change = int((aa1 in sdsc.aa_map_small) != (new_aa in sdsc.aa_map_small))
                tiny_change = int((aa1 in sdsc.aa_map_tiny) != (new_aa in sdsc.aa_map_tiny))
                total_change = aliphatic_change + hydrophobic_change + aromatic_change + positive_change + polar_change + negative_change + charged_change + small_change + tiny_change
                feature_output.add_value('Aliphatic change',aliphatic_change)
                feature_output.add_value('Hydrophobic change',hydrophobic_change)
                feature_output.add_value('Aromatic change',aromatic_change)
                feature_output.add_value('Positive charged change',positive_change)
                feature_output.add_value('Polar change',polar_change)
                feature_output.add_value('Negative charge change',negative_change)
                feature_output.add_value('Charged change',charged_change)
                feature_output.add_value('Small change',small_change)
                feature_output.add_value('Tiny change',tiny_change)
                feature_output.add_value('Total change',total_change)
                feature_output.add_value('B Factor',mappings.weighted_b_factor)

                if centrality_scores != None:
                    feature_output.add_value('AbsoluteCentrality',centrality_scores.AbsoluteCentrality)
                    feature_output.add_value('LengthNormalizedCentrality',centrality_scores.LengthNormalizedCentrality)
                    feature_output.add_value('MinMaxNormalizedCentrality',centrality_scores.MinMaxNormalizedCentrality)
                    feature_output.add_value('AbsoluteCentralityWithNegative',centrality_scores.AbsoluteCentralityWithNegative)
                    feature_output.add_value('LengthNormalizedCentralityWithNegative',centrality_scores.LengthNormalizedCentralityWithNegative)
                    feature_output.add_value('MinMaxNormalizedCentralityWithNegative',centrality_scores.MinMaxNormalizedCentralityWithNegative)
                    feature_output.add_value('AbsoluteComplexCentrality',centrality_scores.AbsoluteComplexCentrality)
                    feature_output.add_value('LengthNormalizedComplexCentrality',centrality_scores.LengthNormalizedComplexCentrality)
                    feature_output.add_value('MinMaxNormalizedComplexCentrality',centrality_scores.MinMaxNormalizedComplexCentrality)
                    feature_output.add_value('AbsoluteComplexCentralityWithNegative',centrality_scores.AbsoluteComplexCentralityWithNegative)
                    feature_output.add_value('LengthNormalizedComplexCentralityWithNegative',centrality_scores.LengthNormalizedComplexCentralityWithNegative)
                    feature_output.add_value('MinMaxNormalizedComplexCentralityWithNegative',centrality_scores.MinMaxNormalizedComplexCentralityWithNegative)

                feature_output.add_value('Intra_SSBOND_Propensity',mappings.weighted_intra_ssbond)
                feature_output.add_value('Inter_SSBOND_Propensity',mappings.weighted_inter_ssbond)
                feature_output.add_value('Intra_Link_Propensity',mappings.weighted_intra_link)
                feature_output.add_value('Inter_Link_Propensity',mappings.weighted_inter_link)
                feature_output.add_value('CIS_Conformation_Propensity',mappings.weighted_cis_conformation)
                feature_output.add_value('CIS_Follower_Propensity',mappings.weighted_cis_follower)
                feature_output.add_value('Inter Chain Median KD',mappings.weighted_inter_chain_median_kd)
                feature_output.add_value('Inter Chain Distance Weighted KD',mappings.weighted_inter_chain_dist_weighted_kd)
                feature_output.add_value('Inter Chain Median RSA',mappings.weighted_inter_chain_median_rsa)
                feature_output.add_value('Inter Chain Distance Weighted RSA',mappings.weighted_inter_chain_dist_weighted_rsa)
                feature_output.add_value('Intra Chain Median KD',mappings.weighted_intra_chain_median_kd)
                feature_output.add_value('Intra Chain Distance Weighted KD',mappings.weighted_intra_chain_dist_weighted_kd)
                feature_output.add_value('Intra Chain Median RSA',mappings.weighted_intra_chain_median_rsa)
                feature_output.add_value('Intra Chain Distance Weighted RSA',mappings.weighted_intra_chain_dist_weighted_rsa)
                
                for chaintype in ['mc','sc']:
                    for interaction_type in ['neighbor','short','long','ligand','ion','metal','Protein','DNA','RNA','Peptide']:
                        feature_name = '%s %s score' % (chaintype,interaction_type)
                        value = rin_profile.getChainSpecificCombiScore(chaintype,interaction_type)
                        feature_output.add_value(feature_name,value)

                        feature_name = '%s %s degree' % (chaintype,interaction_type)
                        value = rin_profile.getChainSpecificCombiDegree(chaintype,interaction_type)
                        feature_output.add_value(feature_name,value)

                        feature_name = '%s %s H-bond score' % (chaintype,interaction_type)
                        value = rin_profile.getScore(chaintype,'hbond',interaction_type)
                        feature_output.add_value(feature_name,value)

                feat_f.write(feature_output.pop_line())

    f.close()
    feat_f.close()

    return classification_output,feature_output

def init_classification_table(class_file):
    classification_output = Output_generator()
    headers = [
            'Uniprot-Ac','Uniprot Id','Refseq','PDB-ID (Input)','Residue-Id','Amino Acid','Position','Tags',
            'Weighted Surface/Core','Class','Simple Class','RIN Class','RIN Simple Class','Individual Interactions',
            'Confidence Value','Secondary Structure','Recommended Structure','Sequence-ID','Coverage','Resolution',
            'Max Seq Id Structure','Max Sequence-ID','Max Seq Id Coverage','Max Seq Id Resolution','Amount of mapped structures'
            ]
    classification_output.add_headers(headers)
    f = open(class_file,'a')
    f.write(classification_output.get_header())

    return classification_output,f

def init_feature_table(feature_file):
    feature_output = Output_generator()
    headers = [
        'Protein (Uniprot-Ac or PDB-Id:Chain-Id)','WT Amino Acid','Position','Mut Amino Acid','AA change','Tags',
        'Distance-based classification','Distance-based simple classification',
        'RIN-based classification','RIN-based simple classification',
        'Classification confidence','Structure location','Amount of mapped structures',
        'Secondary structure assignment','IUPred value','Region structure type','Modres score',
        'Phi','Psi','KD mean',
        'Volume mean','Chemical distance','Blosum62',
        'Aliphatic change','Hydrophobic change','Aromatic change','Positive charged change',
        'Polar change','Negative charge change','Charged change','Small change','Tiny change','Total change',
        'B Factor',
        'AbsoluteCentrality','LengthNormalizedCentrality','MinMaxNormalizedCentrality',
        'AbsoluteCentralityWithNegative','LengthNormalizedCentralityWithNegative','MinMaxNormalizedCentralityWithNegative',
        'AbsoluteComplexCentrality','LengthNormalizedComplexCentrality','MinMaxNormalizedComplexCentrality',
        'AbsoluteComplexCentralityWithNegative','LengthNormalizedComplexCentralityWithNegative','MinMaxNormalizedComplexCentralityWithNegative',
        'Intra_SSBOND_Propensity','Inter_SSBOND_Propensity','Intra_Link_Propensity','Inter_Link_Propensity',
        'CIS_Conformation_Propensity','CIS_Follower_Propensity',
        'Inter Chain Median KD','Inter Chain Distance Weighted KD','Inter Chain Median RSA','Inter Chain Distance Weighted RSA',
        'Intra Chain Median KD','Intra Chain Distance Weighted KD','Intra Chain Median RSA','Intra Chain Distance Weighted RSA'
        ]

    for chaintype in ['mc','sc']:
        for interaction_type in ['neighbor','short','long','ligand','ion','metal','Protein','DNA','RNA','Peptide']:
            feature_name = '%s %s score' % (chaintype,interaction_type)
            headers.append(feature_name)
            feature_name = '%s %s degree' % (chaintype,interaction_type)
            headers.append(feature_name)
            feature_name = '%s %s H-bond score' % (chaintype,interaction_type)
            headers.append(feature_name)

    feature_output.add_headers(headers)

    feat_f = open(feature_file,'a')
    feat_f.write(feature_output.get_header())
    return feature_output,feat_f

def classificationOutput(config,outfolder,session_name,session_id,ligand_filter=None):
    outfile = '%s/%s' % (outfolder,session_name)
    if config.verbosity >= 2:
        t0 = time.time()

    if ligand_filter != None:
        f = open(ligand_filter,'r')
        lines = f.readlines()
        f.close()
        ligand_filter = set()
        for line in lines:
            ligand_filter.add(line.replace('\n','').replace(' ',''))

    if config.verbosity >= 2:
        t1 = time.time()
        print("Time for classificationOutput part 1: ",t1-t0)

    table = 'RS_Mutation_Session'
    rows = ['Mutation','New_AA','Tag']
    eq_rows = {'Session':session_id}

    results = database.select(config,rows,table,equals_rows=eq_rows)

    if config.verbosity >= 2:
        t2 = time.time()
        print("Time for classificationOutput part 2: ",t2-t1)

    tag_map = {}
    for row in results:
        mut_id = row[0]
        new_aa = row[1]
        tags = row[2]
        if not mut_id in tag_map:
            tag_map[mut_id] = {new_aa:tags}
        else:
            tag_map[mut_id][new_aa] = tags

    mutation_dict = database.getMutationDict(set(tag_map.keys()),config)

    gene_id_list = set()
    for m in mutation_dict:
        gene_id_list.add(mutation_dict[m][1])

    gene_score_dict = database.getGeneScoreDict(gene_id_list,session_id,config)

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
        print("Time for classificationOutput part 3: ",t3-t2)

    columns = ['Mutation_Id','Amino_Acid_Change','Gene','Location','Class','Simple_Class','RIN_Class','RIN_Simple_Class','Interactions','Confidence',
            'Secondary_Structure','Recommended_Structure','Max_Seq_Structure','Mapped_Structures','RIN_Profile','IUPRED','IUPRED_Glob',
            'Modres_Score',
            'B_Factor','Weighted_Centrality_Scores','Weighted_Phi','Weighted_Psi','Intra_SSBOND_Propensity',
            'Inter_SSBOND_Propensity','Intra_Link_Propensity','Inter_Link_Propensity','CIS_Conformation_Propensity','CIS_Follower_Propensity',
            'Weighted_Inter_Chain_Median_KD', 'Weighted_Inter_Chain_Dist_Weighted_KD', 'Weighted_Inter_Chain_Median_RSA',
            'Weighted_Inter_Chain_Dist_Weighted_RSA', 'Weighted_Intra_Chain_Median_KD', 'Weighted_Intra_Chain_Dist_Weighted_KD',
            'Weighted_Intra_Chain_Median_RSA', 'Weighted_Intra_Chain_Dist_Weighted_RSA']

    table = 'Mutation'
    results = database.binningSelect(mutation_dict.keys(),columns,table,config)

    if config.verbosity >= 2:
        t4 = time.time()
        print("Time for classificationOutput part 4: ",t4-t3)


    classification_output,f = init_classification_table(class_file)
    feature_output,feat_f = init_feature_table(feature_file)

    stat_dict = {}
    for row in results:
        m = row[0]
        aac = row[1]
        gene_id = row[2]
        weighted_sc = row[3]
        Class =row[4]
        simple_class = row[5]
        rin_class = row[6]
        rin_simple_class = row[7]
        interaction_str = row[8]
        conf = row[9]
        mv_sec_ass = row[10]
        recommended_structure_str = row[11]
        max_seq_structure_str = row[12]
        amount_of_structures = row[13]
        rin_profile_str = row[14]
        iupred = row[15]
        disorder_state = row[16]
        modres_score = row[17]
        b_factor = row[18]
        weighted_centrality_scores = row[19]
        weighted_phi = row[20]
        weighted_psi = row[21]
        intra_ssbond_prop = row[22]
        inter_ssbond_prop = row[23]
        intra_link_prop = row[24]
        inter_link_prop = row[25]
        cis_conformation_prop = row[26]
        cis_follower_prop = row[27]
        inter_chain_median_kd = row[28]
        inter_chain_dist_weighted_kd = row[29]
        inter_chain_median_rsa = row[30]
        inter_chain_dist_weighted_rsa = row[31]
        intra_chain_median_kd = row[32]
        intra_chain_dist_weighted_kd = row[33]
        intra_chain_median_rsa = row[34]
        intra_chain_dist_weighted_rsa = row[35]

        input_res_id = mutation_dict[m][4]
        if input_res_id == None:
            input_res_id = ''

        recommended_structure,seq_id,cov,resolution = process_recommend_structure_str(recommended_structure_str)

        if disorder_state != None:
            if Class == None and (disorder_state == 'disorder' or disorder_state[0] == 'D'):
                Class = 'Disorder'
                simple_class = 'Disorder'
                rin_class = 'Disorder'
                rin_simple_class = 'Disorder'

        max_seq_structure,max_seq_seq_id,max_seq_cov,max_seq_resolution = process_recommend_structure_str(max_seq_structure_str)
        if max_seq_seq_id == '-':
            stat_dict[m] = (Class,0.)
        else:
            stat_dict[m] = (Class,float(max_seq_seq_id))

        (u_ac,refseq,u_id,error_code,error,gene_score) = gene_score_dict[gene_id]

        input_pdb_id = ''
        if len(u_ac) == 6 and u_ac[4] == ':':
            input_pdb_id = u_ac
            u_ac = ''
            prot_name = input_pdb_id
        else:
            prot_name = u_ac
        if config.skipref:
            refseq = ''
        aa1 = aac[0]
        classification_output.add_value('Uniprot-Ac',u_ac)
        classification_output.add_value('Uniprot Id',u_id)
        classification_output.add_value('Refseq',refseq)
        classification_output.add_value('PDB-ID (Input)',input_pdb_id)
        classification_output.add_value('Residue-Id',input_res_id)
        classification_output.add_value('Amino Acid',aac[0])
        classification_output.add_value('Position',aac[1:])
        classification_output.add_value('Tags',tag_map[m])
        classification_output.add_value('Weighted Surface/Core',weighted_sc)
        classification_output.add_value('Class',Class)
        classification_output.add_value('Simple Class',simple_class)
        classification_output.add_value('RIN Class',rin_class)
        classification_output.add_value('RIN Simple Class',rin_simple_class)
        classification_output.add_value('Individual Interactions',interaction_str)
        classification_output.add_value('Confidence Value',conf)
        classification_output.add_value('Secondary Structure',mv_sec_ass)
        classification_output.add_value('Recommended Structure',recommended_structure)
        classification_output.add_value('Sequence-ID',seq_id)
        classification_output.add_value('Coverage',cov)
        classification_output.add_value('Resolution',resolution)
        classification_output.add_value('Max Seq Id Structure',max_seq_structure)
        classification_output.add_value('Max Sequence-ID',max_seq_seq_id)
        classification_output.add_value('Max Seq Id Coverage',max_seq_cov)
        classification_output.add_value('Max Seq Id Resolution',max_seq_resolution)
        classification_output.add_value('Amount of mapped structures',amount_of_structures)

        f.write(classification_output.pop_line())

        centrality_scores = rin.Centrality_scores(code_str = weighted_centrality_scores)
        
        rin_profile = rin.Interaction_profile(profile_str=rin_profile_str)

        for new_aa in tag_map[m]:
            aac = "%s%s" % (aac,new_aa)
            feature_output.add_value('Protein (Uniprot-Ac or PDB-Id:Chain-Id)',prot_name)
            feature_output.add_value('WT Amino Acid',aa1)
            feature_output.add_value('Position',aac[1:-1])
            feature_output.add_value('Mut Amino Acid',new_aa)
            feature_output.add_value('AA change','%s%s' % (aa1,new_aa))
            feature_output.add_value('Tags',tag_map[m][new_aa])
            feature_output.add_value('Distance-based classification',Class)
            feature_output.add_value('Distance-based simple classification',simple_class)
            feature_output.add_value('RIN-based classification',rin_class)
            feature_output.add_value('RIN-based simple classification',rin_simple_class)
            feature_output.add_value('Classification confidence',conf)
            feature_output.add_value('Structure location',weighted_sc)
            feature_output.add_value('Amount of mapped structures',amount_of_structures)
            feature_output.add_value('Secondary structure assignment',mv_sec_ass)
            feature_output.add_value('IUPred value',iupred)
            feature_output.add_value('Region structure type',disorder_state)
            feature_output.add_value('Modres score',modres_score)
            feature_output.add_value('Phi',weighted_phi)
            feature_output.add_value('Psi',weighted_psi)

            KDmean = abs(sdsc.hydropathy[aa1]-sdsc.hydropathy[new_aa])
            feature_output.add_value('KD mean',KDmean)

            d_vol = abs(sdsc.volume[aa1]-sdsc.volume[new_aa])
            feature_output.add_value('Volume mean',d_vol)

            chemical_distance = database.getChemicalDistance(aac)
            feature_output.add_value('Chemical distance',chemical_distance)

            blosum_value = database.getBlosumValue(aac)
            feature_output.add_value('Blosum62',aac)

            aliphatic_change = int((aa1 in sdsc.aa_map_aliphatic) != (new_aa in sdsc.aa_map_aliphatic))
            hydrophobic_change = int((aa1 in sdsc.aa_map_hydrophobic) != (new_aa in sdsc.aa_map_hydrophobic))
            aromatic_change = int((aa1 in sdsc.aa_map_aromatic) != (new_aa in sdsc.aa_map_aromatic))
            positive_change = int((aa1 in sdsc.aa_map_positive) != (new_aa in sdsc.aa_map_positive))
            polar_change = int((aa1 in sdsc.aa_map_polar) != (new_aa in sdsc.aa_map_polar))
            negative_change = int((aa1 in sdsc.aa_map_negative) != (new_aa in sdsc.aa_map_negative))
            charged_change = int((aa1 in sdsc.aa_map_charged) != (new_aa in sdsc.aa_map_charged))
            small_change = int((aa1 in sdsc.aa_map_small) != (new_aa in sdsc.aa_map_small))
            tiny_change = int((aa1 in sdsc.aa_map_tiny) != (new_aa in sdsc.aa_map_tiny))
            total_change = aliphatic_change + hydrophobic_change + aromatic_change + positive_change + polar_change + negative_change + charged_change + small_change + tiny_change
            feature_output.add_value('Aliphatic change',aliphatic_change)
            feature_output.add_value('Hydrophobic change',hydrophobic_change)
            feature_output.add_value('Aromatic change',aromatic_change)
            feature_output.add_value('Positive charged change',positive_change)
            feature_output.add_value('Polar change',polar_change)
            feature_output.add_value('Negative charge change',negative_change)
            feature_output.add_value('Charged change',charged_change)
            feature_output.add_value('Small change',small_change)
            feature_output.add_value('Tiny change',tiny_change)
            feature_output.add_value('Total change',total_change)
            feature_output.add_value('B Factor',b_factor)

            feature_output.add_value('AbsoluteCentrality',centrality_scores.AbsoluteCentrality)
            feature_output.add_value('LengthNormalizedCentrality',centrality_scores.LengthNormalizedCentrality)
            feature_output.add_value('MinMaxNormalizedCentrality',centrality_scores.MinMaxNormalizedCentrality)
            feature_output.add_value('AbsoluteCentralityWithNegative',centrality_scores.AbsoluteCentralityWithNegative)
            feature_output.add_value('LengthNormalizedCentralityWithNegative',centrality_scores.LengthNormalizedCentralityWithNegative)
            feature_output.add_value('MinMaxNormalizedCentralityWithNegative',centrality_scores.MinMaxNormalizedCentralityWithNegative)
            feature_output.add_value('AbsoluteComplexCentrality',centrality_scores.AbsoluteComplexCentrality)
            feature_output.add_value('LengthNormalizedComplexCentrality',centrality_scores.LengthNormalizedComplexCentrality)
            feature_output.add_value('MinMaxNormalizedComplexCentrality',centrality_scores.MinMaxNormalizedComplexCentrality)
            feature_output.add_value('AbsoluteComplexCentralityWithNegative',centrality_scores.AbsoluteComplexCentralityWithNegative)
            feature_output.add_value('LengthNormalizedComplexCentralityWithNegative',centrality_scores.LengthNormalizedComplexCentralityWithNegative)
            feature_output.add_value('MinMaxNormalizedComplexCentralityWithNegative',centrality_scores.MinMaxNormalizedComplexCentralityWithNegative)

            feature_output.add_value('Intra_SSBOND_Propensity',intra_ssbond_prop)
            feature_output.add_value('Inter_SSBOND_Propensity',inter_ssbond_prop)
            feature_output.add_value('Intra_Link_Propensity',intra_link_prop)
            feature_output.add_value('Inter_Link_Propensity',inter_link_prop)
            feature_output.add_value('CIS_Conformation_Propensity',cis_conformation_prop)
            feature_output.add_value('CIS_Follower_Propensity',cis_follower_prop)
            feature_output.add_value('Inter Chain Median KD',inter_chain_median_kd)
            feature_output.add_value('Inter Chain Distance Weighted KD',inter_chain_dist_weighted_kd)
            feature_output.add_value('Inter Chain Median RSA',inter_chain_median_rsa)
            feature_output.add_value('Inter Chain Distance Weighted RSA',inter_chain_dist_weighted_rsa)
            feature_output.add_value('Intra Chain Median KD',intra_chain_median_kd)
            feature_output.add_value('Intra Chain Distance Weighted KD',intra_chain_dist_weighted_kd)
            feature_output.add_value('Intra Chain Median RSA',intra_chain_median_rsa)
            feature_output.add_value('Intra Chain Distance Weighted RSA',intra_chain_dist_weighted_rsa)
            
            for chaintype in ['mc','sc']:
                for interaction_type in ['neighbor','short','long','ligand','ion','metal','Protein','DNA','RNA','Peptide']:
                    feature_name = '%s %s score' % (chaintype,interaction_type)
                    value = rin_profile.getChainSpecificCombiScore(chaintype,interaction_type)
                    feature_output.add_value(feature_name,value)

                    feature_name = '%s %s degree' % (chaintype,interaction_type)
                    value = rin_profile.getChainSpecificCombiDegree(chaintype,interaction_type)
                    feature_output.add_value(feature_name,value)

                    feature_name = '%s %s H-bond score' % (chaintype,interaction_type)
                    value = rin_profile.getScore(chaintype,'hbond',interaction_type)
                    feature_output.add_value(feature_name,value)

            feat_f.write(feature_output.pop_line())

    f.close()
    feat_f.close()

    if config.verbosity >= 2:
        t5 = time.time()
        print("Time for classificationOutput part 5: ",t5-t4)

    writeStatFile(stat_file,mutation_dict,{},tag_map,stat_dict=stat_dict)

    if config.verbosity >= 2:
        t6 = time.time()
        print("Time for classificationOutput part 6: ",t6-t5)

    return class_files,[]

def writeStatFile(out_file,mutation_dict,class_dict,tag_map,stat_dict=None):
    seq_id_threshold = 0.99
    startline = 'Tag\tTotal proteins\tTotal positions\tUnmapped proteins\tEntirely disordered proteins\tProteins mapped to at least one corresponding structure (seq-id > %s%%\tProteins mapped only to structure of homologs (seq-id <= %s%%)\tMapped positions\tMapped into at least one corresponding structure (seq-id > %s%%)\tMapped only in homologs (seq-id <= %s%%)\tUnmapped, Disorder\tUnmapped, Globular' % (str(seq_id_threshold),str(seq_id_threshold),str(seq_id_threshold),str(seq_id_threshold))
    outmap = {'All':[{},0,0,0,0,0,0]}
    
    if stat_dict != None:
        class_dict = stat_dict
        max_seq_key = 1
    else:
        max_seq_key = 21
    for m in tag_map:
        for mut_aa in tag_map[m]:
            if tag_map[m][mut_aa] != None and tag_map[m][mut_aa] != '':
                raw_tags = tag_map[m][mut_aa].split(',')
            else:
                raw_tags = []
        tags = set(['All'])
        for tag in raw_tags:
            if tag[0] == '#':
                tag = tag[1:].split(':')[0]
            tags.add(tag)

        for tag in tags:
            if not tag in outmap:
                outmap[tag] = [{},0,0,0,0,0,0]
            g = mutation_dict[m][1]
            if not g in outmap[tag][0]:
                outmap[tag][0][g] = 1
            outmap[tag][1] += 1
            
            if m in class_dict:
                clas = class_dict[m][0]
                max_seq_id = class_dict[m][max_seq_key]
                if clas != 'Disorder' and clas != None:
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

        prot_numbers = [0,0,0,0]
        for g in outmap[tag][0]:
            prot_numbers[outmap[tag][0][g]] += 1

        if float(tot_pos) == 0.0:
            continue
        line = '%s\t%s\t%s\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)' % (
            tag,
            str(tot_prot),
            str(tot_pos),
            str(prot_numbers[0]),str(100.*float(prot_numbers[0])/float(tot_prot)),
            str(prot_numbers[1]),str(100.*float(prot_numbers[1])/float(tot_prot)),
            str(prot_numbers[2]),str(100.*float(prot_numbers[2])/float(tot_prot)),
            str(prot_numbers[3]),str(100.*float(prot_numbers[3])/float(tot_prot)),
            str(mapped),str(100.*float(mapped)/float(tot_pos)),
            str(mapped_to_corr),str(100.*float(mapped_to_corr)/float(tot_pos)),
            str(mapped_to_homolog),str(100.*float(mapped_to_homolog)/float(tot_pos)),
            str(dis),str(100.*float(dis)/float(tot_pos)),
            str(unmapped),str(100.*float(unmapped)/float(tot_pos)))
        lines.append(line)
    f = open(out_file,'w')
    f.write("\n".join(lines))
    f.close()

#called by structman
def main(sess_id,output_path,config,intertable=False):
    db_name = config.db_name
    db_adress = config.db_adress
    db_password = config.db_password
    db_user_name = config.db_user_name
    go = config.go
    godiff = config.godiff
    classification = config.classification
    path = config.path
    pathdiff = config.pathdiff
    do_modelling = config.do_modelling
    ligand_file = config.ligand_file
    multi_modelling = config.multi_modelling
    mod_per_mut = config.mod_per_gene
    mod_per_gene = config.mod_per_gene
    infile = ''
    tanimoto_cutoff = config.tanimoto_cutoff
    distance_threshold = config.milieu_threshold
    ligand_filter = config.ligand_filter
    proteome = config.proteome
    proc_n = config.proc_n
    verbose = config.verbose

    intertable_conf = config.intertable_conf
    intertable = intertable or intertable_conf


    if sess_id == 0 or sess_id == None:
        session_id = database.getSessionId(infile,config)
    else:
        session_id =sess_id

    db,cursor = config.getDB()

    if infile == '':
        infile = database.getSessionFile(session_id,db,cursor)

    db.close()

    session_name = (infile.rsplit("/",1)[1]).rsplit(".",1)[0]

    t0 = time.time()

    if classification:
        t00 = time.time()
        classfiles,interfiles = classificationOutput(config,output_path,session_name,session_id)
        t01 = time.time()
        if config.verbosity >= 2:
            print("Time for classificationOutput: ",t01-t00)
        for classfile in classfiles:
            classDistributionFromFile(classfile,output_path,session_name,config)
            classDistributionFromFile(classfile,output_path,session_name,config,rin_classes=True)
        t02 = time.time()
        if config.verbosity >= 2:
            print("Time for producing classification distributions: ",t02-t01)
        
        if intertable:
            for interfile in interfiles:
                InteractionScoreAveragesFromFile(interfile,output_path,session_name,by_tag=True)
            t03 = time.time()
            if config.verbosity >= 2:
                print("Time for producing Interaction files: ",t03-t02)
    t1 = time.time()
    if config.verbosity >= 2:
        print("Time for producing classification file: ",t1-t0)    

    db,cursor = config.getDB()

    if go:
        database.goTermAnalysis(session_id,"%s/%s.goterm.tsv" % (output_path,session_name),db,cursor)
        if godiff:
            files = os.listdir(output_path)
            go_files = []
            for f in files:
                if '.goterm.tsv' in f:
                    go_files.append(f)
            print(go_files)
            print(files)
            if len(go_files) == 2:
                fileA = "%s/%s" % (output_path,go_files[0])
                fileB = "%s/%s" % (output_path,go_files[1])
                godiffAna(fileA,fileB)
    if path:
        database.pathwayAnalysis(session_id,"%s/%s.pathway.tsv" % (output_path,session_name),db,cursor)
        if pathdiff:
            files = os.listdir(output_path)
            path_files = []
            for f in files:
                if '.pathway.tsv' in f:
                    path_files.append(f)
            if len(path_files) == 2:
                fileA = "%s/%s" % (output_path,path_files[0])
                fileB = "%s/%s" % (output_path,path_files[1])
                pathdiffAna(fileA,fileB)
    if do_modelling:
        modelling.massModel(session_id,db,cursor,output_path,total_models = 0,model_per_gene = int(mod_per_gene),multiple_mutations = multi_modelling)

    cursor.close()
    db.close()
    if ligand_file == None:
        try:
            ligand_file_names = os.listdir("%s/ligands" % infile.rsplit("/",1)[0])
            ligand_files = []
            for ligand_file_name in ligand_file_names:
                ligand_files.append("%s/ligands/%s" % (infile.rsplit("/",1)[0],ligand_file_name))
        except:
            ligand_files = []
    else:
        ligand_files = [ligand_file]
    for ligand_file in ligand_files:
        t0 = time.time()
        anno_dict = babel.ligandAnalyzer(ligand_file,session_id,db_name,db_adress,db_user_name,db_password,cutoff=tanimoto_cutoff,distance_threshold=distance_threshold)
        t1 = time.time()
        babel.writeReport(anno_dict,"%s/Ligand_Report_%s_%s.tsv" % (output_path,ligand_file.rsplit('/',1)[1].rsplit(".",1)[0],session_name),db_name,db_adress,db_user_name,db_password)
        t2 = time.time()
        if config.verbosity >= 2:
            print("Time for ligandAnalyzer: ",t1-t0)
            print("Time for writeReport: ",t2-t1)
