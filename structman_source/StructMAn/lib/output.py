import pymysql as MySQLdb
import database
import os
import babel
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing

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


def classDistributionFromFile(annotationfile,outfolder,session_name,by_conf=False,rin_classes=False):
    #"Uniprot-Ac\tUniprot Id\tRefseq\tPDB-ID (Input)\tResidue-Id\tAmino Acid\tPosition\tSpecies\tTag\tWeighted Surface/Core\tClass\tSimple Class\tConfidence Value\tSecondary Structure\tRecommended Structure\tSequence-ID\tCoverage\tResolution\tMax Seq Id Structure\tMax Sequence-ID\tMax Seq Id Coverage\tMax Seq Id Resolution\tAmount of mapped structures"
    outfile = '%s/%s' % (outfolder,session_name)

    t0 = time.time()

    f = open(annotationfile,'r')
    lines = f.readlines()
    f.close()

    tag_map = {}
    simple_tag_map = {}
    tag_sizes = {}

    species_map = {}
    simple_species_map = {}
    species_sizes = {}

    tag_species_map = {}
    simple_tag_species_map = {}
    tag_species_sizes = {}

    class_map = {}
    simple_class_map = {}
    size = 0

    simple_high_confidence_map = {}
    hc_size = 0
    hc_threshs = [0.0,0.2,0.4,0.6,0.8]

    violins = {}
    comp_violins = {}

    species_pos = None
    tag_pos = None
    class_pos = None
    simple_class_pos = None
    confidence_pos = None
    rin_class_pos = None
    rin_simple_class_pos = None

    for pos,column_name in enumerate(lines[0].split('\t')):
        if column_name == 'Species':
            species_pos = pos
        if column_name == 'Tag':
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
        species = words[species_pos]
        tag = words[tag_pos]
        if not rin_classes:
            classification = words[class_pos]
            simple_classification = words[simple_class_pos]
        else:
            classification = words[rin_class_pos]
            simple_classification = words[rin_simple_class_pos]
        if words[confidence_pos] != 'None':
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

        if not species in species_map:
            species_map[species] = {}
            simple_species_map[species] = {}
            species_sizes[species] = 0
        if not classification in species_map[species]:
            species_map[species][classification] = 1
        else:
            species_map[species][classification] += 1
        if not simple_classification in simple_species_map[species]:
            simple_species_map[species][simple_classification] = 1
        else:
            simple_species_map[species][simple_classification] += 1
        species_sizes[species] += 1


        tags = tag.split(',')
        for tag in tags:
            if tag[0] == '#':
                violin_tag,violin_value = tag[1:].split(':')

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
            
            if not (species,tag) in tag_species_map:
                tag_species_map[(species,tag)] = {}
                simple_tag_species_map[(species,tag)] = {}
                tag_species_sizes[(species,tag)] = 0
            if not classification in tag_species_map[(species,tag)]:
                tag_species_map[(species,tag)][classification] = 1
            else:
                tag_species_map[(species,tag)][classification] += 1
            if not simple_classification in simple_tag_species_map[(species,tag)]:
                simple_tag_species_map[(species,tag)][simple_classification] = 1
            else:
                simple_tag_species_map[(species,tag)][simple_classification] += 1
            tag_species_sizes[(species,tag)] += 1


        if by_conf: #TODO
            for hc_thresh in hc_threshs:
                if confidence > hc_thresh:
                    if not simple_classification in simple_high_confidence_map:
                        simple_high_confidence_map[simple_classification] = 1
                    else:
                        simple_high_confidence_map[simple_classification] += 1
                    hc_size += 1

    t1 = time.time()
    print(('Time for classDistribution Part1: %s' % str(t1-t0)))

    makeViolins(violins,outfile,session_name)
    makeViolins(comp_violins,outfile,session_name,add='_complex_classes')

    t2 = time.time()
    print(('Time for classDistribution Part2: %s' % str(t2-t1)))

    classes = list(class_map.keys())
    outlines = ['\t%s' % '\t'.join(classes)]

    words = ['total']
    for classification in classes:
        r = float(class_map[classification])/float(size)
        words.append(str(r))
    outlines.append('\t'.join(words))

    if len(species_map) > 1:
        for species in species_map:
            words = [species]
            for classification in classes:
                if not classification in species_map[species]:
                    r = 0.0
                else:
                    r = float(species_map[species][classification])/float(species_sizes[species])
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

    if len(tag_map) > 1 and len(species_map) > 1:
        for (species,tag) in tag_species_map:
            words = ['%s %s' % (species,tag)]
            for classification in classes:
                if not classification in tag_species_map[(species,tag)]:
                    r = 0.0
                else:
                    r = float(tag_species_map[(species,tag)][classification])/float(tag_species_sizes[(species,tag)])
                words.append(str(r))
            outlines.append('\t'.join(words))

    simple_classes = list(simple_class_map.keys())
    simple_outlines = ['\t%s' % '\t'.join(simple_classes)]

    words = ['total']
    for classification in simple_classes:
        r = float(simple_class_map[classification])/float(size)
        words.append(str(r))
    simple_outlines.append('\t'.join(words))

    if len(simple_species_map) > 1:
        for species in simple_species_map:
            words = [species]
            for classification in simple_classes:
                if not classification in simple_species_map[species]:
                    r = 0.0
                else:
                    r = float(simple_species_map[species][classification])/float(species_sizes[species])
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

    if len(simple_tag_map) > 1 and len(simple_species_map) > 1:
        for (species,tag) in simple_tag_species_map:
            words = ['%s %s' % (species,tag)]
            for classification in simple_classes:
                if not classification in simple_tag_species_map[(species,tag)]:
                    r = 0.0
                else:
                    r = float(simple_tag_species_map[(species,tag)][classification])/float(tag_species_sizes[(species,tag)])
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
    print(('Time for classDistribution Part3: %s' % str(t3-t2)))

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


#called by structman
def main(sess_id,output_path,config,overwrite=False,anno=False,intertable=False):
    db_name = config.db_name
    db_adress = config.db_adress
    db_password = config.db_password
    db_user_name = config.db_user_name
    go = config.go
    godiff = config.godiff
    classification = config.classification
    gene = config.gene
    path = config.path
    pathdiff = config.pathdiff
    do_modelling = config.do_modelling
    ligand_file = config.ligand_file
    multi_modelling = config.multi_modelling
    mod_per_mut = config.mod_per_gene
    mod_per_gene = config.mod_per_gene
    infile = ''
    tanimoto_cutoff = config.tanimoto_cutoff
    distance_threshold = config.distance_threshold
    ligand_filter = config.ligand_filter
    proteome = config.proteome
    proc_n = config.proc_n
    verbose = config.verbose

    intertable_conf = config.intertable_conf
    intertable = intertable or intertable_conf

    db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    cursor = db.cursor()
    if sess_id == 0 or sess_id == None:
        session_id = database.getSessionId(infile,db,cursor)
    else:
        session_id =sess_id

    if infile == '':
        infile = database.getSessionFile(session_id,db,cursor)

    session_name = (infile.rsplit("/",1)[1]).rsplit(".",1)[0]

    #database.updateGeneScores(session_id,db,cursor)

    t0 = time.time()

    if classification:
        t00 = time.time()
        if intertable:
            classfiles,interfiles = database.minDistOut(output_path,session_name,session_id,db,cursor,overwrite=overwrite,intertable=intertable,processes=proc_n,verbose=verbose)
        else:
            classfiles = database.classificationOutput(output_path,session_name,session_id,db,cursor,overwrite=overwrite,verbose=verbose)
        t01 = time.time()
        print("Time for minDistOut: ",t01-t00)
        for classfile in classfiles:
            classDistributionFromFile(classfile,output_path,session_name)
            classDistributionFromFile(classfile,output_path,session_name,rin_classes=True)
        t02 = time.time()
        print("Time for producing classification distributions: ",t02-t01)
        
        if intertable:
            for interfile in interfiles:
                InteractionScoreAveragesFromFile(interfile,output_path,session_name,by_tag=True)
            t03 = time.time()
            print("Time for producing Interaction files: ",t03-t02)
    t1 = time.time()
    print("Time for producing classification file: ",t1-t0)    

    t0 = time.time()
    if anno:
        ofile = database.prodAnoOut("%s/%s.anotations.tsv" % (output_path,session_name),session_id,db_name,db_adress,db_user_name,db_password,ligand_filter=ligand_filter,proteome=proteome)
    t1 = time.time()
    print("Time for producing annotation file: ",t1-t0)    
    if gene:
        database.sortGenes(session_id,"%s/%s.protsort.tsv" % (output_path,session_name),db,cursor)
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

        print("Time for ligandAnalyzer: ",t1-t0)
        print("Time for writeReport: ",t2-t1)
