import MySQLdb
import database
import os
import babel
import time

db_adress = ""
db_user_name = ""
db_password = ""
db_name = ""
output_path = ""
conf_path = ""
infile = ""
sess_id = 0
go = False
anno = False
classification=True
gene = False
path = False
godiff = False
pathdiff = False
do_modelling = False
multi_modelling = False
ligand_file = None
mod_per_mut = 0
mod_per_gene = 0
tanimoto_cutoff = 0.05
distance_threshold = 10.0
ligand_filter = None
proteome = False
intertable_conf=False

def classDistributionFromFile(annotationfile,outfolder,session_name,by_conf=False):
    #Uniprot	AAC	Weighted Surface/Core	Complex Class	Simple Class	Confidence Value	Secondary Structure	Chemical Distance	Blosum62 Value	Uniprot Id
    outfile = '%s/%s' % (outfolder,session_name)

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

    for line in lines[1:]:
        words = line.replace('\n','').split('\t')
        species = words[2]
        tag = words[3]
        classification = words[5]
        simple_classification = words[6]
        confidence = float(words[7])

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

    classes = class_map.keys()
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

    simple_classes = simple_class_map.keys()
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
        simple_high_confidence_outlines = ['\t'.join(simple_high_confidence_map.keys())]
        words = []
        for simple_classification in simple_high_confidence_map:
            r = float(simple_high_confidence_map[simple_classification])/float(hc_size)
            words.append(str(r))
        simple_high_confidence_outlines.append('\t'.join(words))

    f = open('%s.class_distribution.tsv' % outfile,'w')
    f.write('\n'.join(outlines))
    f.close()

    f = open('%s.simple_class_distribution.tsv' % outfile,'w')
    f.write('\n'.join(simple_outlines))
    f.close()

    if by_conf:
        f = open('%s.simple_high_confidence_class_distribution.tsv' % outfile,'w')
        f.write('\n'.join(simple_high_confidence_outlines))
        f.close()

def InteractionScoreAveragesFromFile(InteractionProfilesfile,outfile,by_tag=False):
    
    f = open(InteractionProfilesfile,'r')
    lines = f.readlines()
    f.close()


    min_degree = None
    max_degree = None
    min_score = None
    max_score = None

    for line in lines[1:]:
        row = line.replace('\n','').split('\t')
        Ligand_Interaction_Degree = float(row[2])
        Ligand_Interaction_Score = float(row[3])
        Chain_Interaction_Degree = float(row[4])
        Chain_Interaction_Score = float(row[5])
        Short_Interaction_Degree = float(row[6])
        Short_Interaction_Score = float(row[7])
        Medium_Interaction_Degree = float(row[8])
        Medium_Interaction_Score = float(row[9])
        Long_Interaction_Degree = float(row[10])
        Long_Interaction_Score = float(row[11])


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

    for line in lines[1:]:
        row = line.replace('\n','').split('\t')
        Ligand_Interaction_Degree = float(row[2])
        Ligand_Interaction_Score = float(row[3])
        Chain_Interaction_Degree = float(row[4])
        Chain_Interaction_Score = float(row[5])
        Short_Interaction_Degree = float(row[6])
        Short_Interaction_Score = float(row[7])
        Medium_Interaction_Degree = float(row[8])
        Medium_Interaction_Score = float(row[9])
        Long_Interaction_Degree = float(row[10])
        Long_Interaction_Score = float(row[11])

        degrees = {'LI degree':Ligand_Interaction_Degree,'CI degree':Chain_Interaction_Degree,'SI degree':Short_Interaction_Degree,'MI degree':Medium_Interaction_Degree,'LoI degree':Long_Interaction_Degree}
        scores = {'LI score':Ligand_Interaction_Score,'CI score':Chain_Interaction_Score,'SI score':Short_Interaction_Score,'MI score':Medium_Interaction_Score,'LoI score':Long_Interaction_Score}

        if by_tag:
            tag = row[12]
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

    f = open(outfile,'w')
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


def config(config_path):
    global go
    global godiff
    global anno
    global gene
    global path
    global pathdiff
    global do_modelling
    global multi_modelling
    global mod_per_mut
    global mod_per_gene
    global db_name
    global db_adress
    global db_password
    global db_user_name
    global tanimoto_cutoff
    global distance_threshold
    global ligand_filter
    global proteome
    global intertable_conf
    
    f = open(config_path, "r")
    lines = f.read().split('\n')
    f.close()

    out_options = []
    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == '#':
            continue
        words = line.split("=")
        if len(words) > 1:
            opt = words[0]
            arg = words[1].replace("\n","")
            if opt == 'do_anno':
                if arg == "True":
                    anno = True
                else:
                    anno = False
            elif opt == 'do_genesort':
                if arg == "True":
                    gene = True
                else:
                    gene = False
            elif opt == 'do_goterm':
                if arg == "True":
                    go = True
                else:
                    go = False
            elif opt == 'do_godiff':
                if arg == "True":
                    godiff = True
                else:
                    godiff = False
            elif opt == 'do_pathway':
                if arg == "True":
                    path = True
                else:
                    path = False
            elif opt == 'do_pathdiff':
                if arg == "True":
                    pathdiff = True
                else:
                    pathdiff = False
            elif opt == 'multi_modelling':
                if arg == "True":
                    multi_modelling = True
                else:
                    multi_modelling = False
            elif opt == 'mod_per_gene':
                mod_per_gene = arg
            elif opt == 'mod_per_mut':
                mod_per_mut = arg
            elif opt == 'do_modelling':
                if arg == "True":
                    do_modelling = True
                else:
                    do_modelling = False
            elif opt == 'db_adress':
                db_adress = arg
            elif opt == 'db_user_name':
                db_user_name = arg
            elif opt == 'db_password':
                db_password = arg
            elif opt == 'db_name':
                db_name = arg
            elif opt == 'tanimoto_cutoff':
                tanimoto_cutoff = float(arg)
            elif opt == 'lig_dist_thresh':
                distance_threshold = arg
            elif opt == 'ligand_filter':
                ligand_filter = arg
            elif opt == 'proteome':
                if arg == 'True':
                    proteome = True
            elif opt == 'intertable':
                if arg == 'True':
                    intertable_conf = True

#called by structman
def main(db_name,db_adress,db_password,db_user_name,sess_id,output_path,config_path='',overwrite=False,anno=False,intertable=False):
    if config_path != '':
        config(config_path)
    global go
    global godiff
    global classification
    global gene
    global path
    global pathdiff
    global do_modelling
    global ligand_file
    global multi_modelling
    global mod_per_mut
    global mod_per_gene
    global infile
    global tanimoto_cutoff
    global distance_threshold
    global ligand_filter
    global proteome

    global intertable_conf
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

    database.updateGeneScores(session_id,db,cursor)

    if classification:
        classfiles = database.minDistOut(output_path,session_name,session_id,db,cursor,overwrite=overwrite,intertable=intertable)
        for classfile in classfiles:
            classDistributionFromFile(classfile,output_path,session_name)

    t0 = time.time()
    if anno:
        (ofile,reduced) = database.prodAnoOut("%s/%s.anotations.tsv" % (output_path,session_name),session_id,db_name,db_adress,db_user_name,db_password,ligand_filter=ligand_filter,proteome=proteome)
        if reduced:
            os.rename("%s/%s.anotations.tsv" % (output_path,session_name),"%s/%s.R.anotations.tsv" % (output_path,session_name))
    t1 = time.time()
    print "Time for prodAnoOut2: ",t1-t0    
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
            print go_files
            print files
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

        print "Time for ligandAnalyzer: ",t1-t0
        print "Time for writeReport: ",t2-t1
