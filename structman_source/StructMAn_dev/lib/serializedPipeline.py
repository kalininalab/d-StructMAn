import MySQLdb
import os
import sys
import getopt
import time
from multiprocessing import Process, Queue, Manager, Value, Lock
import shutil
import traceback
import math

import pdbParser as pdb
import templateSelection
import templateFiltering
import globalAlignment
import database
import uniprot
import blast
import annovar


#dev options:
anno_anno_calculation = False
full_anno_anno_calculation = False
error_annotations_into_db = True
anno_session_mapping = True
calculate_interaction_profiles=True

proc_n = 32
blast_processes = proc_n
alignment_processes = proc_n
annotation_processes = proc_n
number_of_processes = proc_n

dssp = True

#active options
option_seq_thresh = 35.0
option_res_thresh = 4.5
option_ral_thresh = 0.5
option_seq_wf = 1.0
option_ral_wf = 0.5
option_res_wf = 0.25
option_rval_wf = 0.1
option_lig_wf = 1.0
option_chain_wf = 1.0
option_number_of_templates = 0
compute_surface = True
tax_id = None
ref_genome_id = None
mrna_fasta = None
species = None
config = ""
auto_scale = False
auto_scale_template_threshold = 10
#auto_scale_template_percentage gets set by buildQueue
auto_scale_template_percentage = None

#Paths to Helper Programs     
blast_path = ""
blast_db_path = ""
output_path = ""
pdb_path = ""
annovar_path = ""
smiles_path = ""
inchi_path = ""
human_id_mapping_path = ''
human_proteome_path = ''
dssp_path = ''
rin_db_path = ''
errorlog = ""


#Database
db_adress = ""
db_user_name = ""
db_password = ""
db_name = ""
session = 0

manager = Manager()
lock = manager.Lock()
cwd = ''

def SQLDateTime():
    return time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())

def setBasePaths(main_file_path):
    global blast_db_path
    global smiles_path
    global inchi_path
    global human_id_mapping_path
    global human_proteome_path

    trunk = main_file_path.rsplit('/',1)[0]

    blast_db_path = '%s/lib/base/blast_db/pdbba' % trunk
    smiles_path = '%s/lib/base/ligand_bases/Components-smiles-stereo-oe.smi' % trunk
    inchi_path = '%s/lib/base/ligand_bases/inchi_base.tsv' % trunk
    human_id_mapping_path = '%s/lib/base/id_mapping' % trunk
    human_proteome_path = '%s/lib/base/proteomes/proteom_isoforms.fasta' % trunk

def setPaths(output_path,config):
    global blast_path
    global pdb_path
    global annovar_path
    global dssp_path
    global rin_db_path
    global errorlog

    cwd = os.getcwd()

    f = open(config, "r")
    lines = f.read().split('\n')
    f.close()
    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == '#':
            continue
        words = line.split("=")
        if len(words) < 2:
            continue
        opt = words[0]
        arg = words[1].replace("\n","")
        if opt == 'blast_path':
            blast_path = arg
        elif opt == 'pdb':
            pdb_path = arg
        elif opt == 'annovar_path':
            annovar_path = arg
        elif opt == 'dssp_path':
            dssp_path = arg
        elif opt == 'rin_db_path':
            rin_db_path = arg

    if not os.path.exists("%s/errorlogs" % output_path):
        os.mkdir("%s/errorlogs" % output_path)
    errorlog = "%s/errorlogs/errorlog.txt" % output_path

def parseConfig(config):
    global option_seq_thresh
    global option_res_thresh
    global option_ral_thresh
    global option_seq_wf
    global option_ral_wf
    global option_res_wf
    global option_rval_wf
    global option_lig_wf
    global option_chain_wf
    global option_number_of_templates
    global compute_surface
    global tax_id
    global species
    global ref_genome_id
    global db_adress
    global db_user_name
    global db_password
    global db_name
    global auto_scale
    global error_annotations_into_db
    global anno_session_mapping
    global mrna_fasta
    global anno_anno_calculation
    global calculate_interaction_profiles

    f = open(config, "r")
    lines = f.read().split('\n')
    f.close()
    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == '#':
            continue
        words = line.split("=")
        if len(words) > 1:
            opt = words[0]
            arg = words[1].replace("\n","")
            if opt == 'seq_thresh':
                option_seq_thresh = float(arg)
                if option_seq_thresh <= 1.0:
                    option_seq_thresh *= 100.0
            elif opt == 'res_thresh':
                option_res_thresh = float(arg)
            elif opt == 'cov_thresh':
                option_ral_thresh = float(arg)
                if option_ral_thresh > 1.0:
                    option_ral_thresh *= 0.01
            elif opt == 'seq_wf':
                option_seq_wf = float(arg)
            elif opt == 'cov_wf':
                option_ral_wf = float(arg)
            elif opt == 'res_wf':
                option_res_wf = float(arg)
            elif opt == 'rval_wf':
                option_rval_wf = float(arg)
            elif opt == 'lig_wf':
                option_lig_wf = float(arg)
            elif opt == 'chain_wf':
                option_chain_wf = float(arg)
            elif opt == 'option_number_of_templates':
                option_number_of_templates = int(arg)
            elif opt == 'db_adress':
                db_adress = arg
            elif opt == 'db_user_name':
                db_user_name = arg
            elif opt == 'db_password':
                db_password = arg
            elif opt == 'db_name':
                db_name = arg
            elif opt == 'species':
                species = arg
            elif opt == 'tax_id':
                tax_id = arg
            elif opt == 'ref':
                ref_genome_id = arg
            elif opt == 'compute_surface':
                compute_surface = True
            elif opt == 'auto_scale':
                if arg == 'on':
                    auto_scale = True
                else:
                    auto_scale = False
            elif opt == 'error_annotations_into_db':
                if arg == 'True':
                    error_annotations_into_db = True
                elif arg == 'False':
                    error_annotations_into_db = False
            elif opt == 'anno_session_mapping':
                if arg == 'True':
                    anno_session_mapping = True
                elif arg == 'False':
                    anno_session_mapping = False
            elif opt == 'mrna':
                mrna_fasta = arg
            elif opt == 'anno_anno_calculation':
                if arg == 'True':
                    anno_anno_calculation = True
            elif opt == 'calculate_interaction_profiles':
                if arg == 'True':
                    calculate_interaction_profiles = True
                elif arg == 'False':
                    calculate_interaction_profiles = False

def buildQueue(filename,junksize,mrna_fasta=None):
    global auto_scale_template_percentage
    global species
    global anno_anno_calculation
    global human_id_mapping_path
    global db
    global cursor

    genes = {}
    u_ids = set()
    u_acs = set()
    nps = set()
    id_map = {}
    ac_map = {}
    np_map = {}

    pdb_map = {}

    f = open(filename, "r")
    lines = f.read().split('\n')
    f.close()
    ac_id_map = {}
    tag_map = {}
    species_map = {}
    fasta_map = {}
    global_species = species

    for line in lines:
        #print line
        if line == '':
            continue
        if len(line) < 3:
            print "Skipped input line:\n%s\nToo short.\n" % line
            continue
        line = line.replace(' ','\t')
        words = line.split("\t")
        if len(words) < 2:
            print "Skipped input line:\n%s\nToo few words.\n" % line
            continue
        sp_id = words[0]#.replace("'","\\'")
        try:
            aachange = words[1].replace("\n","")
            if not sp_id.count(':') == 1:
                if ord(aachange[-1]) > 47 and ord(aachange[-1]) < 58:
                    aachange = "%s%s" % (aachange,aachange[0])
                pos = int(aachange.split(',')[0][1:-1])
                if pos < 1:
                    print "Skipped input line:\n%s\nPosition < 1.\n" % line
                    continue
            else:
                aachange = "%s%s" % (aachange,aachange[0])

        except:
            print "File Format Error: ",line
            sys.exit()

        tag = None
        species = global_species

        if len(words) > 2:
            species = words[2].replace(' ','')
            if species == '':
                species = None
        if len(words) > 3:
            tag = words[3]

        if mrna_fasta == None:
            if sp_id[0:3] == "NP_" or sp_id[0:3] == "XP_":
                nps.add(sp_id)
                if not sp_id in np_map:
                    np_map[sp_id] = set([aachange])
                else:
                    np_map[sp_id].add(aachange)
                
            else:
                if sp_id.count('_') > 0:
                    u_ids.add(sp_id)
                    if not sp_id in id_map:
                        id_map[sp_id] = set([aachange])
                    else:
                        id_map[sp_id].add(aachange)
                elif len(sp_id) == 6 and sp_id.count(':') == 1:
                    pdb_chain_tuple = sp_id
                    if not pdb_chain_tuple in pdb_map:
                        pdb_map[pdb_chain_tuple] = set([aachange])
                    else:
                        pdb_map[pdb_chain_tuple].add(aachange)
                else:
                    u_acs.add(sp_id)
                    if not sp_id in ac_map:
                        ac_map[sp_id] = set([aachange])
                    else:
                        ac_map[sp_id].add(aachange)
        else:
            if species != 'multi':
                sp_id = '%s_%s' % (sp_id,species)
            if not sp_id in genes:
                genes[sp_id] = set([aachange])
            else:
                genes[sp_id].add(aachange)
            ac_id_map[sp_id] = '-'
            if not species in fasta_map:
                fasta_map[species] = '%s/%s.fa' % (mrna_fasta,species)

        if tag != None:
            if not sp_id in tag_map:
                tag_map[sp_id] = {}
            if not aachange in tag_map[sp_id]:
                tag_map[sp_id][aachange] = tag
            else:
                old_tags = set(tag_map[sp_id][aachange].split(','))
                new_tags = set(tag.split(','))
                tag_map[sp_id][aachange] = ','.join(old_tags | new_tags)

        if species != None:
            if sp_id in species_map:
                if species != species_map[sp_id]:
                    print "Error: Identical gene ids for different species are not allowed: ",sp_id,species,species_map[sp_id]
            species_map[sp_id] = species

    if species == 'human':
        ac_id_map,u_ac_refseq_map,ac_map = uniprot.humanIdMapping(u_acs,ac_map,human_id_mapping_path,obsolete_check=False)
        id_ac_map = {} #TODO
        np_ac_map = {} #TODO
        np_id_map = {} #TODO

    else:
        u_ac_refseq_map = {}
        u_size = len(u_acs)
        if u_size > 0:
            #print len(u_acs)
            uni_junksize = 5000
            if u_size > 5000:
                n_of_uni_batches = u_size/uni_junksize
                if u_size%uni_junksize != 0:
                    n_of_uni_batches += 1
                uni_batchsize = u_size/n_of_uni_batches
                uni_rest = u_size%n_of_uni_batches
                ac_id_map = {}
                for i in range(0,n_of_uni_batches):
                    batch = set([])
                    if i == n_of_uni_batches - 1: # the last batch
                        uni_batchsize += uni_rest
                    for j in range(0,uni_batchsize):
                        u_ac = u_acs.pop()
                        batch.add(u_ac)
                    ac_id_map.update(uniprot.getUniprotIds(batch,"ACC",target_type="ID"))
            else:
                ac_id_map = uniprot.getUniprotIds(u_acs,"ACC",target_type="ID")
        elif mrna_fasta == None:
            ac_id_map = {}

        if len(u_ids) > 0:
            id_ac_map = uniprot.getUniprotIds(u_ids,"ID",target_type="ACC")
        else:
            id_ac_map = {}

        if len(nps) > 0:
            np_ac_map = uniprot.getUniprotIds(nps,"P_REFSEQ_AC",target_type="ACC")
            np_id_map = uniprot.getUniprotIds(nps,"P_REFSEQ_AC",target_type="ID")
        else:
            np_ac_map = {}
            np_id_map = {}

    #print species,mrna_fasta,u_size
    #print id_ac_map
    #print np_ac_map

    if mrna_fasta == None:
        for u_ac in ac_id_map:
            genes[u_ac] = ac_map[u_ac]

        for pdb_chain_tuple in pdb_map:
            genes[pdb_chain_tuple] = pdb_map[pdb_chain_tuple]
            ac_id_map[pdb_chain_tuple] = pdb_chain_tuple

        for u_id in id_ac_map:
            u_ac = id_ac_map[u_id]
            if not u_ac in genes:
                genes[u_ac] = id_map[u_id]
            else:
                genes[u_ac] = genes[u_ac] | id_map[u_id]
            ac_id_map[u_ac] = u_id

            #replace the entries in the tag listed by u-ids with u-acs
            if u_id in tag_map:
                if not u_ac in tag_map:
                    tag_map[u_ac] = {}
                tag_map[u_ac].update(tag_map[u_id])
                del tag_map[u_id]

            if u_id in species_map:
                species_map[u_ac] = species_map[u_id]
                del species_map[u_id]

        for np in np_ac_map:
            u_ac = np_ac_map[np]
            if not u_ac in genes:
                genes[u_ac] = np_map[np]
            else:
                genes[u_ac] = genes[u_ac] | np_map[np]
            ac_id_map[u_ac] = np_id_map[np]

            #replace the entries in the tag listed by refseqs with u-acs
            if np in tag_map:
                if not u_ac in tag_map:
                    tag_map[u_ac] = {}
                tag_map[u_ac].update(tag_map[np])
                del tag_map[np]

            if np in species_map:
                species_map[u_ac] = species_map[np]
                del species_map[np]


    #detect human datasets and add some entries to the species map
    human_set = True
    for u_ac in ac_id_map:
        u_id = ac_id_map[u_ac]
        if u_id.count('HUMAN') == 0:
            human_set = False
        if not u_ac in species_map:
            if u_id.count('_') == 1:
                species_map[u_ac] = u_id.split('_')[1]

    if human_set:
        species = 'human'

    stored_genes,new_genes = database.geneScan(genes,db,cursor)

    outlist = [stored_genes]
    #computation of auto_scale_template_percentage:
    s = len(new_genes)
    print "Total genes: ",s
    if s > 0:
        auto_scale_template_percentage = 1-0.5**(math.log10(s)-2)
        if auto_scale_template_percentage < 0.0:
            auto_scale_template_percentage = 0.0
        #print auto_scale_template_percentage
    
    #junksize is the maximal number of genes per cycle
    #return [genes],ac_id_map #if this not commented, all genes are taken in one cycle

    if s > junksize:
        if anno_anno_calculation:
            print "CAUTION. Anno-Anno calculation may be wrong for large datasets (>%s Genes)" % str(junksize)

        n_of_batches = s/junksize
        if s%junksize != 0:
            n_of_batches += 1
        batchsize = s/n_of_batches
        rest = s%n_of_batches
        outlist = []

        for i in range(0,n_of_batches):
            new_map = {}
            if i == n_of_batches - 1: # the last batch
                batchsize += rest
            for j in range(0,batchsize):
                (key,value) = new_genes.popitem()
                new_map[key] = value
            outlist.append(new_map)
        new_map = {}
        while len(new_genes) > 0:
            (key,value) = new_genes.popitem()
            new_map[key] = value
        if len(new_map) > 0:
            outlist.append(new_map)
    else:
        outlist.append(new_genes)
    return outlist,ac_id_map,u_ac_refseq_map,tag_map,species_map,fasta_map


def paraGetSequences(gene_queue,out_queue,lock,err_queue):
    with lock:
        gene_queue.put(None)
    while True:
        with lock:
            gene = gene_queue.get()
        if gene == None:
            return
        try:
            (wildtype_sequence,refseqs,go_terms,pathways) = uniprot.getSequence(gene)
            with lock:
                out_queue.put((gene,wildtype_sequence,refseqs,go_terms,pathways))
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            with lock:
                err_queue.put((e,f,g,gene))

def nToAA(seq):
    table={ 
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 
    'GGG': 'G', }
    stop_codons = ('TAG','TGA','TAA')

    i = 1
    triple = ''
    aa_seq = ''
    for char in list(seq):        
        triple += char
        if i%3 == 0:
            if triple in table:
                aa = table[triple]
            elif triple in stop_codons:
                aa = ''
            aa_seq += aa
            triple = ''
        i += 1
    return aa_seq

def getSequences(new_genes,stored_genes,fasta_map,species_map,u_ac_refseq_map,gene_mut_map_new,gene_mut_map_pseudo):
    global mrna_fasta
    global species
    global No_Errors
    global manager
    global lock
    global db
    global cursor
    global number_of_processes
    global human_id_mapping_path
    global human_proteome_path
    global pdb_path

    if species == 'human':
        gene_sequence_map,gene_info_map = uniprot.getHumanSequences(new_genes,u_ac_refseq_map,'%s/human_info_map.tab' % human_id_mapping_path,human_proteome_path)
        database.addGeneInfos(gene_info_map,db,cursor)
        gene_sequence_map = database.getSequences(stored_genes,gene_sequence_map,db,cursor)
        return gene_sequence_map,gene_mut_map_new,gene_mut_map_pseudo

    elif mrna_fasta == None:
        t0 = time.time()
        gene_queue = manager.Queue()
        pdb_dict = {}
        for gene in new_genes:
            if len(gene) == 6 and gene.count(':') == 1:
                pdb_dict[gene] = new_genes[gene]
            else:
                gene_queue.put(gene)
        stored_genes_lookup = {}
        for gene in stored_genes:
            if len(gene) == 6 and gene.count(':') == 1:
                pdb_dict[gene] = stored_genes[gene]
            else:
                stored_genes_lookup[gene] = stored_genes[gene]

        out_queue = manager.Queue()
        err_queue = manager.Queue()
        processes = {}
        for i in range(1,number_of_processes + 1):
            p = Process(target=paraGetSequences, args=(gene_queue,out_queue,lock,err_queue))
            processes[i] = p
            p.start()
        for i in processes:
            processes[i].join()

        t1 = time.time()
        err_queue.put(None)
        while True:
            err = err_queue.get()
            if err == None:
                break
            (e,f,g,gene) = err
            errortext = "getSequence Error: %s\n%s\n\n" % (gene,'\n'.join([str(e),str(f),str(g)]))
            f = open(errorlog,'a')
            f.write(errortext)
            f.close()
            No_Errors = False
        t2 = time.time()
        out_queue.put(None)
        gene_sequence_map = {}
        gene_info_map = {}
        while True:
            out = out_queue.get()
            if out == None:
                break
            (gene,sequence,refseqs,go_terms,pathways) = out
            gene_sequence_map[gene] = sequence
            gene_info_map[new_genes[gene]] = (refseqs,go_terms,pathways,sequence)
        del out_queue
        del gene_queue

        pdb_sequence_map,pdb_pos_map = pdb.getSequences(pdb_dict,pdb_path)
        gene_sequence_map.update(pdb_sequence_map)
        #print pdb_pos_map
        #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
        #structure of gene_mut_map_pseudo: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
        for u_id in gene_mut_map_new:
            if not u_id in pdb_pos_map:
                continue
            (gene_id,aac_map) = gene_mut_map_new[u_id]
            new_aac_map = {}
            for aac_base in aac_map:
                res_id = aac_base[1:]
                if not res_id in pdb_pos_map[u_id]:
                    continue
                pos = pdb_pos_map[u_id][res_id]
                new_aac_base = '%s%s' % (aac_base[0],str(pos))
                new_aac_map[new_aac_base] = aac_map[aac_base]
            gene_mut_map_new[u_id] = (gene_id,new_aac_map)
            #structure of gene_info_map: {gene_id:(refseqs,{go_term_id:go_term_name},{reactome_id:pathway_name}),sequence}
            gene_info_map[gene_id] = (u_id,{},{},pdb_sequence_map[u_id])

        for u_id in gene_mut_map_pseudo:
            if not u_id in pdb_pos_map:
                continue
            (gene_id,aac_map) = gene_mut_map_pseudo[u_id]
            new_aac_map = {}
            for aac_base in aac_map:
                res_id = aac_base[1:]
                pos = pdb_pos_map[u_id][res_id]
                new_aac_base = '%s%s' % (aac_base[0],str(pos))
                new_aac_map[new_aac_base] = aac_map[aac_base]
            gene_mut_map_pseudo[u_id] = (gene_id,new_aac_map)

        t3 = time.time()
        database.addGeneInfos(gene_info_map,db,cursor)
        t4 = time.time()
        gene_sequence_map = database.getSequences(stored_genes_lookup,gene_sequence_map,db,cursor)
        t5 = time.time()

        print "Time for getSequences Part 1: %s" % str(t1-t0)
        print "Time for getSequences Part 2: %s" % str(t2-t1)
        print "Time for getSequences Part 3: %s" % str(t3-t2)
        print "Time for getSequences Part 4: %s" % str(t4-t3)
        print "Time for getSequences Part 5: %s" % str(t5-t4)

        return gene_sequence_map,gene_mut_map_new,gene_mut_map_pseudo
    else:
        gene_sequence_map = {}
        for species in fasta_map:
            mrna_fasta = fasta_map[species]
            f = open(mrna_fasta,'r')
            lines = f.readlines()
            f.close()

            gene_seq_map = {}
            g_gene_map = {}
            for line in lines:
                line = line.replace('\n','')
                #print line
                if line[0] == '>':
                    g_id = line.split()[0][1:]
                    gene_name = line.split()[1]
                    if species == "multi":
                        gene_name = g_id
                    gene_seq_map[gene_name] = ''
                    g_gene_map[g_id] = gene_name
                else:
                    gene_seq_map[gene_name] += line

            gene_info_map = {}

            for gene in new_genes:
                if species_map[gene] == species:
                    if species != "multi":
                        gene_name = gene.replace('_%s' % species,'')
                    else:
                        gene_name = gene
                    if not gene_name in gene_seq_map:
                        raise NameError("Did not found gene: %s in %s" % (gene_name,mrna_fasta))
                    sequence = nToAA(gene_seq_map[gene_name])
                    gene_sequence_map[gene] = sequence
                    gene_info_map[new_genes[gene]] = ([],{},{},sequence)
            database.addGeneInfos(gene_info_map,db,cursor)
            gene_sequence_map.update(database.getSequences(stored_genes,gene_sequence_map,db,cursor))
    return gene_sequence_map,gene_mut_map_new,gene_mut_map_pseudo

def paraBlast(process_queue_blast,out_queue,lock,i,err_queue):
    global blast_path
    global blast_db_path
    global option_number_of_templates
    global option_seq_thresh
    global option_ral_thresh
    global option_res_thresh
    global pdb_path
    global cwd
    with lock:
        process_queue_blast.put(None)
    while True:
        with lock:
            inp = process_queue_blast.get()
        if inp == None:
            return
        try:
            (gene,seq) = inp
            target_name = gene.replace("/","")
            #print  gene
            templates,oligo_map = blast.blast(seq,target_name,blast_path,blast_db_path,option_number_of_templates,option_seq_thresh,option_ral_thresh,cwd = "%s/%d" %(cwd,i))
            #print templates
            templates,del_templates = templateSelection.selectTemplates(templates,pdb_path)
            #print templates
            templates = templateFiltering.filterTemplates(templates,option_seq_thresh,option_res_thresh,option_ral_thresh)
            #print templates
            with lock:
                out_queue.put((gene,templates,oligo_map))
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            with lock:
                err_queue.put((e,f,g,gene))

def autoTemplateSelection(gene_sequence_map,stored_genes,new_genes,gene_mut_map_pseudo,AS_db):
    #structure of stored_genes: {Uniprot-Id:(gene_id,more_restrictive)}
    #structure of new_genes: {Uniprot-Id:gene_id}
    global db
    global cursor
    global manager
    global lock
    global blast_processes
    global option_seq_thresh
    global option_ral_thresh
    global option_res_thresh
    global session
    global cwd
    global option_number_of_templates
    global pdb_path
    global blast_path
    global blast_db_path
    global No_Errors
    global anno_session_mapping
    global anno_anno_calculation

    t0 = time.time()

    process_list_db = set([])
    process_queue_blast = manager.Queue()
    blast_queues = [process_queue_blast]
    n = 0
    gene_error_map = {}

    #seq_map = {}

    for gene in new_genes:
        seq = gene_sequence_map[gene]
        if seq == 0 or seq == 1:
            gene_error_map[new_genes[gene]] = seq
        elif seq == "":
            gene_error_map[new_genes[gene]] = 2
        else:
            process_queue_blast.put((gene,seq))
            n += 1
            if n == 1000:
                process_queue_blast = manager.Queue()
                blast_queues.append(process_queue_blast)
                n = 0
            #seq_map[gene] = seq

    gene_id_gene_map = {}
    for gene in stored_genes:
        gene_id_gene_map[stored_genes[gene][0]] = gene
        if stored_genes[gene][1]:
            process_list_db.add(stored_genes[gene][0])
        else:
            process_list_db.add(stored_genes[gene][0])
            seq = gene_sequence_map[gene]
            process_queue_blast.put((gene,seq))
            n += 1
            if n == 1000:
                process_queue_blast = manager.Queue()
                blast_queues.append(process_queue_blast)
                n = 0
            #seq_map[gene] = seq

    gene_template_map = database.getTemplates(process_list_db,db,cursor)
    #print gene_template_map
    for gene_id in gene_template_map:
        templates = []
        for template_id in gene_template_map[gene_id][0]:
            template = gene_template_map[gene_id][0][template_id]
            templates.append(template)
        templates = templateFiltering.filterTemplates(templates,option_seq_thresh,option_res_thresh,option_ral_thresh)
        del_list = []
        for template_id in gene_template_map[gene_id][0]:
            pdb_id = gene_template_map[gene_id][0][template_id][0]
            filtered = True
            for template in templates:
                if pdb_id == template[0]:
                    gene_template_map[gene_id][0][template_id] = template
                    filtered = False
            if filtered:
                del_list.append(template_id)
        for template_id in del_list:
            del gene_template_map[gene_id][0][template_id]
 
    t1 = time.time()

    #"""paralized blasting of single sequences  
    #print "before parablast"
    t12 = 0.0
    t3 = 0.0
    t2 = 0.0
    t4 = 0.0
    for process_queue_blast in blast_queues:
        t12 += time.time()
        out_queue = manager.Queue()
        err_queue = manager.Queue()
        processes = {}
        for i in range(1,blast_processes + 1):
            try:
                os.stat("%s/%d" %(cwd,i))
            except:
                os.mkdir("%s/%d" %(cwd,i))
            p = Process(target=paraBlast, args=(process_queue_blast,out_queue,lock,i,err_queue))
            processes[i] = p
            p.start()
        for i in processes:
            processes[i].join()

        err_queue.put(None)
        while True:
            err = err_queue.get()
            if err == None:
                break
            (e,f,g,gene) = err
            errortext = "BLAST Error: %s\n%s\n\n" % (gene,'\n'.join([str(e),str(f),str(g)]))
            f = open(errorlog,'a')
            f.write(errortext)
            f.close()
            No_Errors = False

        t2 += time.time()
        #print "after parablast"

        pseudo_template_map = {}
        out_queue.put(None)
        while True:
            out = out_queue.get()
            if out == None:
                break
            (gene,templates,oligo_map) = out
            if gene in stored_genes:
                gene_id = stored_genes[gene][0]
                if gene_id not in gene_template_map:
                    pseudo_template_map[gene] = (templates,[],oligo_map)
                else:
                    pseudo_template_map[gene] = (templates,gene_template_map[gene_id][0],oligo_map.update(gene_template_map[gene_id][1]))
            else:
                pseudo_template_map[gene] = (templates,[],oligo_map)
        #"""
        del out_queue
        del process_queue_blast

        for gene_id in gene_template_map:
            gene = gene_id_gene_map[gene_id]
            if gene not in pseudo_template_map:
                pseudo_template_map[gene] = ([],gene_template_map[gene_id][0],gene_template_map[gene_id][1])
        
        del_list = []
        for gene in pseudo_template_map:
            if len(pseudo_template_map[gene][0]) == 0 and len(pseudo_template_map[gene][1]) == 0:
                if gene in new_genes:
                    gene_error_map[new_genes[gene]] = 3
                del_list.append(gene)
        for gene in del_list:
            del pseudo_template_map[gene]

        t3 += time.time()

        #pseudo mutation/pseudo template combinations can be directly added to results
        #this step is the acutal time saver of the pseudo store system
        #additional do it as a separate process in the background
        stored_annotations,background_process = database.insertMultipleAnnotationSession(session,pseudo_template_map,gene_mut_map_pseudo,db,cursor,AS_db,anno_session_mapping = anno_session_mapping,anno_anno_calculation=anno_anno_calculation)
        #print "after insertmultiple"
        t4 += time.time()

    print "Template Selection Part 1: %s" % (str(t1-t0))
    print "Template Selection Part 2: %s" % (str(t2-t12))
    print "Template Selection Part 3: %s" % (str(t3-t2))
    print "Template Selection Part 4: %s" % (str(t4-t3))

    return pseudo_template_map,gene_error_map,stored_annotations,background_process

def paraAlignment(gene_sequence_map,pseudo_template_map,stored_genes,new_genes,gene_mut_map_new,gene_mut_map_pseudo,stored_annotations):
    #structure of pseudo_template_map: {Uniprot-Id:(template-list,{template-id:stored-template},oligo_map)}
    #structure of gene_sequence_map: {Uniprot_Id:Sequence}
    #structure of stored_genes: {Uniprot-Id:(gene_id,more_restrictive)}
    #structure of new_genes: {Uniprot-Id:gene_id}
    #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #structure of gene_mut_map_pseudo: {Uniprot_Id:(gene_id,{AAC_base:mutation_id})}
    #structure of stored_annotations: {Template_Id:{AAC_Base:Mutation_Id}}
    global manager
    global lock
    global alignment_processes
    global cwd
    global db
    global cursor
    global No_Errors
    global pdb_path
    global smiles_path
    global inchi_path
    global session
    t0 = time.time()
    input_queue = manager.Queue()
    err_queue = manager.Queue()

    in_queues = [input_queue] #partioning the parallized alignment part shall reduce the memory usage, by safing the alignments to the database every 1000 genes.

    stored_map = {}
    gene_id_gene_map = {}
    template_id_template_map = {}

    relax_map = {}
    new_structures = set()

    n = 0

    #temp_amount = 0

    for gene in pseudo_template_map:
        if gene in stored_genes:
            gene_id = stored_genes[gene][0]
        else:
            gene_id = new_genes[gene]
        gene_id_gene_map[gene_id] = gene
        (templates,pseudo_templates,oligo_map) = pseudo_template_map[gene]

        #temp_amount += len(templates)

        if gene in gene_mut_map_new:
            if gene in gene_mut_map_pseudo:
                aaclist = gene_mut_map_new[gene][1].keys() + gene_mut_map_pseudo[gene][1].keys()
            else:
                aaclist = gene_mut_map_new[gene][1].keys()
            aaclist_new = gene_mut_map_new[gene][1].keys()
        else:
            aaclist = gene_mut_map_pseudo[gene][1].keys()
            aaclist_new = []
        for template in templates:
            

            input_queue.put((gene,gene_sequence_map[gene],template,aaclist,oligo_map[template[0]]))
            new_structures.add(template[0])
        if n == 1000:
            input_queue = manager.Queue()
            in_queues.append(input_queue)
            n = 0
        n += 1

        for template_id in pseudo_templates:
            pdb_id = pseudo_templates[template_id][0]
            if len(aaclist_new) > 0 or template_id in stored_annotations:
                
                if template_id in stored_annotations:
                    aaclist_new = aaclist #The concept of saving the stored annotations in order to compute all anno-anno distances may lead to uselessness of the pseudo-template concept, potential removal of the pseudo-template system does not yield performance improvements and is a though quest. => without further reasons to remove the pseudo-template system, leave it be.
                    if pdb_id not in relax_map:
                        relax_map[pdb_id] = set()
                    relax_map[pdb_id].add((gene_id,template_id))
                else:
                    new_structures.add(pdb_id)
                template_id_template_map[template_id] = pdb_id
                if not gene_id in stored_map:
                    stored_map[gene_id] = {template_id:aaclist_new}
                else:
                    stored_map[gene_id][template_id] = aaclist_new
            else:
                relax_map[pdb_id].add((gene_id,template_id))

    print "Size of relax_map: ",len(relax_map)

    for pdb_id in relax_map:
        if pdb_id in new_structures:
            continue
        for (gene_id,template_id) in relax_map[pdb_id]:
            del stored_map[gene_id][template_id] #The stored datastructures can be relaxed, by the structures, which are completely coverd in the stored datastructres
            del stored_annotations[template_id]

    t1 = time.time()
    #structure of gene_template_alignment_map: {gene_id:{template_id:(aln_length,seq_id,sub_infos)}}
    pdb_map = database.getAlignments(stored_map,db,cursor)
    t15 = time.time()
    input_queue = manager.Queue()
    for pdb_id in pdb_map:
        input_queue.put(pdb_id)

    out_queue = manager.Queue()
    processes = {}
    for i in range(1,alignment_processes + 1):
        try:
            os.stat("%s/%d" %(cwd,i))
        except:
            os.mkdir("%s/%d" %(cwd,i))
        p = Process(target=createGTAmap, args=(input_queue,out_queue,pdb_map,lock,err_queue))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join() 

    err_queue.put(None)
    while True:
        err = err_queue.get()
        if err == None:
            break
        (e,f,g,pdb_id) = err
        errortext = "GTAmap Error: %s\n%s\n\n" % (pdb_id,'\n'.join([str(e),str(f),str(g)]))
        f = open(errorlog,'a')
        f.write(errortext)
        f.close()
        No_Errors = False

    out_queue.put(None)
    gene_id_template_id_alignment_map = {}
    
    while True:
        out = out_queue.get()
        if out == None:
            break
        (gene_id,template_id,aln_length,seq_id,sub_infos) = out

        if not gene_id in gene_id_template_id_alignment_map:
            gene_id_template_id_alignment_map[gene_id] = {}
        gene_id_template_id_alignment_map[gene_id][template_id] = (aln_length,seq_id,sub_infos)



    gene_template_alignment_map = {}
    for gene_id in gene_id_template_id_alignment_map:
        gene = gene_id_gene_map[gene_id]
        gene_template_alignment_map[gene] = {}
        for template_id in gene_id_template_id_alignment_map[gene_id]:
            pdb_id = template_id_template_map[template_id]
            gene_template_alignment_map[gene][pdb_id] = gene_id_template_id_alignment_map[gene_id][template_id] #Is this a horrible bug??? Because there are different templates with identical pdb_id and/or different chains - But also with the same gene_id ??? => NO! (currently,) each gene can mapped only on one chain per structure (the one with the longest coverage -> see blast.py)

    #print "before alignment"

    template_sub_amount = 0
    t2 = time.time()
    t21 = 0.0
    t3 = 0.0
    t4 = 0.0
    t5 = 0.0
    update_gene_template_map = {}
    template_id_map = {}
    for input_queue in in_queues:
        t21 += time.time()
        out_queue = manager.Queue()
        error_queue = manager.Queue()
        err_queue = manager.Queue()
        processes = {}
        for i in range(1,alignment_processes + 1):
            try:
                os.stat("%s/%d" %(cwd,i))
            except:
                os.mkdir("%s/%d" %(cwd,i))
            p = Process(target=align, args=(input_queue,out_queue,error_queue,lock,cwd,i,err_queue))
            processes[i] = p
            p.start()
        for i in processes:
            processes[i].join()
            
        err_queue.put(None)
        multiple_errors = set()
        while True:
            err = err_queue.get()
            if err == None:
                break
            (e,f,g,gene,pdb_id,chain) = err
            if (e,gene) in multiple_errors:
                continue
            multiple_errors.add((e,gene))
            errortext = "Alignment Error: %s, %s:%s\n%s\n\n" % (gene,pdb_id,chain,'\n'.join([str(e)]))
            f = open(errorlog,'a')
            f.write(errortext)
            f.close()
            No_Errors = False

        #print "after alignment"
        t3 += time.time()

        error_map = {}
        error_queue.put(None)
        while True:
            err = error_queue.get()
            if err == None:
                break
            (gene,pdb_id,error) = err
            if not gene in error_map:
                error_map[gene] = {}
            error_map[gene][pdb_id] = error

        for gene in error_map:
            for pdb_id in error_map[gene]:
                print error_map[gene][pdb_id]        

        out_queue.put(None)
        template_insertion_list = []

        while True:
            out = out_queue.get()
            if out == None:
                break
            (gene,pdb_id,template,aln_length,seq_id,sub_infos,alignment_pir,new_oligos) = out
           
            if gene in stored_genes:
                gene_id = stored_genes[gene][0]
            else:
                gene_id = new_genes[gene]
            
            if not gene in update_gene_template_map:
                update_gene_template_map[gene] = []
            update_gene_template_map[gene].append(template)

            if not gene in gene_template_alignment_map:
                gene_template_alignment_map[gene] = {}
            gene_template_alignment_map[gene][pdb_id] = (aln_length,seq_id,sub_infos)
            template_sub_amount += len(sub_infos)
            #print pdb_id
            #print pseudo_template_map[gene][2]

            #structure of pseudo_template_map: {Uniprot-Id:(template-list,{template-id:stored-template},oligo_map)}
            
            #If in the alignment process, a structure had to be replaced by asymetric unit, the identifier in the oligo map has to be updated
            if pdb_id.count('_na') > 0:
                old_pdb_id = pdb_id.replace('_na','')
                if old_pdb_id in pseudo_template_map[gene][2]:
                    pseudo_template_map[gene][2][pdb_id] = new_oligos
                    del pseudo_template_map[gene][2][old_pdb_id] # delete the old entry
            #If not, the entry itself has still to be updated
            else:
                pseudo_template_map[gene][2][pdb_id] = new_oligos

            template_insertion_list.append((gene_id,template,pseudo_template_map[gene][2][pdb_id],alignment_pir))

        del input_queue
        del out_queue
        del error_queue

        t4 += time.time()
        #structure of template_id_map[gene_id] = {pdb_id:template_id}
        new_template_id_map = database.insertTemplates(template_insertion_list,session,db,cursor,cwd,pdb_path,smiles_path,inchi_path)
        for gene_id in new_template_id_map:
            if not gene_id in template_id_map:
                template_id_map[gene_id] = new_template_id_map[gene_id]
            else:
                for pdb_id in new_template_id_map[gene_id]:
                    template_id_map[gene_id][pdb_id] = new_template_id_map[gene_id][pdb_id]
        t5 += time.time()            

    print "Amount of template-mutation pairs after alignment: ",template_sub_amount
        
    for gene in update_gene_template_map:
        pseudo_template_map[gene] = (update_gene_template_map[gene],pseudo_template_map[gene][1],pseudo_template_map[gene][2])

    #after_amount = 0

    for gene in pseudo_template_map:
        templates = pseudo_template_map[gene][0]
        if gene in stored_genes:
            gene_id = stored_genes[gene][0]
        else:
            gene_id = new_genes[gene]
        new_template_map = {}
        if gene_id in template_id_map:
            for template in templates:
                new_template_map[template_id_map[gene_id][template[0]]] = template
        pseudo_template_map[gene] = (new_template_map,pseudo_template_map[gene][1],pseudo_template_map[gene][2])
        #after_amount += len(templates)

    #print "Template amounts, before and after: ",temp_amount,after_amount


    print "Alignment Part 1: %s" % (str(t1-t0))
    print "Alignment Part 1.5: %s" % (str(t15-t1))
    print "Alignment Part 2: %s" % (str(t2-t15))
    print "Alignment Part 3: %s" % (str(t3-t21))
    print "Alignment Part 4: %s" % (str(t4-t3))
    print "Alignment Part 5: %s" % (str(t5-t4))

    return pseudo_template_map,gene_template_alignment_map

def createGTAmap(input_queue,out_queue,pdb_map,lock,err_queue):
    with lock:
        input_queue.put(None)
    while True:
        with lock:
            pdb_id = input_queue.get()
        if pdb_id == None:
            return
        try:
            template_page = pdb.standardParsePDB(pdb_id,pdb_path)
            for chain in pdb_map[pdb_id]:
                seq_res_map = globalAlignment.createTemplateFasta(template_page,pdb_id,chain,onlySeqResMap = True)
                for (target_seq,template_seq,aaclist,gene_id,aln_length,seq_id,template_id) in pdb_map[pdb_id][chain]:

                    sub_infos,errors = globalAlignment.getSubPos(target_seq,template_seq,aaclist,seq_res_map)
                    if len(sub_infos) == 0:
                        raise NameError("There are no sub_infos for template: %s" % str(template_id))

                    out_queue.put((gene_id,template_id,aln_length,seq_id,sub_infos))
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            with lock:
                err_queue.put((e,f,g,pdb_id))

def align(input_queue,out_queue,error_queue,lock,cwd,proc_number,err_queue):
    global pdb_path
    t0 = 0.0
    t1 = 0.0
    t2 = 0.0
    aligntimes = [0.0]*5
    parsetimes = [0.0]*3
    with lock:
        input_queue.put(None)
    while True:
        with lock:
            inp = input_queue.get()
        if inp == None:
            #print "Align part1: ",t1-t0
            #print "Detailed part 1 times: "
            #print parsetimes
            #print "Align part2: ",t2-t1
            #print "Detailed part 2 times: "
            #print aligntimes
            return
        (gene,wildtype_sequence,template,aaclist,oligos) = inp
        [pdb_id,seq_id,chain,aln_length] = template[0:4]

        #print "before get pdb file"
        try:
            t0 += time.time()

            (template_page,template,new_oligos,error,part1_times) = pdb.getStandardizedPdbFile(template,pdb_path,oligos=oligos)
            n = 0
            for ti in part1_times:
                parsetimes[n] += ti
                n+=1

            t1 += time.time()
            if error == None:
                """
                try:
                    os.stat("%s/%d" %(cwd,proc_number))
                except:
                    os.mkdir("%s/%d" %(cwd,proc_number))  
                os.chdir("%s/%d" %(cwd,proc_number))
                temp_cwd = "%s/%d" %(cwd,proc_number)
                """
                #print "before alignment"
                #print wildtype_sequence,aaclist
                (aln_length,seq_id,sub_infos,alignment_pir,errors,times) = globalAlignment.alignBioPython(gene,wildtype_sequence,template[0],template_page,template[2],aaclist)
                n = 0
                for ti in times:
                    aligntimes[n] += ti
                    n+=1

                for err in errors:
                    [e,f,g] = sys.exc_info()
                    g = traceback.format_exc(g)
                    with lock:
                        err_queue.put((err,f,g,gene,pdb_id,chain))
                #print "after alignment"

                with lock:
                    out_queue.put((gene,template[0],template,aln_length,seq_id,sub_infos,alignment_pir,new_oligos))
                #os.chdir(cwd)
            else:
                with lock:
                    error_queue.put((gene,pdb_id,error))
            t2 += time.time()
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            with lock:
                err_queue.put((e,f,g,gene,pdb_id,chain))
    

def paraAnnotate(stored_genes,new_genes,gene_mut_map_new,gene_mut_map_pseudo,gene_template_alignment_map,pseudo_template_map,stored_annotations):
    global annotation_processes
    global manager
    global lock
    global db
    global cursor
    global option_lig_wf
    global option_chain_wf
    global cwd
    global session
    global pdb_path
    global full_anno_anno_calculation
    global No_Errors
    global anno_session_mapping
    global error_annotations_into_db
    #new structure of pseudo_template_map: {Uniprot-Id:({template-id:new template},{template-id:stored-template},oligo_map)}
    #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #structure of gene_mut_map_pseudo: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #structure of gene_template_alignment_map: {Uniprot_Id:{PDB_Id:(aln_length,seq_id,sub_infos)}}
    t0 = time.time()
    
    s_amount = 0
    t_amount = 0

    template_mutation_map = {}
    salvation_set = set()
    for gene in pseudo_template_map:
        if gene in gene_template_alignment_map:  
            (new_templates,pseudo_templates,oligo_map) = pseudo_template_map[gene]
            t_amount += len(new_templates)
            aac_map_new = {}
            aac_map_combined = {}
            if gene in gene_mut_map_new:
                for aac_base in gene_mut_map_new[gene][1]:
                    aac_map_new[aac_base] = gene_mut_map_new[gene][1][aac_base]
                    aac_map_combined[aac_base] = gene_mut_map_new[gene][1][aac_base]
            if gene in gene_mut_map_pseudo:
                for aac_base in gene_mut_map_pseudo[gene][1]:
                    aac_map_combined[aac_base] = gene_mut_map_pseudo[gene][1][aac_base]
            for template_id in new_templates:
                template = new_templates[template_id]
                pdb_id = template[0]
                if pdb_id in gene_template_alignment_map[gene]:
                    sub_infos = gene_template_alignment_map[gene][pdb_id][2]
                    s_amount += len(sub_infos)
                    if not pdb_id in template_mutation_map:
                        template_mutation_map[pdb_id] = [{template_id:(sub_infos,template[2],oligo_map[pdb_id])},{template_id:aac_map_combined},template]
                        salvation_set.add(pdb_id)
                    else:
                        template_mutation_map[pdb_id][0][template_id] = (sub_infos,template[2],oligo_map[pdb_id]) 
                        template_mutation_map[pdb_id][1][template_id] = aac_map_combined
                        """
                        for aac in aac_map_combined:
                            if not aac in template_mutation_map[pdb_id][1]:
                                template_mutation_map[pdb_id][1][aac] = aac_map_combined[aac]
                            else:
                                template_mutation_map[pdb_id][1][aac].append(aac_map_combined[aac][0])
                        """
            
            for template_id in pseudo_templates:
                template = pseudo_templates[template_id]
                pdb_id = template[0]
                if pdb_id in gene_template_alignment_map[gene]:
                    sub_infos = gene_template_alignment_map[gene][pdb_id][2]
                    if template_id not in stored_annotations:
                        salvation_set.add(pdb_id)
                    if not pdb_id in template_mutation_map:
                        olis = oligo_map[pdb_id]
                        template_mutation_map[pdb_id] = [{template_id:(sub_infos,template[2],olis)},{template_id:aac_map_new},template]
                    else:
                        olis = oligo_map[pdb_id]
                        template_mutation_map[pdb_id][0][template_id] = (sub_infos,template[2],olis)
                        template_mutation_map[pdb_id][1][template_id] = aac_map_new
                        """
                        for aac in aac_map_new:
                            if not aac in template_mutation_map[pdb_id][1]:
                                template_mutation_map[pdb_id][1][aac] = aac_map_new[aac]
                            else:
                                template_mutation_map[pdb_id][1][aac].append(aac_map_new[aac][0])
                         """
    print "new templates used for template_mutation_map: ",t_amount
    print "template mutation pairs going into template_mutation_map: ",s_amount

    #The deletion set removes all structures, whose templates are all covered in stored_annotations
    deletion_set = []
    for pdb_id in template_mutation_map:
        if not pdb_id in salvation_set:
            deletion_set.append(pdb_id)                
    
    print "Deletion-set-size: %s" % (str(len(deletion_set)))

    for pdb_id in deletion_set:
        if pdb_id in template_mutation_map:
            del template_mutation_map[pdb_id]

    
    input_queue = manager.Queue()
    #input_queues = [input_queue]
    size = 200 #Maximal amount of structures going into the paralellized annotation. The larger, the faster the pipeline, but consumes more memory.
    n = 0
    total_size = len(template_mutation_map)

    t1 = 0.0
    t2 = 0.0
    t3 = 0.0
    t4 = 0.0
    t12 = time.time()
    for pdb_id in template_mutation_map:
        input_queue.put(template_mutation_map[pdb_id])
        n += 1
        if n == size or n == total_size:
            t1 += time.time()
            n = 0
            total_size -= size
            total_annotations = {}
            min_anno_anno_dict = {}
            out_queue = manager.Queue()
            error_queue = manager.Queue()
            err_queue = manager.Queue()
            processes = {}
            for i in range(1,annotation_processes + 1):
                try:
                    os.stat("%s/%d" %(cwd,i))
                except:
                    os.mkdir("%s/%d" %(cwd,i))
                p = Process(target=annotate, args=(input_queue,out_queue,error_queue,lock,cwd,i,stored_annotations,err_queue))
                processes[i] = p
                p.start()
            for i in processes:
                processes[i].join()

            err_queue.put(None)
            while True:
                err = err_queue.get()
                if err == None:
                    break
                (e,f,g,pdb_id) = err
                errortext = "Annotation Error: %s\n%s\n\n" % (pdb_id,'\n'.join([str(e),str(f),str(g)]))
                f = open(errorlog,'a')
                f.write(errortext)
                f.close()
                No_Errors = False

            t2 += time.time()

            
            with lock:
                out_queue.put(None)

            while True:
                out = out_queue.get()
                if out == None:
                    break
                (annotations,anno_anno_dict) = out
                #print len(annotations)
                total_annotations.update(annotations)
                #save only the minimal distance between two mutations
                if not full_anno_anno_calculation:
                    for (template_id,template_id_2,aac_base,aac_base_2) in anno_anno_dict:
                        (min_d,(atom,atom2),chain,chain_2) = anno_anno_dict[(template_id,template_id_2,aac_base,aac_base_2)]
                        pdb_id = annotations[template_id][0]
                        pdb_id_2 = annotations[template_id_2][0]
                        if not aac_base in template_mutation_map[pdb_id][1][template_id]:
                            continue
                        if not aac_base_2 in template_mutation_map[pdb_id_2][1][template_id_2]:
                            continue
                        mut_id = template_mutation_map[pdb_id][1][template_id][aac_base]
                        mut_id_2 = template_mutation_map[pdb_id_2][1][template_id_2][aac_base_2]
                        if int(mut_id) < int(mut_id_2):
                            if not (mut_id,mut_id_2) in min_anno_anno_dict:
                                min_anno_anno_dict[(mut_id,mut_id_2)] = (template_id,template_id_2,min_d,atom,atom2,chain,chain_2)
                            elif min_d < min_anno_anno_dict[(mut_id,mut_id_2)][2]:
                                min_anno_anno_dict[(mut_id,mut_id_2)] = (template_id,template_id_2,min_d,atom,atom2,chain,chain_2)
                        else:
                            if not (mut_id_2,mut_id) in min_anno_anno_dict:
                                min_anno_anno_dict[(mut_id_2,mut_id)] = (template_id_2,template_id,min_d,atom,atom2,chain_2,chain)
                            elif min_d < min_anno_anno_dict[(mut_id_2,mut_id)][2]:
                                min_anno_anno_dict[(mut_id_2,mut_id)] = (template_id_2,template_id,min_d,atom,atom2,chain_2,chain)
                #or save all distances found in all strucutures
                else:
                    for (template_id,template_id_2,aac_base,aac_base_2) in anno_anno_dict:
                        (min_d,(atom,atom2),chain,chain_2) = anno_anno_dict[(template_id,template_id_2,aac_base,aac_base_2)]
                        pdb_id = annotations[template_id][0]
                        pdb_id_2 = annotations[template_id_2][0]
                        if not aac_base in template_mutation_map[pdb_id][1][template_id]:
                            continue
                        if not aac_base_2 in template_mutation_map[pdb_id_2][1][template_id_2]:
                            continue
                        mut_id = template_mutation_map[pdb_id][1][template_id][aac_base]
                        mut_id_2 = template_mutation_map[pdb_id_2][1][template_id_2][aac_base_2]
                        if int(mut_id) < int(mut_id_2):
                            if not (mut_id,mut_id_2) in min_anno_anno_dict:
                                min_anno_anno_dict[(mut_id,mut_id_2)] = [(template_id,template_id_2,min_d,atom,atom2,chain,chain_2)]
                            else:
                                min_anno_anno_dict[(mut_id,mut_id_2)].append((template_id,template_id_2,min_d,atom,atom2,chain,chain_2))
                        else:
                            if not (mut_id_2,mut_id) in min_anno_anno_dict:
                                min_anno_anno_dict[(mut_id_2,mut_id)] = [(template_id_2,template_id,min_d,atom,atom2,chain_2,chain)]
                            else:
                                min_anno_anno_dict[(mut_id_2,mut_id)].append((template_id_2,template_id,min_d,atom,atom2,chain_2,chain))

            t3  += time.time()

            database.insertAnnotations(total_annotations,min_anno_anno_dict,template_mutation_map,option_lig_wf,option_chain_wf,session,db,cursor,full_anno_anno_calculation=full_anno_anno_calculation, error_annotations_into_db=error_annotations_into_db,anno_session_mapping=anno_session_mapping)

            t4  += time.time()
            input_queue = manager.Queue()

    #the partioned annotation process leads to multiple entries in the anno_anno table in the db, if these numbers tend to explode, a filtering done here can be the solution. Or try measuring the space for the min_anno_anno_dict outside the loop.

    print "Annotation Part 1: %s" % (str(t12-t0))
    print "Annotation Part 2: %s" % (str(t2-t1))
    print "Annotation Part 3: %s" % (str(t3-t2))
    print "Annotation Part 4: %s" % (str(t4-t3))

def annotate(input_queue,out_queue,error_queue,lock,cwd,proc_number,stored_annotations,err_queue):
    global compute_surface
    global anno_anno_calculation
    global calculate_interaction_profiles
    global dssp
    global dssp_path
    global pdb_path
    global rin_db_path

    with lock:
        input_queue.put(None)
    while True:
        inp = input_queue.get()
        if inp == None:
            return
        #print inp
        [sub_infos,aac_map_new,template] = inp
        try:
            annotations,anno_anno_dict,errorlist = templateFiltering.structuralAnalysis(sub_infos,template[0],stored_annotations,pdb_path,dssp_path,rin_db_path,anno_anno=anno_anno_calculation,dssp=dssp,calculate_interaction_profiles=calculate_interaction_profiles)
            with lock:
                out_queue.put((annotations,anno_anno_dict))
                if len(errorlist) > 0:
                    for (error,e,f,g) in errorlist:
                        err_queue.put((error,f,g,template[0]))
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            with lock:
                err_queue.put((e,f,g,template[0]))

def main(filename,config,output_path,main_file_path):
    global session
    global db
    global cursor
    global option_seq_thresh
    global option_res_thresh
    global option_ral_thresh
    global cwd
    global species

    setBasePaths(main_file_path)
    setPaths(output_path,config)

    parseConfig(config)

    t0 = time.time()

    if mrna_fasta != None and species == None:
        species = mrna_fasta.split('/')[-1].replace('.fa','').replace('.fasta','')
    if mrna_fasta != None:
        if not os.path.exists(mrna_fasta):
            raise NameError("mRNA path not found: %s" % mrna_fasta)

    #annovar-pipeline in case of vcf-file
    if filename.rsplit(".",1)[1] == "vcf":
        anno_db = "%s_annovar" % db_name.rsplit("_",1)[0]
        #print filename
        nfname = annovar.annovar_pipeline(filename,tax_id,annovar_path,db_adress,db_user_name,db_password,anno_db,mrna_fasta,ref_id=ref_genome_id)
        #print nfname
    else:
        nfname = filename

    #print db_adress,db_user_name,db_password,db_name

    db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    GS_db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    MS_db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    AS_db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    cursor = db.cursor()

    junksize = 500
    gene_aaclist_map_list,ac_id_map,u_ac_refseq_map,tag_map,species_map,fasta_map = buildQueue(nfname,junksize,mrna_fasta=mrna_fasta)
    print "Number of junks: ",len(gene_aaclist_map_list)

    newsession = False
    if session == 0:     
        starttime = SQLDateTime()
        session = database.insertSession(starttime,nfname,option_seq_thresh,option_res_thresh,option_ral_thresh,db,cursor)
        newsession = True

    date = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
    errortext  = "###############################################################################\n%s:\n%s - Session: %s\n" % (date,nfname,str(session))
    No_Errors = True

    f = open(errorlog,'a')
    f.write(errortext)
    f.close()



    errors = []
    print errorlog
    for gene_aaclist_map in gene_aaclist_map_list:

        if len(gene_aaclist_map) == 0:
            #print "empty gene_aaclist_map"
            continue

        try:
            os.stat("%s/tmp_structman_pipeline" %(output_path))
        except:
            os.mkdir("%s/tmp_structman_pipeline" %(output_path))  
        os.chdir("%s/tmp_structman_pipeline" %(output_path))
        cwd = "%s/tmp_structman_pipeline" %(output_path)

        try:
            print "Before geneCheck"
            #check for already stored genes and make all the necessary database interactions
            #structure of stored_genes: {Uniprot-Id:(gene_id,more_restrictive)}
            #structure of new_genes: {Uniprot-Id:gene_id}
            stored_genes,new_genes,stored_gene_ids,new_gene_ids,background_process_GS = database.geneCheck(gene_aaclist_map,ac_id_map,species_map,session,db,cursor,option_seq_thresh,option_res_thresh,option_ral_thresh,GS_db)
            #print stored_genes

            print "Before mutationCheck", len(stored_genes), len(new_genes)
            #check for already stored mutations or position twins and make all the necessary database interactions
            #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
            #structure of gene_mut_map_pseudo: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
            gene_mut_map_new,gene_mut_map_pseudo,background_process_MS = database.mutationCheck(gene_aaclist_map,stored_genes,stored_gene_ids,new_genes,new_gene_ids,session,tag_map,db,cursor,MS_db)
            #print gene_mut_map_new,gene_mut_map_pseudo

            print "Before getSequences"
            #structure of gene_sequence_map: {Uniprot_Id:Sequence}
            gene_sequence_map,gene_mut_map_new,gene_mut_map_pseudo = getSequences(new_genes,stored_genes,fasta_map,species_map,u_ac_refseq_map,gene_mut_map_new,gene_mut_map_pseudo)
            #print gene_mut_map_new,gene_mut_map_pseudo

            print "Before autoTemplateSelection"
            #structure of pseudo_template_map: {Uniprot-Id:(template-list,{template-id:stored-template},oligo_map)}
            pseudo_template_map,gene_error_map,stored_annotations,background_process_AS = autoTemplateSelection(gene_sequence_map,stored_genes,new_genes,gene_mut_map_pseudo,AS_db)
            database.addErrorCodeToGene(gene_error_map,db,cursor)
            print len(stored_annotations)

            print "Before paraAlignment"
            #new structure of pseudo_template_map: {Uniprot-Id:({template-id:new template},{template-id:stored-template},oligo_map)}
            pseudo_template_map,gene_template_alignment_map = paraAlignment(gene_sequence_map,pseudo_template_map,stored_genes,new_genes,gene_mut_map_new,gene_mut_map_pseudo,stored_annotations)

            print len(stored_annotations)
            #print pseudo_template_map
            #print gene_template_alignment_map


            #structure of gene_template_alignment_map: {gene_id:{template_id:(aln_length,seq_id,sub_infos,alignment_pir)}}
            print "Before paraAnnotate"
            paraAnnotate(stored_genes,new_genes,gene_mut_map_new,gene_mut_map_pseudo,gene_template_alignment_map,pseudo_template_map,stored_annotations)
            
            t1 = time.time()
            #join the background inserts
            background_process_GS.join()
            background_process_MS.join()
            if background_process_AS != None:
                background_process_AS.join()

            t2 = time.time()
            print 'Resttime for background inserts: ',t2-t1

        #Error-Handling for a whole input line
        except:
            
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            #print "Pipeline Core Error: ",e,f,g
            errortext = '\n'.join([str(e),str(f),str(g)]) + '\n\n'
            f = open(errorlog,'a')
            f.write(errortext)
            f.close()
            No_Errors = False

        os.chdir(output_path)
        #"""
        try:
            shutil.rmtree(cwd)
        except:
            pass
        #"""

    if No_Errors:
        errortext  = "Finished without any error\n###############################################################################\n"
        f = open(errorlog,'a')
        f.write(errortext)
        f.close()
    else:
        errortext  = "###############################################################################\n"
        f = open(errorlog,'a')
        f.write(errortext)
        f.close()


    if newsession:
        endtime = SQLDateTime()
        database.updateSession(session,endtime,db,cursor)
    db.close()
    GS_db.close()
    MS_db.close()
    AS_db.close()

    tend = time.time()
    print(tend-t0)
    return session
