import MySQLdb
import os
import sys
import getopt
import time
import multiprocessing
from multiprocessing import Process, Queue, Manager, Value, Lock
import subprocess
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
neighborhood_calculation = False
error_annotations_into_db = True
anno_session_mapping = True
calculate_interaction_profiles=False

resources = 'manu'
num_of_cores = None
proc_n = 48
blast_processes = proc_n
alignment_processes = proc_n
annotation_processes = proc_n
number_of_processes = proc_n

dssp = True
pdb_input_asymetric_unit = False

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
dssp_path = ''
rin_db_path = ''
iupred_path = ''
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

    trunk = main_file_path.rsplit('/',1)[0]

    blast_db_path = '%s/lib/base/blast_db/pdbba' % trunk
    smiles_path = '%s/lib/base/ligand_bases/Components-smiles-stereo-oe.smi' % trunk
    inchi_path = '%s/lib/base/ligand_bases/inchi_base.tsv' % trunk
    human_id_mapping_path = '%s/lib/base/id_mapping' % trunk


def setPaths(output_path,config):
    global blast_path
    global pdb_path
    global annovar_path
    global dssp_path
    global rin_db_path
    global iupred_path
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
        elif opt == 'iupred_path':
            iupred_path = arg

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
    global neighborhood_calculation
    global calculate_interaction_profiles
    global pdb_input_asymetric_unit
    global resources
    global proc_n
    global blast_processes
    global alignment_processes
    global annotation_processes
    global number_of_processes

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
            elif opt == 'resources':
                resources = arg
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
            elif opt == 'neighborhood_calculation':
                if arg == 'True':
                    neighborhood_calculation = True
            elif opt == 'calculate_interaction_profiles':
                if arg == 'True':
                    calculate_interaction_profiles = True
                elif arg == 'False':
                    calculate_interaction_profiles = False
            elif opt == 'pdb_input_asymetric_unit':
                if arg == 'True':
                    pdb_input_asymetric_unit = True
                elif arg == 'False':
                    pdb_input_asymetric_unit = False

    if resources == 'auto' and num_of_cores == None:
        proc_n = multiprocessing.cpu_count() -1
        if proc_n <= 0:
            proc_n = 1

        blast_processes = proc_n
        alignment_processes = proc_n
        annotation_processes = proc_n
        number_of_processes = proc_n

    if num_of_cores != None:
        proc_n = num_of_cores
        
        blast_processes = proc_n
        alignment_processes = proc_n
        annotation_processes = proc_n
        number_of_processes = proc_n

    print 'Using %s core(s)' % str(proc_n)

def sequenceScan(genes):
    global pdb_path
    global db_adress
    global db_user_name
    global db_password

    sequenceScanGenes = set()
    sequenceScanPDB = set()
    for u_ac in genes:

        if u_ac.count(':') > 0 :
            if None in genes[u_ac][2]:
                sequenceScanPDB.add(u_ac)
                genes[u_ac][2] = set()
        else:
            if None in genes[u_ac][2]:
                sequenceScanGenes.add(u_ac)
                genes[u_ac][2] = set()

    try:
        db = MySQLdb.connect(db_adress,db_user_name,db_password,'struct_man_db_uniprot')
        cursor = db.cursor()
    except:
        db = None
        cursor = None

    if len(sequenceScanGenes) > 0:
        print "Amount of genes going into sequenceScan: ",len(sequenceScanGenes)

        gene_sequence_map = uniprot.getSequencesPlain(sequenceScanGenes,db,cursor)
        for u_ac in gene_sequence_map:
            if gene_sequence_map[u_ac] == 1 or gene_sequence_map[u_ac] == 0:
                print "Error in sequenceScan with gene: ",u_ac
                continue
            for (pos,aa) in enumerate(gene_sequence_map[u_ac]):
                #print pos,aa
                aac = '%s%s%s' % (aa,str(pos+1),aa)
                genes[u_ac][2].add(aac)

    if len(sequenceScanPDB) > 0:
        pdb_sequence_map = pdb.getSequencesPlain(sequenceScanPDB,pdb_path)
        for u_ac in pdb_sequence_map:
            for (pos,aa) in enumerate(pdb_sequence_map[u_ac]):
                #print pos,aa
                aac = '%s%s%s' % (aa,str(pos+1),aa)
                genes[u_ac][2].add(aac)

    if db != None:
        db.close()
    return genes,sequenceScanPDB

def buildQueue(filename,junksize,mrna_fasta=None):
    global auto_scale_template_percentage
    global species
    global human_id_mapping_path
    global db
    global cursor
    global db_adress
    global db_user_name
    global db_password

    t0 =time.time()

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
        if len(words) < 1:
            print "Skipped input line:\n%s\nToo few words.\n" % line
            continue
        sp_id = words[0]#.replace("'","\\'")
        try:
            if len(words) == 1 or words[1] == '':
                aachange = None
            else:
                aachange = words[1].replace("\n","")
                if not sp_id.count(':') == 1:

                    if ord(aachange[0]) > 47 and ord(aachange[0]) < 58: #if first char is a number
                        aa1 = 'X'
                        aachange = "X%s" % aachange
                    else:
                        aa1 = aachange[0]


                    if ord(aachange[-1]) > 47 and ord(aachange[-1]) < 58: #if last char is a number
                        aa2 = aachange[0]
                        pos = int(aachange[1:])
                        aachange = "%s%s" % (aachange,aa2)
                    else:
                        aa2 = aachange.split(',')[0][-1]
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
            #if len(sp_id) < 2:
            #    print line
            if sp_id[2] == "_":
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

            if not species in fasta_map:
                if species == 'multi':
                    fasta_map[species] = '%s/%s.fa' % (mrna_fasta,species)
                else:
                    fasta_map[species] = mrna_fasta

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

    #print global_species
    #print species
    #print id_map
    #print ac_map

    t1 = time.time()
    print "buildQueue Part 1: ",str(t1-t0)

    try:
        db_uniprot = MySQLdb.connect(db_adress,db_user_name,db_password,'struct_man_db_uniprot')
        cursor_uniprot = db_uniprot.cursor()
    except:
        db_uniprot = None
        cursor_uniprot = None
    genes,tag_map,species_map = uniprot.IdMapping(ac_map,id_map,np_map,db_uniprot,cursor_uniprot,tag_map,species_map)
    if db_uniprot != None:
        db_uniprot.close()

    t2 = time.time()
    print "buildQueue Part 2: ",str(t2-t1)

    if mrna_fasta == None:
        for pdb_chain_tuple in pdb_map:
            genes[pdb_chain_tuple] = pdb_map[pdb_chain_tuple]


    t3 = time.time()
    print "buildQueue Part 3: ",str(t3-t2)

    #detect human datasets and add some entries to the species map
    human_set = True
    for u_ac in genes:
        u_id = genes[u_ac][0]
        if u_id.count('HUMAN') == 0:
            human_set = False
        if not u_ac in species_map:
            if u_id.count('_') == 1:
                species_map[u_ac] = u_id.split('_')[1]

    if human_set:
        species = 'human'
    
    t4 = time.time()
    print "buildQueue Part 4: ",str(t4-t3)

    genes,corrected_input_pdbs = sequenceScan(genes)

    t5 = time.time()
    print "buildQueue Part 5: ",str(t5-t4)

    outlist = []
    #computation of auto_scale_template_percentage:
    s = len(genes)
    print "Total genes: ",s
    if s > 0:
        auto_scale_template_percentage = 1-0.5**(math.log10(s)-2)
        if auto_scale_template_percentage < 0.0:
            auto_scale_template_percentage = 0.0
        #print auto_scale_template_percentage
    
    #junksize is the maximal number of genes per cycle
    #return [genes],ac_id_map #if this not commented, all genes are taken in one cycle

    if s > junksize:

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
                (key,value) = genes.popitem()
                new_map[key] = value
            outlist.append(new_map)
        new_map = {}
        while len(genes) > 0:
            (key,value) = genes.popitem()
            new_map[key] = value
        if len(new_map) > 0:
            outlist.append(new_map)
    else:
        outlist.append(genes)

    t6 = time.time()
    print "buildQueue Part 6: ",str(t6-t5)

    return outlist,tag_map,species_map,fasta_map,corrected_input_pdbs

"""#seems unused, remove?
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
"""

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

def getSequences(new_genes,stored_genes,fasta_map,species_map,gene_mut_map_new,stored_gene_new_pos,IU_db,iupred_path,corrected_input_pdbs):
    global mrna_fasta
    global species
    global No_Errors
    global manager
    global lock
    global db
    global cursor
    global number_of_processes
    global human_id_mapping_path
    global pdb_path
    global pdb_input_asymetric_unit
    global blast_processes
    global cwd

    if mrna_fasta == None: #formerly: species == 'human', now using the full uniprot
        pdb_dict = {}
        for gene in new_genes:
            if len(gene) == 6 and gene.count(':') == 1:
                pdb_dict[gene] = new_genes[gene]
        for gene in stored_genes:
            if len(gene) == 6 and gene.count(':') == 1:
                pdb_dict[gene] = stored_genes[gene]

        t0 = time.time()
        try:
            db_uniprot = MySQLdb.connect(db_adress,db_user_name,db_password,'struct_man_db_uniprot')
            cursor_uniprot = db_uniprot.cursor()
        except:
            db_uniprot = None
            cursor_uniprot = None
        gene_sequence_map,gene_info_map = uniprot.getSequences(new_genes,'%s/human_info_map.tab' % human_id_mapping_path,db_uniprot,cursor_uniprot,pdb_dict)
        if db_uniprot != None:
            db_uniprot.close()
        t1 = time.time()
        print "Time for getSequences Part 1: %s" % str(t1-t0)

        pdb_sequence_map,pdb_pos_map = pdb.getSequences(pdb_dict,pdb_path,AU=pdb_input_asymetric_unit)
        t2 = time.time()
        print "Time for getSequences Part 2: %s" % str(t2-t1)
        gene_sequence_map.update(pdb_sequence_map)
        t3 = time.time()
        print "Time for getSequences Part 3: %s" % str(t3-t2)
        #print pdb_pos_map
        #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
        #structure of gene_mut_map_stored: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
        for u_id in gene_mut_map_new:
            if not u_id in pdb_pos_map:
                continue
            if u_id in corrected_input_pdbs: #Input pdbs for full sequence analysis have already corrected aaclists
                (gene_id,aac_map) = gene_mut_map_new[u_id]
                gene_info_map[gene_id] = ({},{},pdb_sequence_map[u_id])
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
            #structure of gene_info_map: {gene_id:({go_term_id:go_term_name},{reactome_id:pathway_name}),sequence}
            gene_info_map[gene_id] = ({},{},pdb_sequence_map[u_id])

        for u_id in stored_gene_new_pos:
            if not u_id in pdb_pos_map:
                continue
            if u_id in corrected_input_pdbs: #Input pdbs for full sequence analysis have already corrected aaclists
                continue
            (gene_id,aac_map) = stored_gene_new_pos[u_id]
            new_aac_map = {}
            for aac_base in aac_map:
                res_id = aac_base[1:]
                if not res_id in pdb_pos_map[u_id]:
                    print("Given res_id: %s not found in pdb_pos_map for: %s" % (res_id,u_id))
                    continue
                pos = pdb_pos_map[u_id][res_id]
                new_aac_base = '%s%s' % (aac_base[0],str(pos))
                new_aac_map[new_aac_base] = aac_map[aac_base]
            stored_gene_new_pos[u_id] = (gene_id,new_aac_map)

        t4 = time.time()
        print "Time for getSequences Part 4: %s" % str(t4-t3)
        database.addGeneInfos(gene_info_map,db,cursor)
        t5 = time.time()
        print "Time for getSequences Part 5: %s" % str(t5-t4)
        gene_sequence_map = database.getSequences(stored_genes,gene_sequence_map,db,cursor)
        t6 = time.time()
        print "Time for getSequences Part 6: %s" % str(t6-t5)
    else:
        gene_sequence_map = {}
        for species in fasta_map:
            mrna_fasta = fasta_map[species]
            f = open(mrna_fasta,'r')
            lines = f.readlines()
            f.close()

            #print mrna_fasta

            gene_seq_map = {}
            #g_gene_map = {}
            for line in lines:
                line = line.replace('\n','')
                #print line
                if line[0] == '>':
                    g_id = line.split()[0][1:]
                    #gene_name = line.split()[1]
                    #if species == "multi":
                    #    gene_name = g_id
                    gene_seq_map[g_id] = ''
                    #g_gene_map[g_id] = gene_name
                else:
                    gene_seq_map[g_id] += line

            gene_info_map = {}

            #print new_genes

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


    background_iu_process = None
    iupred_map = {}
    if iupred_path != '':
        """
        sys.path.append(iupred_path)
        import iupred2a
        for u_ac in gene_sequence_map:
            seq = gene_sequence_map[u_ac]
            iupred_score,glob_text = iupred2a.iupred(seq, 'long')
            print seq,iupred_score
        """
        
        t0 = time.time()

        in_queue = manager.Queue()
        out_queue = manager.Queue()

        for u_ac in gene_sequence_map:
            #iupred_map[u_ac] = ([],{}) #globular domains and residue->score dict
            seq = gene_sequence_map[u_ac]
            if seq == 0 or seq == 1 or seq == None:
                continue

            in_queue.put((u_ac,seq))
        processes = {}
        for i in range(1,blast_processes + 1):
            try:
                os.stat("%s/%d" %(cwd,i))
            except:
                os.mkdir("%s/%d" %(cwd,i))
            p = Process(target=paraIupred, args=(in_queue,out_queue,lock))
            processes[i] = p
            p.start()
        for i in processes:
            processes[i].join()

        out_queue.put(None)
        while True:
            out = out_queue.get()
            if out == None:
                break
            iupred_parts = out

            iupred_map[iupred_parts[0]] = (iupred_parts[1],iupred_parts[2])

        #print iupred_map

        t1 = time.time()
        print "Time for addIupred Part 1: %s" % str(t1-t0)
        try:
            background_iu_process = database.addIupred(iupred_map,gene_mut_map_new,stored_gene_new_pos,IU_db)
        except:
            background_iu_process = None
        t2 = time.time()
        print "Time for addIupred Part 2: %s" % str(t2-t1)

    return gene_sequence_map,gene_mut_map_new,stored_gene_new_pos,{},background_iu_process

def paraIupred(in_queue,out_queue,lock):
    with lock:
        in_queue.put(None)
    while True:
        with lock:
            inp = in_queue.get()
        if inp == None:
            return

        (u_ac,seq) = inp

        #print seq


        p = subprocess.Popen(['python3','%s/iupred2a.py' % iupred_path,seq,'glob'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err = p.communicate()
        #print out
        after_glob = False
        iupred_parts = [u_ac,[],{}]
        for line in out.split('\n'):
            if len(line) < 3:
                continue
            if line[0] == '#':
                continue
            if line[0] == '1':
                after_glob = True
            words = line.split()
            if not after_glob:
                if len(words) < 4:
                    continue
                if words[0] == 'globular' and words[1] == 'domain':
                    glob_ranges = words[3].split('-')
                    iupred_parts[1].append(glob_ranges)
            if after_glob:
                if len(words) < 3:
                    continue
                iupred_parts[2][words[0]] = (words[1],words[2])

        with lock:
            out_queue.put(iupred_parts)

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
            #print  gene,seq
            structures = blast.blast(seq,target_name,blast_path,blast_db_path,nr=option_number_of_templates,seq_thresh=option_seq_thresh,cov_thresh=option_ral_thresh,cwd = "%s/%d" %(cwd,i))
            #print gene,structures
            #print "Blast",len(structures)
            structures = templateSelection.selectTemplates(structures,pdb_path)
            #print "Select",len(structures)
            structures = templateFiltering.filterTemplates(structures,option_seq_thresh,option_res_thresh,option_ral_thresh)
            #print "Filter",len(structures)
            with lock:
                out_queue.put((gene,structures))
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            with lock:
                err_queue.put((e,f,g,gene))

def autoTemplateSelection(gene_sequence_map,new_genes):
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
    global blast_path
    global blast_db_path
    global No_Errors
    global anno_session_mapping

    t0 = time.time()

    process_list_db = set([])
    process_queue_blast = manager.Queue()
    blast_queues = [process_queue_blast]
    n = 0
    gene_error_map = {}

    #seq_map = {}

    #print gene_sequence_map,db,cursor,blast_path,blast_db_path

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
 
    t1 = time.time()

    #"""paralized blasting of single sequences  
    #print "before parablast"
    t12 = 0.0
    t3 = 0.0
    t2 = 0.0
    t4 = 0.0
    structure_map = {}
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

        template_map = {}
        out_queue.put(None)
        while True:
            out = out_queue.get()
            if out == None:
                break
            (gene,structures) = out
            structure_map[gene] = structures
        #"""
        del out_queue
        del process_queue_blast

        t3 += time.time()

    print "Template Selection Part 1: %s" % (str(t1-t0))
    print "Template Selection Part 2: %s" % (str(t2-t12))
    print "Template Selection Part 3: %s" % (str(t3-t2))

    #print structure_map

    #for gene in structure_map:
    #    print gene
    #    for tup in structure_map[gene]:
    #        print tup,structure_map[gene][tup]

    return structure_map,gene_error_map

def paraAlignment(gene_sequence_map,structure_map,stored_genes,new_genes,gene_mut_map_new,stored_gene_new_pos,pdb_pos_map):
    #structure of template_map: {Uniprot-Id:(template-list,{template-id:stored-template},oligo_map)}
    #structure of gene_sequence_map: {Uniprot_Id:Sequence}
    #structure of stored_genes: {Uniprot-Id:gene_id}
    #structure of new_genes: {Uniprot-Id:gene_id}
    #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #structure of gene_mut_map_stored: {Uniprot_Id:(gene_id,{AAC_base:mutation_id})}
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

    gene_id_gene_map = {}
    template_id_template_map = {}

    #new_structures = set()

    n = 0

    #temp_amount = 0

    for gene in structure_map:
        gene_id = new_genes[gene]
        gene_id_gene_map[gene_id] = gene
        structures = structure_map[gene]

        aaclist = gene_mut_map_new[gene][1]

        for (pdb_id,chain) in structures:
            input_queue.put((gene,gene_sequence_map[gene],pdb_id,chain,structures[(pdb_id,chain)],aaclist))
            #new_structures.add(pdb_id)
        if n == 1000:
            input_queue = manager.Queue()
            in_queues.append(input_queue)
            n = 0
        n += 1

    t1 = time.time()
    print "Alignment Part 1: %s" % (str(t1-t0))
    gene_structure_alignment_map,structure_id_map = database.getAlignments(stored_genes,db,cursor)
    t2 = time.time()
    print "Alignment Part 2: %s" % (str(t2-t1))

    input_queue = manager.Queue()
    out_queue = manager.Queue()

    #print stored_gene_new_pos

    for u_ac in stored_gene_new_pos:
        gene_id,aaclist = stored_gene_new_pos[u_ac]
        if not gene_id in gene_structure_alignment_map: #This happens for genes, for which we cannot find a structure
            continue
        #print u_ac
        if u_ac in pdb_pos_map:
            pos_res_map = dict((v,k) for k,v in pdb_pos_map[gene].iteritems())
        else:
            pos_res_map = {}
        #print u_ac
        for structure_id in gene_structure_alignment_map[gene_id]:
            structure_id_map[(pdb_id,chain)] = structure_id
            (target_seq,template_seq,coverage,seq_id) = gene_structure_alignment_map[gene_id][structure_id]
            #print structure_id,aaclist
            input_queue.put((aaclist,structure_id,pos_res_map,pdb_id,chain,target_seq,template_seq,gene_id))
    processes = {}
    for i in range(1,alignment_processes + 1):
        try:
            os.stat("%s/%d" %(cwd,i))
        except:
            os.mkdir("%s/%d" %(cwd,i))
        p = Process(target=paraMap, args=(input_queue,out_queue,lock))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()


    total_mappings = []
    mutation_updates = []

    out_queue.put(None)
    while True:
        out = out_queue.get()
        if out == None:
            break
        (mutation_updates_part,mappings) = out
        mutation_updates += mutation_updates_part
        total_mappings += mappings
        #print mappings

    t3 = time.time()
    print "Alignment Part 3: %s" % (str(t3-t2))

    #print "before alignment"

    template_sub_amount = 0
    t31 = 0.0
    t4 = 0.0
    t5 = 0.0
    t6 = 0.0
    t7 = 0.0
    t8 = 0.0
    t9 = 0.0
    gene_template_alignment_map = {}
    template_id_map = {}
    alignment_insertion_list = []
    structure_insertion_list = {}

    new_gene_stored_structure_mappings = []

    for input_queue in in_queues:
        t31 += time.time()
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
            errortext = "Alignment Error: %s, %s:%s\n%s\n%s\n%s\n\n" % (gene,pdb_id,chain,'\n'.join([str(e)]),str(f),str(g))
            f = open(errorlog,'a')
            f.write(errortext)
            f.close()
            No_Errors = False

        #print "after alignment"
        t4 += time.time()

        error_map = {}
        error_queue.put(None)
        while True:
            err = error_queue.get()
            if err == None:
                break
            (error,gene,pdb_id,chain) = err
            if not gene in error_map:
                error_map[gene] = {}
                #Only 1 error per gene
                error_map[gene][(pdb_id,chain)] = error

        for gene in error_map:
            for (pdb_id,chain) in error_map[gene]:
                errortext = "Parse error during alignment: %s - %s\n\n" % (gene,str(error_map[gene][(pdb_id,chain)]))
                f = open(errorlog,'a')
                f.write(errortext)
                f.close()

        out_queue.put(None)


        while True:
            out = out_queue.get()
            if out == None:
                break
            (gene,pdb_id,chain,structure,coverage,seq_id,sub_infos,alignment_pir,aaclist,update_map) = out
            structure_map[gene][(pdb_id,chain)] = structure
            #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
            #structure of gene_mut_map_stored: {Uniprot_Id:(gene_id,{AAC_base:mutation_id})}
            #structure of pdb_pos_map: pdb_pos_map[pdb_chain_tuple(=gene)][res_nr (in the pdb file)] = pos (in the sequence)

            if gene in pdb_pos_map:
                pos_res_map = dict((v,k) for k,v in pdb_pos_map[gene].iteritems())
            else:
                pos_res_map = {}
            #print pos_res_map
            if gene in gene_mut_map_new:
                gene_mut_map_new[gene] = (gene_mut_map_new[gene][0],aaclist)
            for old_aac_base in update_map:
                pos = old_aac_base[1:]
                database_pos = pos
                if int(pos) in pos_res_map:
                    database_pos = pos_res_map[int(pos)]
                new_aac_base = '%s%s' % (update_map[old_aac_base][1],pos)
                database_aac_base = '%s%s' % (update_map[old_aac_base][1],database_pos)
                mutation_id = aaclist[new_aac_base]
                mutation_updates.append((mutation_id,database_aac_base))

            if gene in stored_genes:
                gene_id = stored_genes[gene]
            else:
                gene_id = new_genes[gene]
            
            for aacbase in aaclist:
                if not aacbase in sub_infos: #this can happen for given positions > length of the sequence
                    continue
                m_id = aaclist[aacbase]
                res_id,t_aa = sub_infos[aacbase]
                new_gene_stored_structure_mappings.append((m_id,(pdb_id,chain),res_id,gene_id))

            if not gene in gene_template_alignment_map:
                gene_template_alignment_map[gene] = {}
            gene_template_alignment_map[gene][(pdb_id,chain)] = sub_infos

            #print pdb_id
            #print template_map[gene][2]

            #structure of template_map: {Uniprot-Id:(template-list,{template-id:stored-template},oligo_map)}


            if not (pdb_id,chain) in structure_insertion_list:
                structure_insertion_list[(pdb_id,chain)] = (structure['Resolution'],structure['Oligo'],structure['IAP'])
            alignment_insertion_list.append((gene_id,pdb_id,chain,structure['Seq_Id'],structure['Coverage'],alignment_pir))

        del input_queue
        del out_queue
        del error_queue

        t5 += time.time()

    database.updateMutations(mutation_updates,db,cursor)

    t6 += time.time()
    #structure of template_id_map[gene_id] = {pdb_id:template_id}
    #print structure_id_map
    structure_id_map,stored_structures = database.insertStructures(structure_insertion_list,db,cursor,smiles_path,inchi_path,pdb_path,structure_id_map)
    #print structure_id_map
    t7 += time.time()

    database.insertMapping(total_mappings,new_gene_stored_structure_mappings,stored_structures,gene_template_alignment_map,gene_mut_map_new,db,cursor)

    t8 += time.time()
    database.insertAlignments(alignment_insertion_list,structure_id_map,stored_structures,db,cursor)

    t9 += time.time()

    print "Amount of template-mutation pairs after alignment: ",template_sub_amount
        

    #after_amount = 0

    #print "Template amounts, before and after: ",temp_amount,after_amount
    print "Alignment Part 4: %s" % (str(t4-t31))
    print "Alignment Part 5: %s" % (str(t5-t4))
    print "Alignment Part 6: %s" % (str(t6-t5))
    print "Alignment Part 7: %s" % (str(t7-t6))
    print "Alignment Part 8: %s" % (str(t8-t7))
    print "Alignment Part 9: %s" % (str(t9-t8))

    return structure_map,gene_template_alignment_map,structure_id_map

def paraMap(input_queue,out_queue,lock):
    global pdb_path
    with lock:
        input_queue.put(None)
    while True:
        with lock:
            inp = input_queue.get()
        if inp == None:
            return

        (aaclist,structure_id,pos_res_map,pdb_id,chain,target_seq,template_seq,gene_id) = inp

        template_page = pdb.standardParsePDB(pdb_id,pdb_path)
        seq_res_map = globalAlignment.createTemplateFasta(template_page,pdb_id,chain,onlySeqResMap = True)
        sub_infos,errors,aaclist,update_map = globalAlignment.getSubPos(target_seq,template_seq,aaclist,seq_res_map)

        #structure of sub_infos: {aacbase:(res_id,template_aa)}
        #structure of aaclist: {aacbase:mutation_id}

        mutation_updates = []

        for old_aac_base in update_map:
            pos = old_aac_base[1:]
            database_pos = pos
            if int(pos) in pos_res_map:
                database_pos = pos_res_map[int(pos)]
            new_aac_base = '%s%s' % (update_map[old_aac_base][1],pos)
            database_aac_base = '%s%s' % (update_map[old_aac_base][1],database_pos)
            mutation_id = aaclist[new_aac_base]
            mutation_updates.append((mutation_id,database_aac_base))

        mappings = []
        for aacbase in aaclist:
            if not aacbase in sub_infos: #this can happen for given positions > length of the sequence
                continue
            m_id = aaclist[aacbase]
            res_id,t_aa = sub_infos[aacbase]
            mappings.append((m_id,structure_id,res_id,gene_id))

        with lock:
            out_queue.put((mutation_updates,mappings))

    return

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

        (gene,wildtype_sequence,pdb_id,chain,structure,aaclist) = inp

        try:
            t0 += time.time()

            (template_page,structure,error,part1_times) = pdb.getStandardizedPdbFile(pdb_id,chain,structure,pdb_path)
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
                #if pdb_id == '4JAN' and chain == 'I':
                #    print template_page
                (coverage,seq_id,sub_infos,alignment_pir,errors,times,aaclist,update_map) = globalAlignment.alignBioPython(gene,wildtype_sequence,pdb_id,template_page,chain,aaclist)
                structure['Seq_Id'] = seq_id
                structure['Coverage'] = coverage
                n = 0
                for ti in times:
                    aligntimes[n] += ti
                    n+=1

                for err in errors:
                    #[e,f,g] = sys.exc_info()
                    #g = traceback.format_exc(g)
                    with lock:
                        error_queue.put((err,gene,pdb_id,chain))
                #print "after alignment"

                with lock:
                    out_queue.put((gene,pdb_id,chain,structure,coverage,seq_id,sub_infos,alignment_pir,aaclist,update_map))
                #os.chdir(cwd)
            else:
                with lock:
                    error_queue.put((error,gene,pdb_id,chain))
            t2 += time.time()
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            with lock:
                err_queue.put((e,f,g,gene,pdb_id,chain))
    

def paraAnnotate(gene_mut_map_new,gene_template_alignment_map,structure_map,structure_id_map):
    global annotation_processes
    global manager
    global lock
    global db
    global cursor
    global option_lig_wf
    global option_chain_wf
    global cwd
    global session
    global No_Errors
    global anno_session_mapping
    global error_annotations_into_db

    global smiles_path
    global inchi_path
    global pdb_path

    #new structure of template_map: {Uniprot-Id:({template-id:new template},{template-id:stored-template},oligo_map)}
    #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #structure of gene_mut_map_stored: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #structure of gene_template_alignment_map: {Uniprot_Id:{PDB_Id:(coverage,seq_id,sub_infos)}}
    t0 = time.time()
    
    input_queue = manager.Queue()

    #new structure for paralell annotation: give each annotation process a pdb and dict of chain and structure information,
    #all chains without structure information shall go into the database as fully annotated structures, which are not aligned (yet)
    #these are unmapped chains, or protein interaction partners, their presence in the database is necesary for the computation of the structural neighborhood
    #a new structure-structure entry has to be created and has to be inserted into the database before the residue insertion
    #homooligomer information is not present for such cases, this information has to be updated, if the chain is mapped at a later iteration or run
    pdb_structure_map = {}
    size_sorted = {}
    max_size = 0

    for gene in structure_map:
        structures = structure_map[gene]
        for (pdb_id,chain) in structures:
            if not (pdb_id,chain) in structure_id_map:
                continue
            structure = structures[(pdb_id,chain)]
            s_id = structure_id_map[(pdb_id,chain)]
            if not pdb_id in pdb_structure_map:
                pdb_structure_map[pdb_id] = {}
                
                size = 0
                for iap in structure['IAP']:
                    if iap[0] != 'Ligand':
                        size += 1
                if size > max_size:
                    max_size = size
                if not size in size_sorted:
                    size_sorted[size] = []
                size_sorted[size].append(pdb_id)
            pdb_structure_map[pdb_id][chain] = (s_id,structure)

    #print size_sorted

    n_of_structs = 0
    while max_size > 0:
        if max_size in size_sorted:
            for pdb_id in size_sorted[max_size]:
                input_queue.put((pdb_id,pdb_structure_map[pdb_id]))
                n_of_structs += 1
        max_size -= 1

    t01 = time.time()
    print "Annotation Part 0.1: %s" % (str(t01-t0))
    

    total_annotations = {}
    interacting_structure_dict = {}

    out_queue = manager.Queue()
    error_queue = manager.Queue()
    err_queue = manager.Queue()
    processes = {}
    if annotation_processes > n_of_structs:
        annotation_processes = n_of_structs
    for i in range(1,annotation_processes + 1):
        try:
            os.stat("%s/%d" %(cwd,i))
        except:
            os.mkdir("%s/%d" %(cwd,i))
        p = Process(target=annotate, args=(input_queue,out_queue,error_queue,lock,cwd,i,err_queue,t01,i))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()

    t02 = time.time()
    print "Annotation Part 0.2: %s" % (str(t02-t01))

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

    t1 = time.time()
    print "Annotation Part 1: %s" % (str(t1-t02))
    
    with lock:
        out_queue.put(None)

    while True:
        out = out_queue.get()
        if out == None:
            break
        (annotation_chain_dict,pdb_id,chain_structure_map,interacting_chain_map,residue_residue_dict) = out
        #print len(annotations)

        for chain in interacting_chain_map:
            interacting_structure_dict[(pdb_id,chain)] = interacting_chain_map[chain]
        
        for chain in annotation_chain_dict:
            annotations = annotation_chain_dict[chain]

            total_annotations[(pdb_id,chain)] = annotations,residue_residue_dict[chain]
            

    t11 = time.time()
    print 'Annotation Part 1.1: %s' % (str(t11-t1))
    interacting_structure_ids = database.insertInteractingChains(interacting_structure_dict,db,cursor,smiles_path,inchi_path,pdb_path)


    t2  = time.time()
    print "Annotation Part 2: %s" % (str(t2-t11))
    structure_residue_map = database.insertResidues(total_annotations,db,cursor,structure_id_map,interacting_structure_ids)
    #print structure_residue_map

    t3  = time.time()

    print "Annotation Part 3: %s" % (str(t3-t2))
    #print gene_mut_map_new
    #print gene_template_alignment_map
    database.insertNewMappings(gene_mut_map_new,gene_template_alignment_map,structure_residue_map,structure_id_map,db,cursor)

    t4 = time.time()

    print "Annotation Part 4: %s" % (str(t4-t3))

def annotate(input_queue,out_queue,error_queue,lock,cwd,proc_number,err_queue,t01,i):
    global compute_surface
    global neighborhood_calculation
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
            #t02 = time.time()
            #print "Annotation Process %s: %s" % (str(i),str(t02-t01))
            return
        #print inp
        (pdb_id,chain_structure_map) = inp
        #print 'Proc %s - annotate: %s' % (str(i),pdb_id)
        try:
            annotation_chain_dict,interacting_chain_map,residue_residue_dict,errorlist = templateFiltering.structuralAnalysis(pdb_id,chain_structure_map,pdb_path,dssp_path,rin_db_path,neighborhood_calculation=neighborhood_calculation,dssp=dssp,calculate_interaction_profiles=calculate_interaction_profiles)
            with lock:
                out_queue.put((annotation_chain_dict,pdb_id,chain_structure_map,interacting_chain_map,residue_residue_dict))
                if len(errorlist) > 0:
                    for (error,e,f,g) in errorlist:
                        err_queue.put((error,f,g,pdb_id))
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            with lock:
                err_queue.put((e,f,g,pdb_id))

def main(filename,config,output_path,main_file_path,n_of_cores):
    global session
    global db
    global cursor
    global option_seq_thresh
    global option_res_thresh
    global option_ral_thresh
    global cwd
    global species
    global mrna_fasta
    global num_of_cores

    num_of_cores = n_of_cores

    t0 = time.time()

    setBasePaths(main_file_path)
    setPaths(output_path,config)

    parseConfig(config)



    if mrna_fasta != None and species == None:
        species = mrna_fasta.split('/')[-1].replace('.fa','').replace('.fasta','')
    if mrna_fasta != None:
        if not os.path.exists(mrna_fasta):
            raise NameError("mRNA path not found: %s" % mrna_fasta)

    #annovar-pipeline in case of vcf-file
    if filename.rsplit(".",1)[1] == "vcf":
        anno_db = "%s_annovar" % db_name.rsplit("_",1)[0]
        print 'Convert vcf file format using Annovar'
        nfname = annovar.annovar_pipeline(filename,tax_id,annovar_path,db_adress,db_user_name,db_password,anno_db,mrna_fasta,ref_id=ref_genome_id)
        #print nfname
    else:
        nfname = filename

    mrna_fasta_for_annovar = True #make this a option, for the case of mrna fasta files with non-standard protein-ids (may include further debugging)

    if mrna_fasta_for_annovar:
        mrna_fasta = None

    #print db_adress,db_user_name,db_password,db_name

    db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    MS_db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    #AS_db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    IU_db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    cursor = db.cursor()

    t01 = time.time()

    print "Time for preparation before buildQueue: %s" % (str(t01-t0))

    junksize = 500

    print "Call buildQueue with junksize: %s and file: %s" % (str(junksize),nfname)
    gene_aaclist_map_list,tag_map,species_map,fasta_map,corrected_input_pdbs = buildQueue(nfname,junksize,mrna_fasta=mrna_fasta)

    t02 = time.time()
    print "Time for buildQueue: %s" % (str(t02-t01))

    #sys.exit()

    #print gene_aaclist_map_list
    #print tag_map
    #print species_map
    #print fasta_map
    
    print "Number of junks: ",len(gene_aaclist_map_list)

    newsession = False
    if session == 0:     
        starttime = SQLDateTime()
        session = database.insertSession(starttime,nfname,db,cursor)
        newsession = True

    date = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
    errortext  = "###############################################################################\n%s:\n%s - Session: %s\n" % (date,nfname,str(session))
    No_Errors = True

    f = open(errorlog,'a')
    f.write(errortext)
    f.close()



    errors = []
    print errorlog
    junk_nr = 1 
    for gene_aaclist_map in gene_aaclist_map_list:
        if len(gene_aaclist_map) == 0:
            #print "empty gene_aaclist_map"
            continue

        print "Junk %s/%s" % (str(junk_nr),str(len(gene_aaclist_map_list)))
        junk_nr+=1
        try:
            os.stat("%s/tmp_structman_pipeline" %(output_path))
        except:
            os.mkdir("%s/tmp_structman_pipeline" %(output_path))  
        os.chdir("%s/tmp_structman_pipeline" %(output_path))
        cwd = "%s/tmp_structman_pipeline" %(output_path)

        try:
            t1 = time.time()
            print "Before geneCheck"
            #check for already stored genes and make all the necessary database interactions
            #structure of stored_genes: {Uniprot-Id:gene_id}
            #structure of new_genes: {Uniprot-Id:gene_id}
            stored_genes,new_genes,stored_gene_ids,new_gene_ids = database.geneCheck(gene_aaclist_map,species_map,session,db,cursor)
            #print stored_genes
            #print new_genes

            t2 = time.time()
            print "Time for geneCheck: %s" % (str(t2-t1))

            print "Before mutationCheck", len(stored_genes), len(new_genes)
            #check for already stored mutations or position twins and make all the necessary database interactions
            #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
            #structure of gene_mut_map_stored: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
            gene_mut_map_new,stored_gene_new_pos,background_process_MS = database.mutationCheck(gene_aaclist_map,stored_genes,stored_gene_ids,new_genes,new_gene_ids,session,tag_map,db,cursor,MS_db)

            #print gene_mut_map_new

            t3 = time.time()
            print "Time for mutationCheck: %s" % (str(t3-t2))

            print "Before getSequences"
            #structure of gene_sequence_map: {Uniprot_Id:Sequence}
            gene_sequence_map,gene_mut_map_new,stored_gene_new_pos,pdb_pos_map,background_iu_process = getSequences(new_genes,stored_genes,fasta_map,species_map,gene_mut_map_new,stored_gene_new_pos,IU_db,iupred_path,corrected_input_pdbs)

            #print gene_mut_map_new

            t4 = time.time()
            print "Time for getSequences: %s" % (str(t4-t3))

            print "Before autoTemplateSelection"
            #structure of template_map: {Uniprot-Id:(template-list,{template-id:stored-template},oligo_map)}
            structure_map,gene_error_map = autoTemplateSelection(gene_sequence_map,new_genes)
            database.addErrorCodeToGene(gene_error_map,db,cursor)
            #print structure_map

            t5 = time.time()
            print "Time for Template Selection: %s" % (str(t5-t4))

            print "Before paraAlignment"
            #new structure of template_map: {Uniprot-Id:({template-id:new template},{template-id:stored-template},oligo_map)}
            structure_map,gene_template_alignment_map,structure_id_map = paraAlignment(gene_sequence_map,structure_map,stored_genes,new_genes,gene_mut_map_new,stored_gene_new_pos,pdb_pos_map)

            t6 = time.time()
            print "Time for Alignment: %s" % (str(t6-t5))

            #structure of gene_template_alignment_map: {gene_id:{template_id:(coverage,seq_id,sub_infos,alignment_pir)}}
            print "Before paraAnnotate"
            paraAnnotate(gene_mut_map_new,gene_template_alignment_map,structure_map,structure_id_map)
            
            t7 = time.time()
            print "Time for Annotation: %s" % (str(t7-t6))

            t1 = time.time()
            #join the background inserts
            if background_process_MS != None:
                background_process_MS.join()
            #if background_process_AS != None:
            #    background_process_AS.join()

            if background_iu_process != None:
                background_iu_process.join()

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
    MS_db.close()
    #AS_db.close()
    IU_db.close()

    tend = time.time()
    print(tend-t0)
    return session
