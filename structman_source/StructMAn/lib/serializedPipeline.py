import pymysql as MySQLdb
import os
import sys
import getopt
import time
import multiprocessing
import subprocess
import shutil
import traceback
import math
import json

import pdbParser as pdb
import templateSelection
import templateFiltering
import globalAlignment
import database
import uniprot
import blast
import annovar
import MMseqs2


def SQLDateTime():
    return time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())

def sequenceScan(config,genes,tag_map):
    pdb_path = config.pdb_path
    db_adress = config.db_adress
    db_user_name = config.db_user_name
    db_password = config.db_password

    sequenceScanGenes = set()
    sequenceScanPDB = set()
    for u_ac in genes:

        if u_ac.count(':') > 0:
            if None in genes[u_ac][2]:
                genes[u_ac][2] = set()
            sequenceScanPDB.add(u_ac) #PDB inputs are always processed by the sequence scan
        else:
            if None in genes[u_ac][2]:
                genes[u_ac][2] = set()
                sequenceScanGenes.add(u_ac)

    try:
        db = MySQLdb.connect(db_adress,db_user_name,db_password,'struct_man_db_uniprot')
        cursor = db.cursor()
    except:
        db = None
        cursor = None

    if len(sequenceScanGenes) > 0:
        print("Amount of genes going into sequenceScan: ",len(sequenceScanGenes))

        gene_sequence_map = uniprot.getSequencesPlain(sequenceScanGenes,db,cursor,debug=config.verbose)
        for u_ac in gene_sequence_map:
            if gene_sequence_map[u_ac][0] == 1 or gene_sequence_map[u_ac][0] == 0:
                print("Error in sequenceScan with gene: ",u_ac)
                continue
            for (pos,aa) in enumerate(gene_sequence_map[u_ac][0]):
                aac = '%s%s%s' % (aa,str(pos+1),aa)
                genes[u_ac][2].add(aac)

    residue_id_backmap = {}
    if len(sequenceScanPDB) > 0:
        pdb_sequence_map,pdb_pos_map = pdb.getSequences(sequenceScanPDB,pdb_path)
        for u_ac in pdb_sequence_map:
            residue_id_backmap[u_ac] = {}
            res_pos_map = pdb_pos_map[u_ac]
            for res_id in res_pos_map:
                residue_id_backmap[u_ac][res_pos_map[res_id]] = res_id
            if len(genes[u_ac][2]) == 0: #The full set of positions is only taking, if the pdb-id is given in the input file without a specific position
                for (pos,aa) in enumerate(pdb_sequence_map[u_ac][0]):
                    aac = '%s%s%s' % (aa,str(pos+1),aa)
                    genes[u_ac][2].add(aac)
                
            else: #The given residue-ids are mapped onto their sequence position
                new_aacs = set()
                new_tags = {}
                for aac in genes[u_ac][2]:
                    res_id = aac[1:-1]
                    if not res_id in res_pos_map:
                        print('Skipped residue ',res_id,' from ',u_ac,' since illegal alternative location id')
                        continue
                    seq_pos = res_pos_map[res_id]
                    new_aac = '%s%s%s' % (aac[0],seq_pos,aac[-1])
                    new_aacs.add(new_aac)
                    new_tags[new_aac] = tag_map[u_ac][aac]
                genes[u_ac][2] = new_aacs
                tag_map[u_ac] = new_tags

    if db != None:
        db.close()
    return genes,sequenceScanPDB,residue_id_backmap,tag_map

def buildQueue(config,db,cursor,filename,junksize,mrna_fasta=None):
    auto_scale_template_percentage = config.auto_scale_template_percentage
    species = config.species

    db_adress = config.db_adress
    db_user_name = config.db_user_name
    db_password = config.db_password

    verbose = config.verbose

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
        if line == '':
            continue
        if len(line) < 3:
            print("Skipped input line:\n%s\nToo short.\n" % line)
            continue
        line = line.replace(' ','\t')
        words = line.split("\t")
        if len(words) < 1:
            print("Skipped input line:\n%s\nToo few words.\n" % line)
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

                else:
                    if aachange.count('_') == 1:
                        aachange = aachange.replace('_','')
                    else:
                        aachange = "%s%s" % (aachange,aachange[0])

        except:
            print("File Format Error: ",line)
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
                    pdb_chain_tuple = '%s:%s' % (sp_id[:4].upper(),sp_id[-1]) #enforce uppercase pdb-id 
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
                    print("Error: Identical gene ids for different species are not allowed: ",sp_id,species,species_map[sp_id])
            species_map[sp_id] = species

    t1 = time.time()
    if verbose:
        print("buildQueue Part 1: ",str(t1-t0))

    try:
        db_uniprot = MySQLdb.connect(db_adress,db_user_name,db_password,'struct_man_db_uniprot')
        cursor_uniprot = db_uniprot.cursor()
    except:
        db_uniprot = None
        cursor_uniprot = None
    genes,tag_map,species_map = uniprot.IdMapping(ac_map,id_map,np_map,db_uniprot,cursor_uniprot,tag_map,species_map,pdb_map,verbose=verbose)
    if db_uniprot != None:
        db_uniprot.close()

    t2 = time.time()
    if verbose:
        print("buildQueue Part 2: ",str(t2-t1))

    if mrna_fasta == None:
        for pdb_chain_tuple in pdb_map:
            genes[pdb_chain_tuple] = ['',set([]),pdb_map[pdb_chain_tuple]]


    t3 = time.time()
    if verbose:
        print("buildQueue Part 3: ",str(t3-t2))

    #detect human datasets and add some entries to the species map
    human_set = True
    for u_ac in genes:
        if u_ac in pdb_map:
            continue
        u_id = genes[u_ac][0]
        if u_id.count('HUMAN') == 0:
            human_set = False
        if not u_ac in species_map:
            if u_id.count('_') == 1:
                species_map[u_ac] = u_id.split('_')[1]

    if human_set:
        species = 'human'
    
    t4 = time.time()
    if verbose:
        print("buildQueue Part 4: ",str(t4-t3))

    genes,corrected_input_pdbs,residue_id_backmap,tag_map = sequenceScan(config,genes,tag_map)

    t5 = time.time()
    if verbose:
        print("buildQueue Part 5: ",str(t5-t4))

    outlist = []
    #computation of auto_scale_template_percentage:
    s = len(genes)
    print("Total proteins: ",s)
    if s > 0:
        auto_scale_template_percentage = 1-0.5**(math.log10(s)-2)
        if auto_scale_template_percentage < 0.0:
            auto_scale_template_percentage = 0.0
        #print auto_scale_template_percentage
    
    #junksize is the maximal number of genes per cycle
    #return [genes],ac_id_map #if this not commented, all genes are taken in one cycle

    if s > junksize:

        n_of_batches = s//junksize
        if s%junksize != 0:
            n_of_batches += 1
        batchsize = s//n_of_batches
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
    if verbose:
        print("buildQueue Part 6: ",str(t6-t5))

    return outlist,tag_map,species_map,fasta_map,corrected_input_pdbs,residue_id_backmap

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

def getSequences(config,new_genes,stored_genes,fasta_map,species_map,gene_mut_map_new,stored_gene_new_pos,corrected_input_pdbs,manager,lock,db,cursor):
    mrna_fasta = config.mrna_fasta
    species = config.species
    No_Errors = True

    number_of_processes = config.proc_n
    human_id_mapping_path = config.human_id_mapping_path
    pdb_path = config.pdb_path
    pdb_input_asymetric_unit = config.pdb_input_asymetric_unit
    blast_processes = config.blast_processes
    cwd = os.getcwd()
    verbose = config.verbose

    if mrna_fasta == None: #formerly: species == 'human', now using the full uniprot
        pdb_dict = {}
        for gene in new_genes:
            if len(gene) == 6 and gene.count(':') == 1:
                pdb_dict[gene] = new_genes[gene]
        for gene in stored_genes:
            if len(gene) == 6 and gene.count(':') == 1:
                pdb_dict[gene] = stored_genes[gene]


        if verbose:
            print('Get Sequences for ',len(new_genes)-len(pdb_dict),' unstored proteins')

        t0 = time.time()
        try:
            db_uniprot = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,'struct_man_db_uniprot')
            cursor_uniprot = db_uniprot.cursor()
        except:
            db_uniprot = None
            cursor_uniprot = None
        
        if db_uniprot != None:
            used_local_uniprot_db = True
            new_genes_c = new_genes.copy()
            new_genes_c.update(stored_genes)
            gene_sequence_map,gene_info_map = uniprot.getSequences(new_genes_c,'%s/human_info_map.tab' % human_id_mapping_path,db_uniprot,cursor_uniprot,pdb_dict)
            db_uniprot.close()
        else:
            used_local_uniprot_db = False
            gene_sequence_map,gene_info_map = uniprot.getSequences(new_genes,'%s/human_info_map.tab' % human_id_mapping_path,db_uniprot,cursor_uniprot,pdb_dict)
        t1 = time.time()
        if verbose:
            print("Time for getSequences Part 1: %s" % str(t1-t0))

        pdb_sequence_map,pdb_pos_map = pdb.getSequences(pdb_dict,pdb_path,AU=pdb_input_asymetric_unit)
        #print pdb_pos_map
        t2 = time.time()
        if verbose:
            print("Time for getSequences Part 2: %s" % str(t2-t1))
        gene_sequence_map.update(pdb_sequence_map)
        t3 = time.time()
        if verbose:
            print("Time for getSequences Part 3: %s" % str(t3-t2))
        #print pdb_pos_map
        #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
        #structure of gene_mut_map_stored: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
        for u_id in gene_mut_map_new:
            if not u_id in pdb_pos_map:
                continue
            if u_id in corrected_input_pdbs: #Input pdbs for full sequence analysis have already corrected aaclists
                (gene_id,aac_map) = gene_mut_map_new[u_id]
                gene_info_map[gene_id] = ({},{},pdb_sequence_map[u_id][0])
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
            gene_info_map[gene_id] = ({},{},pdb_sequence_map[u_id][0])

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
                    print(("Given res_id: %s not found in pdb_pos_map for: %s" % (res_id,u_id)))
                    continue
                pos = pdb_pos_map[u_id][res_id]
                new_aac_base = '%s%s' % (aac_base[0],str(pos))
                new_aac_map[new_aac_base] = aac_map[aac_base]
            stored_gene_new_pos[u_id] = (gene_id,new_aac_map)

        t4 = time.time()
        if verbose:
            print("Time for getSequences Part 4: %s" % str(t4-t3))
        database.addGeneInfos(gene_info_map,db,cursor)
        t5 = time.time()
        if verbose:
            print("Time for getSequences Part 5: %s" % str(t5-t4))
        if not used_local_uniprot_db:
            gene_sequence_map = database.getSequences(stored_genes,gene_sequence_map,db,cursor)
        t6 = time.time()
        if verbose:
            print("Time for getSequences Part 6: %s" % str(t6-t5))
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
                    gene_sequence_map[gene] = sequence,None,None
                    gene_info_map[new_genes[gene]] = ([],{},{},sequence)
            database.addGeneInfos(gene_info_map,db,cursor)
            gene_sequence_map.update(database.getSequences(stored_genes,gene_sequence_map,db,cursor))


    background_iu_process = None
    iupred_map = {}
    if config.iupred_path != '':
        iupred_path = config.iupred_path
        if iupred_path.count('mobidb-lite') > 0:
            mobi_lite = True
        else:
            mobi_lite = False

        t0 = time.time()

        if not mobi_lite:
            in_queue = manager.Queue()
            out_queue = manager.Queue()
        else:
            mobi_list = []

        stored_disorder_ids = []
        for u_ac in gene_sequence_map:
            #iupred_map[u_ac] = ([],{}) #globular domains and residue->score dict
            seq,disorder_scores,disorder_regions = gene_sequence_map[u_ac]
            if seq == 0 or seq == 1 or seq == None:
                continue
            if disorder_scores == None:
                if not mobi_lite:
                    in_queue.put((u_ac,seq))
                else:
                    mobi_list.append('>%s\n%s\n' % (u_ac,seq))
            elif disorder_scores != 'Stored':
                disorder_scores_datastruct = {}
                for pos,score in enumerate(disorder_scores):
                    seq_pos = pos + 1
                    if pos >= len(seq):
                        if verbose:
                            print('Warning: illegal position ',pos,' for sequence of ',u_ac)
                        continue
                    disorder_scores_datastruct[seq_pos] = (seq[pos],score)
                iupred_map[u_ac] = (disorder_regions,disorder_scores_datastruct,'MobiDB3.0')

        if not mobi_lite:
            processes = {}
            for i in range(1,blast_processes + 1):
                try:
                    os.stat("%s/%d" %(cwd,i))
                except:
                    os.mkdir("%s/%d" %(cwd,i))
                p = multiprocessing.Process(target=paraIupred, args=(config,in_queue,out_queue,lock))
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

                iupred_map[iupred_parts[0]] = (iupred_parts[1],iupred_parts[2],'IUpred')

        else:
            mobi_tmp_file = 'mobi_tmp_file.fasta'
            f = open(mobi_tmp_file,'w')
            f.write(''.join(mobi_list))
            f.close()
            mobi_bin_path = '%s/binx/' % iupred_path.rsplit('/',1)[0]
            mobi_threads = min([7,config.proc_n])
            p = subprocess.Popen([iupred_path,mobi_tmp_file,'-t',str(mobi_threads),'-bin',mobi_bin_path,'-l'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
            out,err = p.communicate()
            os.remove(mobi_tmp_file)
            if err != '':
                print('Warning: mobidb-lite threw an error: ',err)
            else:
                entries = []
                for line in out.split('}'):
                    if line.count('{') == 0:
                        continue
                    line = line + '}'
                    #print(line)
                    mobi_results = json.loads(line)
                    u_ac = mobi_results['acc']
                    raw_scores = mobi_results['p']
                    raw_regions = mobi_results['regions']
                    regions = []
                    for a,b in raw_regions:
                        regions.append([a,b,'disorder'])
                    scores = {}
                    seq,disorder_scores,disorder_regions = gene_sequence_map[u_ac]
                    for pos,score in enumerate(raw_scores):
                        scores[pos+1] = (seq[pos],score)
                    iupred_map[u_ac] = (regions,scores,'mobidb-lite')

        t1 = time.time()
        if verbose:
            print("Time for addIupred Part 1: %s" % str(t1-t0))

        background_iu_process = database.addIupred(iupred_map,gene_mut_map_new,stored_gene_new_pos,config)

        t2 = time.time()
        if verbose:
            print("Time for addIupred Part 2: %s" % str(t2-t1))

    return gene_sequence_map,gene_mut_map_new,stored_gene_new_pos,{},background_iu_process,No_Errors

def paraIupred(config,in_queue,out_queue,lock):
    iupred_path = config.iupred_path
    with lock:
        in_queue.put(None)
    while True:
        with lock:
            inp = in_queue.get()
        if inp == None:
            return

        (u_ac,seq) = inp

        #print seq


        p = subprocess.Popen(['python3','%s/iupred2a.py' % iupred_path,seq,'glob'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
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
                    lower_bound,upper_bound = words[3].split('-')
                    iupred_parts[1].append((int(lower_bound),int(upper_bound),'globular'))
            if after_glob:
                if len(words) < 3:
                    continue
                iupred_parts[2][int(words[0])] = (words[1],words[2])

        with lock:
            out_queue.put(iupred_parts)

def paraBlast(config,process_queue_blast,out_queue,lock,i,err_queue):
    blast_path = config.blast_path
    blast_db_path = config.blast_db_path
    option_number_of_templates = config.option_number_of_templates
    option_seq_thresh = config.option_seq_thresh
    option_ral_thresh = config.option_ral_thresh
    option_res_thresh = config.option_res_thresh
    pdb_path = config.pdb_path
    cwd = os.getcwd()
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
            g = traceback.format_exc()
            with lock:
                err_queue.put((e,f,g,gene))

def autoTemplateSelection(config,db,cursor,manager,lock,No_Errors,gene_sequence_map,new_genes,search_tool='MMseqs2'):
    #structure of stored_genes: {Uniprot-Id:(gene_id,more_restrictive)}
    #structure of new_genes: {Uniprot-Id:gene_id}

    blast_processes = config.blast_processes
    option_seq_thresh = config.option_seq_thresh
    option_res_thresh = config.option_res_thresh
    cwd = os.getcwd()
    pdb_path = config.pdb_path
    errorlog = config.errorlog_path

    mmseqs2_path = config.mmseqs2_path
    mmseqs2_db_path = config.mmseqs2_db_path

    verbose = config.verbose

    if verbose:
        print('Sequence search with: ',search_tool)

    if search_tool == 'Blast':

        t0 = time.time()

        process_list_db = set([])
        process_queue_blast = manager.Queue()
        blast_queues = [process_queue_blast]
        n = 0
        gene_error_map = {}

        #seq_map = {}

        #print gene_sequence_map,db,cursor,blast_path,blast_db_path

        for gene in new_genes:
            seq = gene_sequence_map[gene][0]
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
                p = multiprocessing.Process(target=paraBlast, args=(config,process_queue_blast,out_queue,lock,i,err_queue))
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

        if verbose:
            print("Template Selection Part 1: %s" % (str(t1-t0)))
            print("Template Selection Part 2: %s" % (str(t2-t12)))
            print("Template Selection Part 3: %s" % (str(t3-t2)))

        #print structure_map

        #for gene in structure_map:
        #    print gene
        #    for tup in structure_map[gene]:
        #        print tup,structure_map[gene][tup]
    elif search_tool == 'MMseqs2':
        t0 = time.time()
        gene_error_map = {}
        filtered_gene_seq_map = {}
        for gene in new_genes:
            filtered_gene_seq_map[gene] = gene_sequence_map[gene]
        if len(filtered_gene_seq_map) > 0:
            mmseqs_tmp_folder = config.mmseqs_tmp_folder
            raw_structure_map,pdb_ids = MMseqs2.search(filtered_gene_seq_map,mmseqs2_db_path,mmseqs2_path,mmseqs_tmp_folder,option_seq_thresh,verbose=verbose)

            t1 = time.time()

            structure_map = templateSelection.filterRawStructureMap(raw_structure_map,pdb_ids,pdb_path,option_res_thresh,blast_processes)

            t2 = time.time()
            if verbose:
                print("Template Selection Part 1: %s" % (str(t1-t0)))
                print("Template Selection Part 2: %s" % (str(t2-t1)))
        else:
            structure_map = {}

    return structure_map,gene_error_map,No_Errors

def paraAlignment(config,db,cursor,manager,lock,No_Errors,gene_sequence_map,structure_map,stored_genes,stored_gene_ids,new_genes,gene_mut_map_new,stored_gene_new_pos,pdb_pos_map):
    #structure of template_map: {Uniprot-Id:(template-list,{template-id:stored-template},oligo_map)}
    #structure of gene_sequence_map: {Uniprot_Id:Sequence}
    #structure of stored_genes: {Uniprot-Id:gene_id}
    #structure of new_genes: {Uniprot-Id:gene_id}
    #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #structure of gene_mut_map_stored: {Uniprot_Id:(gene_id,{AAC_base:mutation_id})}
    
    alignment_processes = config.alignment_processes
    pdb_path = config.pdb_path
    smiles_path = config.smiles_path
    inchi_path = config.inchi_path
    verbose = config.verbose
    errorlog = config.errorlog_path
    t0 = time.time()
    input_queue = manager.Queue()
    err_queue = manager.Queue()

    in_queues = [input_queue] #partioning the parallized alignment part shall reduce the memory usage, by safing the alignments to the database every 1000 genes.

    gene_id_gene_map = {}
    template_id_template_map = {}

    #new_structures = set()

    n = 0

    for gene in structure_map:
        gene_id = new_genes[gene]
        gene_id_gene_map[gene_id] = gene
        structures = structure_map[gene]

        aaclist = gene_mut_map_new[gene][1]

        for (pdb_id,chain) in structures:
            input_queue.put((gene,gene_sequence_map[gene][0],pdb_id,chain,structures[(pdb_id,chain)],aaclist))
            #new_structures.add(pdb_id)
        if n == 1000:
            input_queue = manager.Queue()
            in_queues.append(input_queue)
            n = 0
        n += 1

    t1 = time.time()
    if verbose:
        print("Alignment Part 1: %s" % (str(t1-t0)))
    gene_structure_alignment_map,structure_id_map,id_structure_map = database.getAlignments(stored_gene_ids,db,cursor,verbose=verbose)
    t2 = time.time()
    if verbose:
        print("Alignment Part 2: %s" % (str(t2-t1)))

    input_queue = manager.Queue()
    out_queue = manager.Queue()

    for u_ac in stored_gene_new_pos:
        gene_id,aaclist = stored_gene_new_pos[u_ac]
        if not gene_id in gene_structure_alignment_map: #This happens for genes, for which we cannot find a structure
            continue
        #print u_ac
        if u_ac in pdb_pos_map:
            pos_res_map = dict((v,k) for k,v in pdb_pos_map[gene].items())
        else:
            pos_res_map = {}
        #print u_ac
        if not u_ac in structure_map:
            structure_map[u_ac] = {}
        for structure_id in gene_structure_alignment_map[gene_id]:
            pdb_id,chain = id_structure_map[structure_id]
            (target_seq,template_seq,coverage,seq_id) = gene_structure_alignment_map[gene_id][structure_id]
            structure_map[u_ac][(pdb_id,chain)] = {'Coverage':coverage,'Seq_Id':seq_id}
            input_queue.put((aaclist,structure_id,pos_res_map,pdb_id,chain,target_seq,template_seq,gene_id))
    processes = {}
    for i in range(1,alignment_processes + 1):
        p = multiprocessing.Process(target=paraMap, args=(config,input_queue,out_queue,lock))
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

    t3 = time.time()
    if verbose:
        print("Alignment Part 3: %s" % (str(t3-t2)))

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
            p = multiprocessing.Process(target=align, args=(config,input_queue,out_queue,error_queue,lock,err_queue))
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
                pos_res_map = dict((v,k) for k,v in pdb_pos_map[gene].items())
            else:
                pos_res_map = {}

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
    structure_id_map,stored_structures,structure_dict = database.insertStructures(structure_insertion_list,db,cursor,smiles_path,inchi_path,pdb_path,structure_id_map)

    t7 += time.time()

    if verbose:
        print("Amount of mappings based on stored structures: ",len(total_mappings))

    m_r_map,residue_dict = database.getStoredResidues(total_mappings,new_gene_stored_structure_mappings,stored_structures,gene_template_alignment_map,gene_mut_map_new,db,cursor,verbose=verbose)

    t8 += time.time()
    database.insertAlignments(alignment_insertion_list,structure_id_map,stored_structures,db,cursor)

    t9 += time.time()

    #after_amount = 0

    #print "Template amounts, before and after: ",temp_amount,after_amount
    if verbose:
        print("Alignment Part 4: %s" % (str(t4-t31)))
        print("Alignment Part 5: %s" % (str(t5-t4)))
        print("Alignment Part 6: %s" % (str(t6-t5)))
        print("Alignment Part 7: %s" % (str(t7-t6)))
        print("Alignment Part 8: %s" % (str(t8-t7)))
        print("Alignment Part 9: %s" % (str(t9-t8)))

    return structure_map,gene_template_alignment_map,structure_id_map,No_Errors,structure_dict,m_r_map,residue_dict,stored_structures

def paraMap(config,input_queue,out_queue,lock):
    pdb_path = config.pdb_path
    with lock:
        input_queue.put(None)
    while True:
        with lock:
            inp = input_queue.get()
        if inp == None:
            return

        (aaclist,structure_id,pos_res_map,pdb_id,chain,target_seq,template_seq,gene_id) = inp

        template_page = pdb.standardParsePDB(pdb_id,pdb_path,obsolete_check=True)
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
            if res_id == None: #This happens, when the position is mapped to a gap in the alignment
                continue
            mappings.append((m_id,structure_id,res_id,gene_id))

        with lock:
            out_queue.put((mutation_updates,mappings))

    return

def align(config,input_queue,out_queue,error_queue,lock,err_queue):
    pdb_path = config.pdb_path
    option_seq_thresh = config.option_seq_thresh
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

        #print pdb_id,chain

        try:
            t0 += time.time()

            (template_page,structure,error,part1_times) = pdb.getStandardizedPdbFile(pdb_id,chain,structure,pdb_path)
            n = 0
            for ti in part1_times:
                parsetimes[n] += ti
                n+=1

            t1 += time.time()
            if error == None:

                (coverage,seq_id,sub_infos,alignment_pir,errors,times,aaclist,update_map) = globalAlignment.alignBioPython(gene,wildtype_sequence,pdb_id,template_page,chain,aaclist)

                structure['Seq_Id'] = seq_id
                structure['Coverage'] = coverage
                n = 0
                for ti in times:
                    aligntimes[n] += ti
                    n+=1

                for err in errors:
                    #[e,f,g] = sys.exc_info()
                    #g = traceback.format_exc()
                    with lock:
                        error_queue.put((err,gene,pdb_id,chain))
                #print "after alignment"

                if 100.0*seq_id >= option_seq_thresh:
                    with lock:
                        out_queue.put((gene,pdb_id,chain,structure,coverage,seq_id,sub_infos,alignment_pir,aaclist,update_map))
            else:
                with lock:
                    error_queue.put((error,gene,pdb_id,chain))
            t2 += time.time()
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc()
            with lock:
                err_queue.put((e,f,g,gene,pdb_id,chain))
    

def paraAnnotate(config,db,cursor,manager,lock,No_Errors,gene_mut_map_new,gene_mut_map_stored,gene_template_alignment_map,structure_map,structure_id_map,structure_dict,m_r_map,stored_residue_dict,stored_structures):
    annotation_processes = config.annotation_processes
    anno_session_mapping = config.anno_session_mapping
    error_annotations_into_db = config.error_annotations_into_db

    smiles_path = config.smiles_path
    inchi_path = config.inchi_path
    pdb_path = config.pdb_path
    verbose = config.verbose
    errorlog = config.errorlog_path

    #new structure of template_map: {Uniprot-Id:({template-id:new template},{template-id:stored-template},oligo_map)}
    #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #structure of gene_mut_map_stored: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #structure of gene_template_alignment_map: {Uniprot_Id:{(PDB_Id,Chain):(csub_infos)}}
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
    gene_structure_map = {}
    pdb_ids = set()
    for gene in structure_map:
        structures = structure_map[gene]
        if gene in gene_mut_map_new:
            g_id = gene_mut_map_new[gene][0]
        else:
            g_id = gene_mut_map_stored[gene][0]
        gene_structure_map[g_id] = {}
        deletion_list = []
        for (pdb_id,chain) in structures:
            if not (pdb_id,chain) in structure_id_map: #Structures, that where filtered after the alignment
                deletion_list.append((pdb_id,chain))
                continue
            pdb_ids.add(pdb_id)
            s_id = structure_id_map[(pdb_id,chain)]
            structure = structures[(pdb_id,chain)]
            cov = structure['Coverage']
            seq_id = structure['Seq_Id']
            gene_structure_map[g_id][s_id] = (seq_id,cov)

            if 'IAP' not in structure:
                #del structure_map[gene][(pdb_id,chain)]
                continue
            if (pdb_id,chain) in stored_structures:
                continue

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
        for (pdb_id,chain) in deletion_list:
            del structure_map[gene][(pdb_id,chain)]

    '''
    sys.path.append("/TL/sin/work/agress/RIN")
    import createRINdb
    createRINdb.calculateRINsFromPdbList(pdb_structure_map.keys(),fromScratch=True,forceCentrality=True)
    '''

    n_of_structs = 0
    while max_size > 0:
        if max_size in size_sorted:
            for pdb_id in size_sorted[max_size]:
                input_queue.put((pdb_id,pdb_structure_map[pdb_id]))
                n_of_structs += 1
        max_size -= 1

    if verbose:
        t01 = time.time()
        print("Annotation Part 0.1: %s" % (str(t01-t0)))
    

    total_annotations = {}
    interacting_structure_dict = {}
    complex_profiles = {}

    out_queue = manager.Queue()
    error_queue = manager.Queue()
    err_queue = manager.Queue()
    processes = {}
    if annotation_processes > n_of_structs:
        annotation_processes = n_of_structs
    if verbose:
        print("Going into Annotation with ",annotation_processes,' processes')
    for i in range(1,annotation_processes + 1):
        p = multiprocessing.Process(target=annotate, args=(config,input_queue,out_queue,error_queue,lock,err_queue,t01))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()

    if verbose:
        t02 = time.time()
        print("Annotation Part 0.2: %s" % (str(t02-t01)))

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

    if verbose:
        t1 = time.time()
        print("Annotation Part 1: %s" % (str(t1-t02)))
    
    with lock:
        out_queue.put(None)

    while True:
        out = out_queue.get()
        if out == None:
            break
        (annotation_chain_dict,pdb_id,chain_structure_map,interacting_chain_map,residue_residue_dict,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles) = out
        for chain in interacting_chain_map:
            interacting_structure_dict[(pdb_id,chain)] = interacting_chain_map[chain]
        
        for chain in annotation_chain_dict:
            annotations = annotation_chain_dict[chain]

            total_annotations[(pdb_id,chain)] = annotations,residue_residue_dict[chain]

        complex_profiles[pdb_id] = ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles
            

    if verbose:
        t11 = time.time()
        print('Annotation Part 1.1: %s' % (str(t11-t1)))
    interacting_structure_ids = database.insertInteractingChains(interacting_structure_dict,db,cursor,smiles_path,inchi_path,pdb_path)

    if verbose:
        t12 = time.time()
        print('Annotation Part 1.2: %s' % (str(t12-t11)))
    complex_profile = database.insertComplexes(complex_profiles,pdb_ids,db,cursor)

    if verbose:
        t2  = time.time()
        print("Annotation Part 2: %s" % (str(t2-t12)))

    structure_residue_map,residue_dict = database.insertResidues(total_annotations,db,cursor,structure_id_map,interacting_structure_ids,stored_structures)
    residue_dict.update(stored_residue_dict)

    if verbose:
        t21  = time.time()
        print("Annotation Part 2.1: %s" % (str(t21-t2)))

    if verbose:
        t3  = time.time()
        print("Annotation Part 3: %s" % (str(t3-t21)))

        print('Res_dict size before updateResidues: ',len(residue_dict))

    m_r_map,residue_dict = updateResidues(gene_mut_map_new,gene_template_alignment_map,structure_residue_map,structure_id_map,db,cursor,m_r_map,residue_dict)

    residue_dict.update(stored_residue_dict)

    if verbose:
        print('Res_dict size after: ',len(residue_dict))

    if verbose:
        t4 = time.time()
        print("Annotation Part 4: %s" % (str(t4-t3)))

    database.insertClassifications(db,cursor,residue_dict,structure_dict,gene_structure_map,m_r_map,config,complex_profiles)

    if verbose:
        t5 = time.time()
        print("Annotation Part 5: %s" % (str(t5-t4)))

    return No_Errors

def updateResidues(gene_mut_map_new,gene_template_alignment_map,structure_residue_map,structure_id_map,db,cursor,m_r_map,residue_dict):
    #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #structure of gene_template_alignment_map: {Uniprot_Id:{(PDB_Id,chain):sub_infos}}
    #structute of structure_residue_map: {s_id:{res_nr:r_id}}

    value_strs = []
    condensed_residue_dict = {}

    for u_ac in gene_template_alignment_map:
        for (pdb_id,chain) in gene_template_alignment_map[u_ac]:
            if not (pdb_id,chain) in structure_id_map:
                continue
            s_id = structure_id_map[(pdb_id,chain)]
            g_id = gene_mut_map_new[u_ac][0]

            sub_infos = gene_template_alignment_map[u_ac][(pdb_id,chain)]
            for aacbase in sub_infos:

                res_nr,t_aa = sub_infos[aacbase]
                if res_nr == None:
                    continue
                m_id = gene_mut_map_new[u_ac][1][aacbase]
                if not s_id in structure_residue_map:
                    print('Warning: Structure ',s_id,'not in structure_residue_map, u_ac: ',u_ac,',PDB,Chain: ',pdb_id,chain)
                    continue
                if not res_nr in structure_residue_map[s_id]:
                    print('Warning: res_nr ',res_nr,'not in structure_residue_map, u_ac,s_id,aacb,m_id:',u_ac,s_id,aacbase,m_id)
                    continue
                r_id = structure_residue_map[s_id][res_nr]
                if not m_id in m_r_map:
                    m_r_map[m_id] = set()
                m_r_map[m_id].add(r_id)
                condensed_residue_dict[r_id] = residue_dict[r_id]
    return m_r_map,condensed_residue_dict

def annotate(config,input_queue,out_queue,error_queue,lock,err_queue,t01):
    neighborhood_calculation = config.neighborhood_calculation
    calculate_interaction_profiles = config.calculate_interaction_profiles
    dssp = config.dssp
    dssp_path = config.dssp_path
    pdb_path = config.pdb_path
    rin_db_path = config.rin_db_path

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
            annotation_chain_dict,interacting_chain_map,residue_residue_dict,errorlist,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles = templateFiltering.structuralAnalysis(pdb_id,chain_structure_map,pdb_path,dssp_path,rin_db_path,neighborhood_calculation=neighborhood_calculation,dssp=dssp,calculate_interaction_profiles=calculate_interaction_profiles)
            with lock:
                out_queue.put((annotation_chain_dict,pdb_id,chain_structure_map,interacting_chain_map,residue_residue_dict,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles))
                if len(errorlist) > 0:
                    for (error,e,f,g) in errorlist:
                        err_queue.put((error,f,g,pdb_id))
            #print 'finished Proc %s - annotate: %s' % (str(i),pdb_id)
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc()
            with lock:
                err_queue.put((e,f,g,pdb_id))

def main(filename,config,output_path,main_file_path):
    n_of_cores = config.proc_n
    species = config.species
    mrna_fasta = config.mrna_fasta
    num_of_cores = config.proc_n
    verbose = config.verbose
    search_tool = config.search_tool
    errorlog = config.errorlog_path
    session = 0 #This can later be used, structman.py could give a specific session id and the pipeline can then expand that session

    background_process_MS = None

    manager = multiprocessing.Manager()
    lock = manager.Lock()

    t0 = time.time()

    if mrna_fasta != None and species == None:
        species = mrna_fasta.split('/')[-1].replace('.fa','').replace('.fasta','')
    if mrna_fasta != None:
        if not os.path.exists(mrna_fasta):
            raise NameError("mRNA path not found: %s" % mrna_fasta)

    #annovar-pipeline in case of vcf-file
    if filename.rsplit(".",1)[1] == "vcf":
        anno_db = "%s_annovar" % db_name.rsplit("_",1)[0]
        print('Convert vcf file format using Annovar')
        if mrna_fasta != None:
            '... and using mrna file: ',mrna_fasta
        nfname = annovar.annovar_pipeline(filename,config.tax_id,config.annovar_path,config.db_adress,config.db_user_name,config.db_password,anno_db,mrna_fasta,ref_id=config.ref_genome_id)
    else:
        nfname = filename

    mrna_fasta_for_annovar = True #make this a option, for the case of mrna fasta files with non-standard protein-ids (may include further debugging)

    if mrna_fasta_for_annovar:
        mrna_fasta = None

    db = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,config.db_name)
    cursor = db.cursor()

    t01 = time.time()

    if verbose:
        print("Time for preparation before buildQueue: %s" % (str(t01-t0)))

    junksize = 500

    if verbose:
        print("Call buildQueue with chunksize: %s and file: %s" % (str(junksize),nfname))
    gene_aaclist_map_list,tag_map,species_map,fasta_map,corrected_input_pdbs,residue_id_backmap = buildQueue(config,db,cursor,nfname,junksize,mrna_fasta=mrna_fasta)

    t02 = time.time()
    if verbose:
        print("Time for buildQueue: %s" % (str(t02-t01)))
    
    print("Number of chunks: ",len(gene_aaclist_map_list))

    newsession = False
    if session == 0:     
        starttime = SQLDateTime()
        session = database.insertSession(starttime,nfname,db,cursor)
        newsession = True

    date = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
    errortext  = "###############################################################################\n%s:\n%s - Session: %s\n" % (date,nfname,str(session))

    f = open(errorlog,'a')
    f.write(errortext)
    f.close()

    errors = []
    print(errorlog)
    junk_nr = 1 
    for gene_aaclist_map in gene_aaclist_map_list:
        if len(gene_aaclist_map) == 0:
            continue

        print("Chunk %s/%s" % (str(junk_nr),str(len(gene_aaclist_map_list))))
        junk_nr+=1
        try:
            os.stat("%s/tmp_structman_pipeline" %(output_path))
        except:
            os.mkdir("%s/tmp_structman_pipeline" %(output_path))  
        os.chdir("%s/tmp_structman_pipeline" %(output_path))
        cwd = "%s/tmp_structman_pipeline" %(output_path)

        try:
            t1 = time.time()
            print("Before geneCheck")
            #check for already stored genes and make all the necessary database interactions
            #structure of stored_genes: {Uniprot-Id:gene_id}
            #structure of new_genes: {Uniprot-Id:gene_id}
            stored_genes,new_genes,stored_gene_ids,new_gene_ids = database.geneCheck(gene_aaclist_map,species_map,session,db,cursor)
            #print stored_genes
            #print new_genes

            t2 = time.time()
            if verbose:
                print("Time for geneCheck: %s" % (str(t2-t1)))

            print("Before mutationCheck", len(stored_genes), len(new_genes))
            #check for already stored mutations or position twins and make all the necessary database interactions
            #structure of gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
            #structure of stored_gene_new_pos: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
            gene_mut_map_new,stored_gene_new_pos,background_process_MS = database.mutationCheck(gene_aaclist_map,stored_genes,stored_gene_ids,new_genes,new_gene_ids,session,tag_map,db,cursor,config,residue_id_backmap)

            #print gene_mut_map_new

            t3 = time.time()
            if verbose:
                print("Time for mutationCheck: %s" % (str(t3-t2)))

            print("Before getSequences")
            #structure of gene_sequence_map: {Uniprot_Id:Sequence}
            gene_sequence_map,gene_mut_map_new,stored_gene_new_pos,pdb_pos_map,background_iu_process,No_Errors = getSequences(config,new_genes,stored_genes,fasta_map,species_map,gene_mut_map_new,
                                                                                                                                stored_gene_new_pos,corrected_input_pdbs,manager,lock,db,cursor)

            #print gene_mut_map_new

            t4 = time.time()
            if verbose:
                print("Time for getSequences: %s" % (str(t4-t3)))

            print("Before autoTemplateSelection")
            #structure of template_map: {Uniprot-Id:(template-list,{template-id:stored-template},oligo_map)}
            structure_map,gene_error_map,No_Errors = autoTemplateSelection(config,db,cursor,manager,lock,No_Errors,gene_sequence_map,new_genes,search_tool=search_tool)
            database.addErrorCodeToGene(gene_error_map,db,cursor)
            #print structure_map

            t5 = time.time()
            if verbose:
                print("Time for Template Selection: %s" % (str(t5-t4)))

            print("Before paraAlignment")
            #new structure of template_map: {Uniprot-Id:({template-id:new template},{template-id:stored-template},oligo_map)}
            structure_map,gene_template_alignment_map,structure_id_map,No_Errors,structure_dict,m_r_map,residue_dict,stored_structures = paraAlignment(config,db,cursor,manager,lock,No_Errors,gene_sequence_map,structure_map,stored_genes,stored_gene_ids,new_genes,gene_mut_map_new,stored_gene_new_pos,pdb_pos_map)

            t6 = time.time()
            if verbose:
                print("Time for Alignment: %s" % (str(t6-t5)))

            if background_iu_process != None:
                background_iu_process.join() #Disorder values have to finished before the classification happens

            #structure of gene_template_alignment_map: {gene_id:{template_id:(coverage,seq_id,sub_infos,alignment_pir)}}
            print("Before paraAnnotate")
            No_Errors = paraAnnotate(config,db,cursor,manager,lock,No_Errors,gene_mut_map_new,stored_gene_new_pos,gene_template_alignment_map,structure_map,structure_id_map,structure_dict,m_r_map,residue_dict,stored_structures)
            t7 = time.time()
            if verbose:
                print("Time for Annotation: %s" % (str(t7-t6)))

            t1 = time.time()
            #join the background inserts
            if background_process_MS != None:
                background_process_MS.join()
            #if background_process_AS != None:
            #    background_process_AS.join()


            t2 = time.time()
            if verbose:
                print('Resttime for background inserts: ',t2-t1)

        #Error-Handling for a whole input line
        except:
            
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc()
            #print "Pipeline Core Error: ",e,f,g
            errortext = '\n'.join([str(e),str(f),str(g)]) + '\n\n'
            f = open(errorlog,'a')
            f.write(errortext)
            f.close()
            No_Errors = False
            if background_process_MS != None:
                background_process_MS.join()

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

        print("At least one error occured, please check the errorlog.")

    if newsession:
        endtime = SQLDateTime()
        database.updateSession(session,endtime,db,cursor)
    db.close()

    tend = time.time()
    print((tend-t0))
    return session
