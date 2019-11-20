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

def sequenceScan(config,genes):
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
                    #print pos,aa
                    aac = '%s%s%s' % (aa,str(pos+1),aa)
                    genes[u_ac][2].add(aac)
                
            else: #The given residue-ids are mapped onto their sequence position
                new_aacs = set()
                for aac in genes[u_ac][2]:
                    res_id = aac[1:-1]
                    if not res_id in res_pos_map:
                        print('Skipped residue ',res_id,' from ',u_ac,' since illegal alternative location id')
                        continue
                    seq_pos = res_pos_map[res_id]
                    new_aacs.add('%s%s%s' % (aac[0],seq_pos,aac[-1]))
                genes[u_ac][2] = new_aacs

    if db != None:
        db.close()
    return genes,sequenceScanPDB,residue_id_backmap

def buildQueue(config,filename,junksize,mrna_fasta=None):
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

        if tag == None:
            tag = ''

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

    genes,corrected_input_pdbs,residue_id_backmap = sequenceScan(config,genes)

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

def getSequences(config,gene_aaclist_map,fasta_map,species_map,corrected_input_pdbs,manager,lock):
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
        u_acs = set()
        for u_ac in gene_aaclist_map:
            if len(u_ac) == 6 and u_ac.count(':') == 1:
                pdb_dict[u_ac] = gene_aaclist_map[u_ac]
            else:
                u_acs.add(u_ac)

        t0 = time.time()
        try:
            db_uniprot = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,'struct_man_db_uniprot')
            cursor_uniprot = db_uniprot.cursor()
        except:
            db_uniprot = None
            cursor_uniprot = None

        gene_sequence_map = uniprot.getSequencesPlain(u_acs,db_uniprot,cursor_uniprot)

        if db_uniprot != None:
            db_uniprot.close()

        
        t1 = time.time()
        if verbose:
            print("Time for getSequences Part 1: %s" % str(t1-t0))

        pdb_sequence_map,pdb_pos_map = pdb.getSequences(pdb_dict,pdb_path,AU=pdb_input_asymetric_unit)

        t2 = time.time()
        if verbose:
            print("Time for getSequences Part 2: %s" % str(t2-t1))

        gene_sequence_map.update(pdb_sequence_map)

        t3 = time.time()
        if verbose:
            print("Time for getSequences Part 3: %s" % str(t3-t2))

    else:
        gene_sequence_map = {}
        for species in fasta_map:
            mrna_fasta = fasta_map[species]
            f = open(mrna_fasta,'r')
            lines = f.readlines()
            f.close()

            gene_seq_map = {}
            for line in lines:
                line = line.replace('\n','')
                if line[0] == '>':
                    g_id = line.split()[0][1:]

                    gene_seq_map[g_id] = ''
                else:
                    gene_seq_map[g_id] += line

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
                    gene_sequence_map[gene] = sequence,None,None
                    gene_info_map[new_genes[gene]] = ([],{},{},sequence)


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

    return gene_sequence_map,iupred_map,No_Errors

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

def autoTemplateSelection(config,manager,lock,No_Errors,gene_sequence_map,search_tool='MMseqs2'):
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

        for gene in gene_sequence_map:
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
     
        t1 = time.time()

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

    elif search_tool == 'MMseqs2':
        t0 = time.time()
        gene_error_map = {}

        if len(gene_sequence_map) > 0:
            mmseqs_tmp_folder = config.mmseqs_tmp_folder
            raw_structure_map,pdb_ids = MMseqs2.search(gene_sequence_map,mmseqs2_db_path,mmseqs2_path,mmseqs_tmp_folder,option_seq_thresh,verbose=verbose)

            t1 = time.time()

            structure_map = templateSelection.filterRawStructureMap(raw_structure_map,pdb_ids,pdb_path,option_res_thresh,blast_processes)
            t2 = time.time()
            if verbose:
                print("Template Selection Part 1: %s" % (str(t1-t0)))
                print("Template Selection Part 2: %s" % (str(t2-t1)))
        else:
            structure_map = {}

    return structure_map,gene_error_map,No_Errors

def paraAlignment(config,manager,lock,No_Errors,gene_sequence_map,structure_map,gene_aaclist_map):
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

    n = 0
    for gene in gene_aaclist_map:
        
        aaclist = gene_aaclist_map[gene][2]
        aacbases = {}
        for aac in aaclist:
            aacbase = aac[:-1]
            aa2 = aac[-1]
            if not aacbase in aacbases:
                aacbases[aacbase] = set(aa2)
            else:
                aacbases[aacbase].add(aa2)

        gene_aaclist_map[gene][2] = aacbases

        if not gene in structure_map:
            continue

        structures = structure_map[gene]

        for (pdb_id,chain) in structures:
            input_queue.put((gene,gene_sequence_map[gene][0],pdb_id,chain,structures[(pdb_id,chain)],aacbases))

        if n == 1000:
            input_queue = manager.Queue()
            in_queues.append(input_queue)
            n = 0
        n += 1

    t1 = time.time()
    if verbose:
        print("Alignment Part 1: %s" % (str(t1-t0)))

    t2 = time.time()
    if verbose:
        print("Alignment Part 2: %s" % (str(t2-t1)))

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

            if gene in gene_aaclist_map:
                gene_aaclist_map[gene][2] = aaclist

            if not gene in gene_template_alignment_map:
                gene_template_alignment_map[gene] = {}
            gene_template_alignment_map[gene][(pdb_id,chain)] = sub_infos

        del input_queue
        del out_queue
        del error_queue

        t5 += time.time()

    t6 += time.time()

    t7 += time.time()

    t8 += time.time()

    t9 += time.time()

    if verbose:
        print("Alignment Part 4: %s" % (str(t4-t31)))
        print("Alignment Part 5: %s" % (str(t5-t4)))
        print("Alignment Part 6: %s" % (str(t6-t5)))
        print("Alignment Part 7: %s" % (str(t7-t6)))
        print("Alignment Part 8: %s" % (str(t8-t7)))
        print("Alignment Part 9: %s" % (str(t9-t8)))

    return structure_map,gene_template_alignment_map,No_Errors


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

                (coverage,seq_id,sub_infos,alignment_pir,errors,times,aaclist,update_map) = globalAlignment.alignBioPython(gene,wildtype_sequence,pdb_id,template_page,chain,aaclist,ignore_gaps=True)
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
    

def paraAnnotate(config,manager,lock,No_Errors,gene_aaclist_map,gene_template_alignment_map,structure_map,residue_id_backmap,tag_map,species_map,outpath,session_name,iupred_map):
    annotation_processes = config.annotation_processes
    anno_session_mapping = config.anno_session_mapping
    error_annotations_into_db = config.error_annotations_into_db

    smiles_path = config.smiles_path
    inchi_path = config.inchi_path
    pdb_path = config.pdb_path
    verbose = config.verbose
    errorlog = config.errorlog_path

    t0 = time.time()
    
    input_queue = manager.Queue()

    pdb_structure_map = {}
    size_sorted = {}
    max_size = 0
    gene_structure_map = {}
    pdb_ids = set()
    for gene in structure_map:
        structures = structure_map[gene]
        for (pdb_id,chain) in structures:
            pdb_ids.add(pdb_id)
            structure = structures[(pdb_id,chain)]
            cov = structure['Coverage']
            seq_id = structure['Seq_Id']

            if 'IAP' not in structure:
                #del structure_map[gene][(pdb_id,chain)]
                continue

            sub_infos = gene_template_alignment_map[gene][(pdb_id,chain)]
            if len(sub_infos) == 0:
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
            
            if not chain in pdb_structure_map[pdb_id]:
                pdb_structure_map[pdb_id][chain] = [structure,sub_infos]
            else:
                pdb_structure_map[pdb_id][chain][1].update(sub_infos)


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
        (annotation_chain_dict,pdb_id,chain_structure_map,residue_residue_dict,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles) = out
        
        for chain in annotation_chain_dict:
            annotations = annotation_chain_dict[chain]

            total_annotations[(pdb_id,chain)] = annotations,residue_residue_dict[chain]

        complex_profiles[pdb_id] = ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles

    min_l_dict,min_m_dict,min_i_dict,min_c_dict,min_r_dict,min_d_dict,min_homo_dict,mutation_surface_dict,mutation_sec_dict,mutation_inter_dict,structure_classification_map,mutation_oligo_chain_map,centrality_dict,mutation_dict,aa1_dict = createStructureDicts(manager,lock,gene_aaclist_map,gene_template_alignment_map,total_annotations,iupred_map,pdb_structure_map,ligand_filter=config.ligand_filter,n_of_processes=config.proc_n)

    class_dict,structure_classification_map,interaction_dict = database.createClassDict(min_l_dict,min_m_dict,min_i_dict,min_c_dict,min_r_dict,min_d_dict,min_homo_dict,mutation_surface_dict,mutation_sec_dict,centrality_dict,mutation_dict,structure_classification_map,complex_profiles,config.verbose)


    outfile = '%s/%s' % (outpath,session_name)
    outfile = "%s.classification.tsv" % (outfile)
    appendOutput(class_dict,structure_classification_map,interaction_dict,gene_aaclist_map,residue_id_backmap,tag_map,species_map,outfile,aa1_dict,mutation_sec_dict)

    return No_Errors

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
            annotation_chain_dict,residue_residue_dict,errorlist,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles = templateFiltering.liteAnalysis(pdb_id,chain_structure_map,pdb_path,dssp_path,rin_db_path,neighborhood_calculation=neighborhood_calculation,dssp=dssp,calculate_interaction_profiles=calculate_interaction_profiles)
            with lock:
                out_queue.put((annotation_chain_dict,pdb_id,chain_structure_map,residue_residue_dict,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles))
                if len(errorlist) > 0:
                    for (error,e,f,g) in errorlist:
                        err_queue.put((error,f,g,pdb_id))
            #print 'finished Proc %s - annotate: %s' % (str(i),pdb_id)
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc()
            with lock:
                err_queue.put((e,f,g,pdb_id))

def appendOutput(class_dict,structure_classification_map,interaction_dict,gene_aaclist_map,residue_id_backmap,tag_map,species_map,outfile,aa1_dict,mutation_sec_dict):
    lines = []

    for u_ac in gene_aaclist_map:
        u_id,refseqs,aaclist = gene_aaclist_map[u_ac]
        gpan = ';'.join(refseqs)
        for aac_base in aaclist:
            aa2 = aaclist[aac_base]
            pos = aac_base[1:]
            m = (u_ac,pos)
            aa1 = aa1_dict[m]
            if not m in class_dict:
                continue
            aac = '%s%s' % (aac_base,aa2)

            if aac in tag_map[u_ac]:
                tag = tag_map[u_ac][aac]
            else:
                tag = 'None'

            (Class,conf,weighted_sc,conf_sc,best_res,max_seq_res,amount_of_structures,
            weighted_c,conf_c,
            weighted_d,conf_d,
            weighted_r,conf_r,
            weighted_l,conf_l,
            weighted_m,conf_m,
            weighted_i,conf_i,
            weighted_h,conf_h,
            max_seq_id,
            weighted_raw,weighted_cent,weighted_norm,
            weighted_lig_degree,weighted_lig_score,
            weighted_metal_degree,weighted_metal_score,
            weighted_ion_degree,weighted_ion_score,
            weighted_prot_degree,weighted_prot_score,
            weighted_rna_degree,weighted_rna_score,
            weighted_dna_degree,weighted_dna_score,
            weighted_modres,modres_prop) = class_dict[m]

            if m in mutation_sec_dict:
                mv_sec_ass = database.majority_vote(mutation_sec_dict[m])
            else:
                mv_sec_ass = None

            simple_class = database.simplifyClass(Class,weighted_sc)

            if best_res != None:
                [r_id,qual,res_aa,res_nr,pdb_id,chain,resolution,cov,seq_id,rsa,min_lig,min_metal,min_ion,iacs] = best_res

                recommended_structure = '%s:%s %s:%s;%s;%s;%s' % (pdb_id,chain,res_nr,res_aa,str(seq_id),str(cov),str(resolution))
            else:
                resolution = '-'
                cov = '-'
                seq_id = '-'
                recommended_structure = '-'

            if max_seq_res != None:
                [max_seq_r_id,max_seq_qual,max_seq_res_aa,max_seq_res_nr,max_seq_pdb_id,max_seq_chain,max_seq_resolution,max_seq_cov,max_seq_seq_id,max_seq_rsa,max_min_lig,max_min_metal,max_min_ion,max_iacs] = max_seq_res

                max_seq_structure = '%s:%s %s:%s;%s;%s;%s' % (max_seq_pdb_id,max_seq_chain,max_seq_res_nr,max_seq_res_aa,str(max_seq_seq_id),str(max_seq_cov),str(max_seq_resolution))
            else:
                max_seq_resolution = '-'
                max_seq_cov = '-'
                max_seq_seq_id = '-'
            max_seq_structure = '-'

            input_res_id = ''
            input_pdb_id = ''
            if u_ac in residue_id_backmap:
                input_res_id = residue_id_backmap[u_ac][pos]
                input_pdb_id = u_ac
                u_ac = ''
                u_id = ''
                gpan = ''

            if m in interaction_dict:
                interaction_str = str(interaction_dict[m])
            else:
                interaction_str = None

            lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (u_ac,u_id,gpan,input_pdb_id,input_res_id,aa1,pos,species_map[u_ac],tag,weighted_sc,Class,simple_class,interaction_str,str(conf),mv_sec_ass,recommended_structure,seq_id,cov,resolution,max_seq_structure,max_seq_seq_id,max_seq_cov,max_seq_resolution,str(amount_of_structures)))

    f = open(outfile,'a')
    f.write('\n'.join(lines) + '\n')
    f.close()
    return

def paraStructureDicts(inqueue,outqueue,lock,gene_template_alignment_map,annotations,pdb_structure_map,ligand_filter):
    residue_calc = {}

    with lock:
        inqueue.put(None)

    outs = []
    outs_package_size = 100000
    n = 0

    while True:
        with lock:
            m_ids = inqueue.get()
        if m_ids == None:
            break
        for m_id in m_ids:
            (u_ac,aacbase) = m_id
            m_id = (u_ac,aacbase[1:])
            if not u_ac in gene_template_alignment_map:
                continue
            for (pdb_id,chain) in gene_template_alignment_map[u_ac]:
                sub_infos = gene_template_alignment_map[u_ac][(pdb_id,chain)]
                if not aacbase in sub_infos:
                    continue
                res_id = sub_infos[aacbase][0]
                res_aa,lig_dists,chain_dists,raac,ssa,homomer_map,profile_str,centrality_scores,avg_b_factor,modres = annotations[(pdb_id,chain)][0][res_id]
                
                lig_dist_strings = []
                dists = []
                for lig_id in lig_dists:
                    (dist,atom_pair) = lig_dists[lig_id]
                    if not dist == None:
                        dists.append(dist)
                        lig_dist_strings.append("%s:%s:(%s-%s)" % (str(lig_id),dist,str(atom_pair[0]),str(atom_pair[1])))
                if len(dists) > 0:
                    min_lig_dist = min(dists)
                else:
                    min_lig_dist = -1
                dists = []
                lig_dist_string = ",".join(lig_dist_strings)
                chain_dist_strings = []
                for chain_id in chain_dists:
                    (dist,atom_pair,min_resi) = chain_dists[chain_id]
                    if not dist == None:
                        dists.append(dist)
                        chain_dist_strings.append("%s.%s:%s:(%s-%s)" % (chain_id,min_resi,dist,str(atom_pair[0]),str(atom_pair[1])))
                if len(dists) > 0:
                    min_chain_dist = min(dists)
                else:
                    min_chain_dist = -1
                chain_dist_string = ",".join(chain_dist_strings)

                if len(lig_dists) == 0:
                    lig_dist_string = "-"
                if len(chain_dists) == 0:
                    chain_dist_string = "-"

                n += 1

                if centrality_scores != None and centrality_scores != 'NULL':
                    Raw_Centrality,Centrality,Norm_Centrality,Cleaned_Raw_Centrality,Cleaned_Centrality,Cleaned_Norm_Centrality = centrality_scores
                else:
                    Raw_Centrality = None
                    Centrality = None
                    Norm_Centrality = None
                    Cleaned_Raw_Centrality = None
                    Cleaned_Centrality = None
                    Cleaned_Norm_Centrality = None
                structure = pdb_structure_map[pdb_id][chain][0]
                
                resolution = structure['Resolution']
                cov = structure['Coverage']
                seq_id = structure['Seq_Id']

                oligos = structure['Oligo']

                chains = []

                for iap in structure['IAP']:
                    if iap[0] == 'Ligand':
                        continue
                    chains.append('%s:%s' % (iap[1],iap[0]))
                    
                chains = ','.join(chains)

                qual = templateFiltering.qualityScore(resolution,cov,seq_id)

                if not (pdb_id,chain,res_id) in residue_calc:
                    try:
                        rel_sur_acc = raac
                        if rel_sur_acc > 1.0:
                            rel_sur_acc = 0.01*rel_sur_acc
                    except:
                        rel_sur_acc = None
                    try:
                        sec_str_ass = ssa
                    except:
                        sec_str_ass = None
                    
                    min_ld,min_md,min_id,min_cd,min_rd,min_dd,min_lig,min_metal,min_ion,iacs = database.processClassData(lig_dist_string,chain_dist_string,chains,ligand_filter=ligand_filter)

                    minimal_distances = []
                    if min_cd != None:
                        minimal_distances.append(min_cd)
                    if min_dd != None:
                        minimal_distances.append(min_dd)
                    if min_rd != None:
                        minimal_distances.append(min_rd)
                    if min_ld != None:
                        minimal_distances.append(min_ld)
                    if min_md != None:
                        minimal_distances.append(min_md)
                    if min_id != None:
                        minimal_distances.append(min_id)

                    if len(minimal_distances) == 0:
                        min_minimal_distances = 2.0
                    else:
                        min_minimal_distances = min(minimal_distances)

                    if homomer_map == {}:
                        homomer_dist = None
                    else:
                        h_dists = []
                        for h_chain in homomer_map:
                            h_dist = homomer_map[h_chain]
                            h_dists.append(h_dist)
                        homomer_dist = min(h_dists)

                    residue_calc[(pdb_id,chain,res_id)] = (min_minimal_distances,rel_sur_acc,sec_str_ass,homomer_dist,min_ld,min_md,min_id,min_cd,min_rd,min_dd,
                                        profile_str,Raw_Centrality,Centrality,Norm_Centrality,min_lig,min_metal,min_ion,iacs)
                    

                else:
                    (min_minimal_distances,rel_sur_acc,sec_str_ass,homomer_dist,min_ld,min_md,min_id,min_cd,min_rd,min_dd,
                                        profile_str,Raw_Centrality,Centrality,Norm_Centrality,min_lig,min_metal,min_ion,iacs) = residue_calc[(pdb_id,chain,res_id)]

                if rel_sur_acc > database.surface_threshold:
                    sc = "Surface"
                else:
                    sc = "Core"
                Class,conf = database.getWeightedClass(sc,1.0,min_cd,1.0,min_dd,1.0,min_rd,1.0,min_ld,1.0,min_md,1.0,min_id,1.0)
                simpleClass =  database.simplifyClass(Class,sc)

                out = (m_id,(pdb_id,chain,res_id),Class,simpleClass,qual,res_aa,res_id,pdb_id,chain,resolution,cov,seq_id,
                                profile_str,
                                Raw_Centrality,Centrality,Norm_Centrality,
                                oligos,chains,
                                rel_sur_acc,
                                sec_str_ass,
                                homomer_dist,
                                min_ld,
                                min_md,
                                min_id,
                                min_cd,
                                min_rd,
                                min_dd,
                                min_minimal_distances,
                                min_lig,min_metal,min_ion,iacs,modres
                        )
                outs.append(out)
                if n == outs_package_size:
                    with lock:
                        outqueue.put(outs)
                    outs = []
                    n = 0
    if outs != []:
        with lock:
            outqueue.put(outs)
    

def createStructureDicts(manager,lock,gene_aaclist_map,gene_template_alignment_map,annotations,iupred_map,pdb_structure_map,ligand_filter=None,n_of_processes=1):
    min_l_dict = {}
    min_m_dict = {}
    min_i_dict = {}
    min_c_dict = {}
    min_r_dict = {}
    min_d_dict = {}
    min_homo_dict = {}

    mutation_surface_dict = {}

    mutation_sec_dict = {}

    mutation_inter_dict = {}

    mutation_oligo_chain_map = {}

    centrality_dict = {}

    residue_calc = {}

    aa1_dict = {}

    #needed in order to pick the recommended structure
    structure_classification_map = {}

    inqueue = manager.Queue()
    outqueue = manager.Queue()

    package_size = 100
    package = []

    i = 0
    n = 0

    mutation_dict = {}

    for u_ac in gene_aaclist_map:
        regions = iupred_map[u_ac][0]
        method = iupred_map[u_ac][2]
        if method == 'MobiDB3.0'or method == 'mobidb-lite':
            pos_region_type = 'globular'
        else:
            pos_region_type = 'disorder'
        for aacbase in gene_aaclist_map[u_ac][2]:
            i += 1

            aac_pos = int(aacbase[1:])
            for [a,b,region_type] in regions:
                if aac_pos > int(a) and aac_pos < int(b):
                    pos_region_type = region_type
            if aac_pos in iupred_map[u_ac][1]:
                aa1,iupred_score = iupred_map[u_ac][1][aac_pos]
            else:
                iupred_score = None

            mutation_dict[(u_ac,aacbase[1:])] = (aacbase,None,iupred_score,pos_region_type,None)
            structure_classification_map[(u_ac,aacbase[1:])] = {}
            mutation_inter_dict[(u_ac,aacbase[1:])] = []
            mutation_oligo_chain_map[(u_ac,aacbase[1:])] = [0,0]
            mutation_surface_dict[(u_ac,aacbase[1:])] = []
            mutation_sec_dict[(u_ac,aacbase[1:])] = []
            centrality_dict[(u_ac,aacbase[1:])] = []
            min_l_dict[(u_ac,aacbase[1:])] = []
            min_m_dict[(u_ac,aacbase[1:])] = []
            min_i_dict[(u_ac,aacbase[1:])] = []
            min_c_dict[(u_ac,aacbase[1:])] = []
            min_r_dict[(u_ac,aacbase[1:])] = []
            min_d_dict[(u_ac,aacbase[1:])] = []
            min_homo_dict[(u_ac,aacbase[1:])] = []

            aa1_dict[(u_ac,aacbase[1:])] = aacbase[0]

            package.append((u_ac,aacbase))
            if i == package_size:
                inqueue.put(package)
                package = []
                i = 0
                n += 1

    if package != []:
        inqueue.put(package)
        n += 1

    print('Amount of packages: ',n)

    processes = {}
    for i in range(1,min(n+1,n_of_processes + 1)):
        p = multiprocessing.Process(target=paraStructureDicts, args=(inqueue,outqueue,lock,gene_template_alignment_map,annotations,pdb_structure_map,ligand_filter))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()
    with lock:
        outqueue.put(None)

    while True:
        outs = outqueue.get()
        if outs == None:
            break
        for out in outs:
            (m_id,r_id,Class,simpleClass,qual,res_aa,residue,pdb_id,chain,resolution,cov,seq_id,
                                profile_str,
                                Raw_Centrality,Centrality,Norm_Centrality,
                                oligos,chains,
                                rel_sur_acc,
                                sec_str_ass,
                                homomer_dist,
                                min_ld,
                                min_md,
                                min_id,
                                min_cd,
                                min_rd,
                                min_dd,
                                min_minimal_distances,
                                min_lig,min_metal,min_ion,iacs,modres
                        ) = out

            structure_classification_map[m_id][r_id] = (Class,simpleClass,qual,res_aa,residue,pdb_id,chain,resolution,cov,seq_id,rel_sur_acc,min_lig,min_metal,min_ion,iacs,modres)
            if profile_str != 'NULL':
                mutation_inter_dict[m_id].append((profile_str,qual))
            if Centrality != None and Centrality != 'NULL':
                centrality_dict[m_id].append((Raw_Centrality,Centrality,Norm_Centrality,qual)) 

            if len(oligos) > mutation_oligo_chain_map[m_id][0]:
                mutation_oligo_chain_map[m_id][0] = len(oligos)
            if len(chains.split(',')) > mutation_oligo_chain_map[m_id][1]:
                mutation_oligo_chain_map[m_id][1] = len(chains.split(','))

            if not min_minimal_distances < 1.2:
                #print rel_sur_acc
                if rel_sur_acc != None and cov > 0.0:
                    mutation_surface_dict[m_id].append((rel_sur_acc,qual,cov))

                mutation_sec_dict[m_id].append((sec_str_ass,qual))

                if not min_ld == None:
                    min_l_dict[m_id].append((min_ld,qual))
                if not min_md == None:
                    min_m_dict[m_id].append((min_md,qual))
                if not min_id == None:
                    min_i_dict[m_id].append((min_id,qual))
                if not min_cd == None:
                    min_c_dict[m_id].append((min_cd,qual))
                if not min_rd == None:
                    min_r_dict[m_id].append((min_rd,qual))
                if not min_dd == None:
                    min_d_dict[m_id].append((min_dd,qual))

                if not homomer_dist == None:
                    min_homo_dict[m_id].append((homomer_dist,qual))

    return min_l_dict,min_m_dict,min_i_dict,min_c_dict,min_r_dict,min_d_dict,min_homo_dict,mutation_surface_dict,mutation_sec_dict,mutation_inter_dict,structure_classification_map,mutation_oligo_chain_map,centrality_dict,mutation_dict,aa1_dict

def main(filename,config,output_path,main_file_path):
    n_of_cores = config.proc_n
    species = config.species
    mrna_fasta = config.mrna_fasta
    num_of_cores = config.proc_n
    verbose = config.verbose
    search_tool = config.search_tool
    errorlog = config.errorlog_path

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

    session_name = (nfname.rsplit("/",1)[1]).rsplit(".",1)[0]

    outfile = '%s/%s' % (output_path,session_name)
    outfile = "%s.classification.tsv" % (outfile)

    if verbose:
        print('Writing outputs into: ',outfile)

    f = open(outfile,'w')
    f.write("Uniprot-Ac\tUniprot Id\tRefseq\tPDB-ID (Input)\tResidue-Id\tAmino Acid\tPosition\tSpecies\tTag\tWeighted Surface/Core\tClass\tSimple Class\tIndividual Interactions\tConfidence Value\tSecondary Structure\tRecommended Structure\tSequence-ID\tCoverage\tResolution\tMax Seq Id Structure\tMax Sequence-ID\tMax Seq Id Coverage\tMax Seq Id Resolution\tAmount of mapped structures\n")
    f.close()

    mrna_fasta_for_annovar = True #make this a option, for the case of mrna fasta files with non-standard protein-ids (may include further debugging)

    if mrna_fasta_for_annovar:
        mrna_fasta = None

    t01 = time.time()

    if verbose:
        print("Time for preparation before buildQueue: %s" % (str(t01-t0)))

    junksize = 500

    if verbose:
        print("Call buildQueue with chunksize: %s and file: %s" % (str(junksize),nfname))
    gene_aaclist_map_list,tag_map,species_map,fasta_map,corrected_input_pdbs,residue_id_backmap = buildQueue(config,nfname,junksize,mrna_fasta=mrna_fasta)

    t02 = time.time()
    if verbose:
        print("Time for buildQueue: %s" % (str(t02-t01)))
    
    print("Number of chunks: ",len(gene_aaclist_map_list))

    date = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
    errortext  = "###############################################################################\n%s:\n%s\n" % (date,nfname)

    f = open(errorlog,'a')
    f.write(errortext)
    f.close()

    errors = []
    print(errorlog)
    junk_nr = 1 
    for gene_aaclist_map in gene_aaclist_map_list:
        if len(gene_aaclist_map) == 0:
            #print "empty gene_aaclist_map"
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
            t3 = time.time()
            print("Before getSequences")

            gene_sequence_map,iupred_map,No_Errors = getSequences(config,gene_aaclist_map,fasta_map,species_map,corrected_input_pdbs,manager,lock)

            t4 = time.time()
            if verbose:
                print("Time for getSequences: %s" % (str(t4-t3)))

            print("Before autoTemplateSelection")
            structure_map,gene_error_map,No_Errors = autoTemplateSelection(config,manager,lock,No_Errors,gene_sequence_map,search_tool=search_tool)

            t5 = time.time()
            if verbose:
                print("Time for Template Selection: %s" % (str(t5-t4)))

            print("Before paraAlignment")
            structure_map,gene_template_alignment_map,No_Errors = paraAlignment(config,manager,lock,No_Errors,gene_sequence_map,structure_map,gene_aaclist_map)

            t6 = time.time()
            if verbose:
                print("Time for Alignment: %s" % (str(t6-t5)))

            #structure of gene_template_alignment_map: {gene_id:{template_id:(coverage,seq_id,sub_infos,alignment_pir)}}
            print("Before paraAnnotate")
            No_Errors = paraAnnotate(config,manager,lock,No_Errors,gene_aaclist_map,gene_template_alignment_map,structure_map,residue_id_backmap,tag_map,species_map,output_path,session_name,iupred_map)
            t7 = time.time()
            if verbose:
                print("Time for Annotation: %s" % (str(t7-t6)))


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

    
    tend = time.time()
    print((tend-t0))
    return
