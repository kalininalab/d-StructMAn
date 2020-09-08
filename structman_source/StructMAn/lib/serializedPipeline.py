import pymysql as MySQLdb
import os
import sys
import getopt
import time
import multiprocessing
from multiprocessing.managers import SyncManager
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
import sdsc

import cProfile
import pstats
import io
try:
    from memory_profiler import profile
except:
    pass

def SQLDateTime():
    return time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())

#@profile
def sequenceScan(config,proteins):
    pdb_path = config.pdb_path

    sequenceScanProteins = {}
    sequenceScanPDB = {}

    for u_ac in proteins:
        uni_pos,tags = proteins[u_ac].popNone()
        
        if u_ac.count(':') > 0:
            sequenceScanPDB[u_ac] = tags,uni_pos #PDB inputs are always processed by the sequence scan, but the positions are only added if uni_pos is true
        else:
            sequenceScanProteins[u_ac] = tags,uni_pos #New: process everything to filter input by sanity checks


    if len(sequenceScanProteins) > 0:
        if config.verbosity >= 1:
            print("Amount of proteins going into sequenceScan: ",len(sequenceScanProteins))

        gene_sequence_map = uniprot.getSequencesPlain(sequenceScanProteins.keys(),config,debug=config.verbose)
        for u_ac in gene_sequence_map:
            if gene_sequence_map[u_ac][0] == 1 or gene_sequence_map[u_ac][0] == 0:
                config.errorlog.add_error("Error in sequenceScan with gene: %s" % u_ac)
                continue
            seq,disorder_scores,disorder_regions_datastruct = gene_sequence_map[u_ac]
            proteins[u_ac].sequence = seq
            proteins[u_ac].set_disorder_scores(disorder_scores)
            proteins[u_ac].set_disorder_regions(disorder_regions_datastruct)
            tags,uni_pos = sequenceScanProteins[u_ac]

            for (pos,aa) in enumerate(seq):
                seq_pos = pos+1
                if seq_pos in proteins[u_ac].positions:
                    proteins[u_ac].positions[seq_pos].check(aa,overwrite = config.overwrite_incorrect_wt_aa)
                    proteins[u_ac].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    position = sdsc.Position(pos=seq_pos,wt_aa=aa,tags=tags,checked=True)
                    proteins[u_ac].positions[seq_pos] = position
            #sanity check filter
            del_list = []
            for pos in proteins[u_ac].positions:
                if not proteins[u_ac].positions[pos].checked:
                    del_list.append(pos)
            for pos in del_list:
                if config.verbosity >= 2:
                    print('Filtered position through sanity check:',u_ac,pos)
                del proteins[u_ac].positions[pos]


    if len(sequenceScanPDB) > 0:
        pdb_sequence_map,pdb_pos_map = pdb.getSequences(sequenceScanPDB.keys(),pdb_path)
        for u_ac in pdb_sequence_map:
            proteins[u_ac].sequence = pdb_sequence_map[u_ac][0]
            tags,uni_pos = sequenceScanPDB[u_ac]

            residue_id_backmap = {}
            res_pos_map = pdb_pos_map[u_ac]
            for res_id in res_pos_map:
                residue_id_backmap[res_pos_map[res_id]] = res_id

            for (pos,aa) in enumerate(pdb_sequence_map[u_ac][0]):
                seq_pos = pos+1
                if not seq_pos in residue_id_backmap:
                    continue
                res_id = residue_id_backmap[seq_pos]
                if res_id in proteins[u_ac].res_id_map:
                    proteins[u_ac].res_id_map[res_id].pos = seq_pos
                    proteins[u_ac].positions[seq_pos] = proteins[u_ac].res_id_map[res_id]
                    if tags == None:
                        continue
                    proteins[u_ac].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    if tags == None:
                        continue
                    position = sdsc.Position(pos=seq_pos,pdb_res_nr=res_id,wt_aa=aa,tags=tags)
                    proteins[u_ac].positions[seq_pos] = position
                    proteins[u_ac].res_id_map[res_id] = position

    return proteins

def parseFasta(config,nfname):
    f = open(nfname,'r')
    lines = f.readlines()
    f.close()

    seq_map = {}

    for line in lines:
        line = line[:-1]
        if len(line) == 0:
            continue
        if line[0] == '>':
            entry_id = line[1:].split()[0]
            seq_map[entry_id] = ''

        else:
            seq_map[entry_id] += line

    proteins = {}
    for prot_id in seq_map:
        positions = set()
        protein = sdsc.Protein(u_ac=prot_id,positions = positions)
        seq = seq_map[prot_id]
        protein.sequence = seq
        for (pos,aa) in enumerate(seq):
            seq_pos = pos+1
            if  not seq_pos in protein.positions:
                position = sdsc.Position(pos=seq_pos,wt_aa=aa)
                protein.positions[seq_pos] = position

        proteins[prot_id] = protein
    return [proteins]

#@profile
def buildQueue(config,filename,junksize,mrna_fasta=None):
    verbose = config.verbose

    t0 =time.time()

    proteins = {}
    u_ids = set()
    u_acs = set()
    id_map = {}
    ac_map = {}
    np_map = {}

    pdb_map = {}

    f = open(filename, "r")
    lines = f.read().split('\n')
    f.close()
    ac_id_map = {}
    tag_map = {}
    fasta_map = {}

    for line in lines:
        if line == '':
            continue
        if len(line) < 3:
            if config.verbosity >= 1:
                print("Skipped input line:\n%s\nToo short.\n" % line)
            continue
        line = line.replace(' ','\t')
        words = line.split("\t")
        if len(words) < 1:
            if config.verbosity >= 1:
                print("Skipped input line:\n%s\nToo few words.\n" % line)
            continue
        sp_id = words[0]#.replace("'","\\'")

        tags = set()

        if len(words) > 2:
            tags = set(words[2].split(','))

        try:
            if len(words) == 1 or words[1] == '':
                position = sdsc.Position(tags=tags)

            else:

                aachange = words[1].replace("\n","")
                if not sp_id.count(':') == 1: #this means sp_id is not a pdb-id

                    if ord(aachange[0]) > 47 and ord(aachange[0]) < 58: #if first char is a number
                        aa1 = 'X'
                        aachange = "X%s" % aachange
                    else:
                        aa1 = aachange[0]


                    if ord(aachange[-1]) > 47 and ord(aachange[-1]) < 58: #if last char is a number

                        pos = int(aachange[1:])


                        position = sdsc.Position(pos=pos,wt_aa=aa1,tags=tags)
                    else:
                        aa2 = aachange.split(',')[0][-1]
                        pos = int(aachange.split(',')[0][1:-1])

                        position = sdsc.Position(pos=pos,wt_aa=aa1,mut_aas = set(aa2),tags=tags)

                else: #this means sp_id is a pdb-id
                    if aachange.count('_') == 1: #This means a mutant amino acid type is given
                        aachange = aachange.replace('_','')
                        position = sdsc.Position(pdb_res_nr=aachange[1:-1],wt_aa=aachange[0],mut_aas = set(aachange[-1]),tags=tags)
                    else:
                        position = sdsc.Position(pdb_res_nr=aachange[1:],wt_aa=aachange[0],tags=tags)

        except:
            print("File Format Error: ",line)
            sys.exit()

        if mrna_fasta == None:

            if sp_id[2] == "_":
                if not sp_id in np_map:
                    np_map[sp_id] = [position]
                else:
                    np_map[sp_id].append(position)
                
            else:
                if sp_id.count('_') > 0:
                    u_ids.add(sp_id)
                    if not sp_id in id_map:
                        id_map[sp_id] = [position]
                    else:
                        id_map[sp_id].append(position)
                elif len(sp_id) == 6 and sp_id.count(':') == 1:
                    pdb_chain_tuple = '%s:%s' % (sp_id[:4].upper(),sp_id[-1]) #enforce uppercase pdb-id 
                    if not pdb_chain_tuple in pdb_map:
                        pdb_map[pdb_chain_tuple] = [position]
                    else:
                        pdb_map[pdb_chain_tuple].append(position)
                else:
                    u_acs.add(sp_id)
                    if not sp_id in ac_map:
                        ac_map[sp_id] = [position]
                    else:
                        ac_map[sp_id].append(position)
        else:
            #fasta input needs an update
            """
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
            """

    t1 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 1: ",str(t1-t0))

    try:
        db_adress = config.db_adress
        db_user_name = config.db_user_name
        db_password = config.db_password
        db_uniprot = MySQLdb.connect(db_adress,db_user_name,db_password,config.mapping_db)
        cursor_uniprot = db_uniprot.cursor()
    except:
        db_uniprot = None
        cursor_uniprot = None

    proteins = uniprot.IdMapping(ac_map,id_map,np_map,db_uniprot,cursor_uniprot,pdb_map,verbose=verbose)

    if db_uniprot != None:
        db_uniprot.close()

    t2 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 2: ",str(t2-t1))

    proteins = sequenceScan(config,proteins)

    t3 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 3: ",str(t3-t2))

    t4 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 4: ",str(t4-t3))


    outlist = []

    s = len(proteins)
    if config.verbosity >= 1:
        print("Total proteins: ",s)

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
                (key,value) = proteins.popitem()
                new_map[key] = value
            outlist.append(new_map)
        new_map = {}
        while len(proteins) > 0:
            (key,value) = proteins.popitem()
            new_map[key] = value
        if len(new_map) > 0:
            outlist.append(new_map)
    else:
        outlist.append(proteins)

    t5 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 5: ",str(t5-t4))

    return outlist

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

#@profile
def getSequences(proteins,config,manager,lock,skip_db = False):

    number_of_processes = config.proc_n
    pdb_path = config.pdb_path
    pdb_input_asymetric_unit = config.pdb_input_asymetric_unit
    blast_processes = config.blast_processes
    cwd = os.getcwd()
    verbose = config.verbose

    t0 = time.time()

    uniprot.getSequences(proteins,config)

    t1 = time.time()
    if config.verbosity >= 2:
        print("Time for getSequences Part 1: %s" % str(t1-t0))

    #pdb_sequence_map,pdb_pos_map = pdb.getSequences(pdb_dict,pdb_path,AU=pdb_input_asymetric_unit)

    if not skip_db:
        db,cursor = config.getDB()
        database.addProtInfos(proteins,db,cursor)
        db.close()

    t2 = time.time()
    if config.verbosity >= 2:
        print("Time for getSequences Part 2: %s" % str(t2-t1))

    u_acs = proteins.get_protein_u_acs()

    background_iu_process = None
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

        n_disorder = 0
        stored_disorder_ids = []
        for u_ac in u_acs:
            seq = proteins.get_sequence(u_ac)
            if seq == 0 or seq == 1 or seq == None:
                continue
            disorder_scores = proteins.get_disorder_scores(u_ac)
            disorder_regions = proteins.get_disorder_regions(u_ac)
            if disorder_scores == None:
                if not mobi_lite:
                    in_queue.put((u_ac,seq))
                    n_disorder += 1
                else:
                    mobi_list.append('>%s\n%s\n' % (u_ac,seq))
            elif disorder_scores != 'Stored':
                disorder_scores_datastruct = {}
                for pos,score in enumerate(disorder_scores):
                    seq_pos = pos + 1
                    if pos >= len(seq):
                        config.errorlog.add_warning('Warning: illegal position %s for sequence of %s' % (str(pos),u_ac))
                        continue
                    disorder_scores_datastruct[seq_pos] = score

                proteins.set_disorder_scores(u_ac,disorder_scores_datastruct)
                proteins.set_disorder_regions(u_ac,disorder_regions)
                proteins.set_disorder_tool(u_ac,'MobiDB3.0')

        if not mobi_lite:
            processes = {}
            for i in range(1,min(blast_processes,n_disorder) + 1):
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
                proteins.set_disorder_scores(iupred_parts[0],iupred_parts[2])
                proteins.set_disorder_regions(iupred_parts[0],iupred_parts[1])
                proteins.set_disorder_tool(iupred_parts[0],'IUpred')

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
                config.errorlog.add_warning('Warning: mobidb-lite threw an error: %s' % err)
            else:
                entries = []
                for line in out.split('}'):
                    if line.count('{') == 0:
                        continue
                    line = line + '}'
                    mobi_results = json.loads(line)
                    u_ac = mobi_results['acc']
                    raw_scores = mobi_results['p']
                    raw_regions = mobi_results['regions']
                    regions = []
                    for a,b in raw_regions:
                        regions.append([a,b,'disorder'])
                    scores = {}
                    seq = proteins.get_sequence(u_ac)
                    disorder_scores = proteins.get_disorder_scores(u_ac)
                    disorder_regions = proteins.get_disorder_regions(u_ac)
                    for pos,score in enumerate(raw_scores):
                        scores[pos+1] = score
                    proteins.set_disorder_scores(u_ac,scores)
                    proteins.set_disorder_regions(u_ac,regions)
                    proteins.set_disorder_tool(u_ac,'mobidb-lite')


        t1 = time.time()
        if config.verbosity >= 2:
            print("Time for addIupred Part 1: %s" % str(t1-t0))

        if not skip_db:
            background_iu_process = database.addIupred(proteins,config)

            t2 = time.time()
            if config.verbosity >= 2:
                print("Time for addIupred Part 2: %s" % str(t2-t1))
        else:
            background_iu_process = None

    return background_iu_process

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

        p = subprocess.Popen(['python3','%s/iupred2a.py' % iupred_path,seq,'glob'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        out,err = p.communicate()

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
                iupred_parts[2][int(words[0])] = words[2]

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
            structures = blast.blast(seq,target_name,blast_path,blast_db_path,nr=option_number_of_templates,seq_thresh=option_seq_thresh,cov_thresh=option_ral_thresh,cwd = "%s/%d" %(cwd,i))

            structures = templateSelection.selectTemplates(structures,pdb_path)
            structures = templateFiltering.filterTemplates(structures,option_seq_thresh,option_res_thresh,option_ral_thresh)

            with lock:
                out_queue.put((gene,structures))
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc()
            with lock:
                err_queue.put((e,f,g,gene))

#@profile
def autoTemplateSelection(config,manager,lock,proteins,search_tool='MMseqs2'):

    blast_processes = config.blast_processes
    option_seq_thresh = config.option_seq_thresh
    option_res_thresh = config.option_res_thresh
    cwd = os.getcwd()
    pdb_path = config.pdb_path
    errorlog = config.errorlog

    verbose = config.verbose

    if config.verbosity >= 1:
        print('Sequence search with: ',search_tool)

    if search_tool == 'Blast': #Legacy

        t0 = time.time()

        process_list_db = set([])
        process_queue_blast = manager.Queue()
        blast_queues = [process_queue_blast]
        n = 0
        gene_error_map = {}

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
                errorlog.add_error(error_text)

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

        if config.verbosity >= 2:
            print("Template Selection Part 1: %s" % (str(t1-t0)))
            print("Template Selection Part 2: %s" % (str(t2-t12)))
            print("Template Selection Part 3: %s" % (str(t3-t2)))

    elif search_tool == 'MMseqs2':
        t0 = time.time()

        raw_structure_map,pdb_ids = MMseqs2.search(proteins,config)

        t1 = time.time()

        templateSelection.filterRawStructureMap(raw_structure_map,pdb_ids,pdb_path,option_res_thresh,blast_processes,proteins)

        t2 = time.time()
        if config.verbosity >= 2:
            print("Template Selection Part 1: %s" % (str(t1-t0)))
            print("Template Selection Part 2: %s" % (str(t2-t1)))

    return

#@profile
def paraAlignment(config,manager,lock,proteins,skip_db = False):

    alignment_processes = config.alignment_processes
    pdb_path = config.pdb_path
    smiles_path = config.smiles_path
    inchi_path = config.inchi_path
    verbose = config.verbose
    errorlog = config.errorlog
    t0 = time.time()
    input_queue = manager.Queue()
    err_queue = manager.Queue()

    t1 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 1: %s" % (str(t1-t0)))

    if not skip_db:
        db,cursor = config.getDB()
        database.getAlignments(proteins,db,cursor,config)
        db.close()

    in_queues = [input_queue] #partioning the parallized alignment part shall reduce the memory usage, by safing the alignments to the database every 1000 proteins.

    n = 0
    align_sizes = []
    u_acs = proteins.get_protein_u_acs()

    for u_ac in u_acs:
        annotation_list = proteins.get_protein_annotation_list(u_ac)
        for (pdb_id,chain) in annotation_list:

            if proteins.is_annotation_stored(pdb_id,chain,u_ac):
                continue

            oligo = proteins.get_oligo(pdb_id,chain)
            seq = proteins.get_sequence(u_ac)
            aaclist = proteins.getAACList(u_ac)

            input_queue.put((u_ac,pdb_id,chain,oligo,seq,aaclist))

        if n == 1000:
            input_queue = manager.Queue()
            in_queues.append(input_queue)
            align_sizes.append(n)
            n = 0
        n += 1

    align_sizes.append(n)

    t2 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 2: %s" % (str(t2-t1)))

    input_queue = manager.Queue()
    out_queue = manager.Queue()

    n_mappings = 0

    u_acs = proteins.get_protein_u_acs()

    for u_ac in u_acs:
        if not proteins.is_protein_stored(u_ac):
            continue

        all_stored = True
        positions = proteins.get_position_ids(u_ac)

        for pos in positions:
            if not proteins.is_position_stored(u_ac,pos):
                all_stored = False
                break
        if all_stored:
            continue

        annotation_list = proteins.get_protein_annotation_list(u_ac)

        for (pdb_id,chain) in annotation_list:
            if not proteins.is_annotation_stored(pdb_id,chain,u_ac):
                continue
            structure_id = proteins.get_structure_db_id(pdb_id,chain)
            target_seq,template_seq = proteins.get_alignment(u_ac,pdb_id,chain)
            aaclist = proteins.getAACList(u_ac)
            prot_id = proteins.get_protein_db_id(u_ac)
            input_queue.put((u_ac,pdb_id,chain,structure_id,target_seq,template_seq,aaclist,prot_id))
            n_mappings += 1

    processes = {}
    for i in range(1,min(n_mappings,alignment_processes) + 1):
        p = multiprocessing.Process(target=paraMap, args=(config,input_queue,out_queue,lock))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()

    if config.profiling:
        profile = cProfile.Profile()
        profile.enable()

    out_queue.put(None)
    while True:
        out = out_queue.get()
        if out == None:
            break
        (u_ac,pdb_id,chain,sub_infos) = out
        proteins.set_sub_infos(u_ac,pdb_id,chain,sub_infos)

    t3 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 3: %s" % (str(t3-t2)))

    t31 = 0.0
    t4 = 0.0
    t5 = 0.0
    t6 = 0.0
    t7 = 0.0
    t8 = 0.0
    t9 = 0.0

    if config.profiling:
        align_profile = ['u_ac\tpdb_id\tchain\ttime\tseq_len\tlen(alignment_pir)\tpid']

    for pos,input_queue in enumerate(in_queues):

        t31 += time.time()
        alignment_insertion_list = []
        structure_insertion_list = set()
        out_queue = manager.Queue()
        error_queue = manager.Queue()
        err_queue = manager.Queue()
        processes = {}
        align_size = align_sizes[pos] #number of todo alignments in input_queue

        if config.profiling:
            profile.disable()
        for i in range(1,min(align_size,alignment_processes) + 1):
            p = multiprocessing.Process(target=align, args=(config,input_queue,out_queue,error_queue,lock,err_queue))
            processes[i] = p
            p.start()
        for i in processes:
            processes[i].join()
        if config.profiling:
            profile.enable()
            
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
            errorlog.add_error(errortext)

        t4 += time.time()

        error_map = {}
        error_queue.put(None)
        while True:
            err = error_queue.get()
            if err == None:
                break
            (error,u_ac,pdb_id,chain) = err

            proteins.remove_annotation(u_ac,pdb_id,chain)

            if not u_ac in error_map:
                error_map[u_ac] = {}
                #Only 1 error per gene
                error_map[u_ac][(pdb_id,chain)] = error

        for u_ac in error_map:
            for (pdb_id,chain) in error_map[u_ac]:
                errortext = "Parse error during alignment: %s - %s\n\n" % (u_ac,str(error_map[u_ac][(pdb_id,chain)]))
                errorlog.add_error(errortext)

        out_queue.put(None)

        while True:
            out = out_queue.get()
            if out == None:
                break
            (u_ac,pdb_id,chain,alignment_pir,seq_id,coverage,interaction_partners,chain_type_map,oligo,sub_infos,profile_entry) = out

            if config.profiling:
                align_profile.append('\t'.join(profile_entry))

            proteins.set_coverage(u_ac,pdb_id,chain,coverage)
            proteins.set_sequence_id(u_ac,pdb_id,chain,seq_id)
            proteins.set_sub_infos(u_ac,pdb_id,chain,sub_infos)

            proteins.set_interaction_partners(pdb_id,interaction_partners)
            proteins.set_chain_type_map(pdb_id,chain_type_map)

            proteins.set_oligo(pdb_id,chain,oligo)

            structure_insertion_list.add((pdb_id,chain))
            prot_id = proteins.get_protein_db_id(u_ac)
            alignment_insertion_list.append((u_ac,prot_id,pdb_id,chain,alignment_pir))

        del input_queue
        del out_queue
        del error_queue

        t5 += time.time()

        if not skip_db:
            t6 += time.time()
            db,cursor = config.getDB()
            database.insertStructures(structure_insertion_list,db,cursor,smiles_path,inchi_path,pdb_path,proteins)

            t7 += time.time()

            t8 += time.time()
            database.insertAlignments(alignment_insertion_list,proteins,db,cursor,config)
            db.close()
            t9 += time.time()

    if not skip_db:
        db,cursor = config.getDB()
        database.getStoredResidues(proteins,db,cursor,config)
        db.close()

    if config.profiling:
        logfile = '/'.join(config.errorlog_path.split('/')[:-1]) + '/align_profile.tsv'
        f = open(logfile,'w')
        f.write('\n'.join(align_profile))
        f.close()

        profile.disable()
        cum_stats = pstats.Stats(profile)
        cum_stats.sort_stats('cumulative')
        cum_stats.print_stats()

    if config.verbosity >= 2:
        print("Alignment Part 4: %s" % (str(t4-t31)))
        print("Alignment Part 5: %s" % (str(t5-t4)))
        if not skip_db:
            print("Alignment Part 6: %s" % (str(t6-t5)))
            print("Alignment Part 7: %s" % (str(t7-t6)))
            print("Alignment Part 8: %s" % (str(t8-t7)))
            print("Alignment Part 9: %s" % (str(t9-t8)))

    return

def paraMap(config,input_queue,out_queue,lock):
    pdb_path = config.pdb_path
    with lock:
        input_queue.put(None)
    while True:
        with lock:
            inp = input_queue.get()
        if inp == None:
            return

        (u_ac,pdb_id,chain,structure_id,target_seq,template_seq,aaclist,prot_id) = inp

        template_page = pdb.standardParsePDB(pdb_id,pdb_path,obsolete_check=True)
        if template_page == None:
            continue
        seq_res_map = globalAlignment.createTemplateFasta(template_page,pdb_id,chain,onlySeqResMap = True)

        sub_infos,errors,aaclist,update_map = globalAlignment.getSubPos(target_seq,template_seq,aaclist,seq_res_map)

        with lock:
            out_queue.put((u_ac,pdb_id,chain,sub_infos))

    return

#@profile
def align(config,input_queue,out_queue,error_queue,lock,err_queue):
    pdb_path = config.pdb_path
    option_seq_thresh = config.option_seq_thresh

    with lock:
        input_queue.put(None)
    while True:
        with lock:
            inp = input_queue.get()
        if inp == None:

            return
        if config.profiling:
            t0 = time.time()
        (u_ac,pdb_id,chain,oligo,seq,aaclist) = inp

        try:
            (template_page,interaction_partners,chain_type_map,oligo) = pdb.getStandardizedPdbFile(pdb_id,pdb_path,oligo=oligo)

            align_out = globalAlignment.alignBioPython(config,u_ac,seq,pdb_id,template_page,chain,aaclist)

            if isinstance(align_out, str):
                with lock:
                    error_queue.put((align_out,u_ac,pdb_id,chain))
                continue

            (coverage,seq_id,sub_infos,alignment_pir,errors,times,aaclist,update_map) = align_out

            if sub_infos == None or seq_id == None:
                for err in errors:
                    with lock:
                        error_queue.put((err,u_ac,pdb_id,chain))
                continue
            if len(errors) > 0:
                for err in errors:
                    with lock:
                       error_queue.put((err,u_ac,pdb_id,chain))
                continue
            if 100.0*seq_id >= option_seq_thresh:
                with lock:
                    if config.profiling:
                        t1 = time.time()
                        profile_entry = (u_ac,pdb_id,chain,str(t1-t0),str(len(seq)),str(len(alignment_pir)),str(multiprocessing.current_process().pid))
                    else:
                        profile_entry = None
                    out_queue.put((u_ac,pdb_id,chain,alignment_pir,seq_id,coverage,interaction_partners,chain_type_map,oligo,sub_infos,profile_entry))

        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc()
            with lock:
                err_queue.put((e,f,g,u_ac,pdb_id,chain))
    

def classification(proteins,config):

    u_acs = proteins.get_protein_u_acs()

    for u_ac in u_acs:
        package = []

        positions = proteins.get_position_ids(u_ac)

        annotation_list = proteins.get_protein_annotation_list(u_ac)

        iupred_map = proteins.get_disorder_scores(u_ac)

        for (pdb_id,chain) in annotation_list:

            sub_infos = proteins.get_sub_infos(u_ac,pdb_id,chain)

            resolution = proteins.get_resolution(pdb_id)
            chains = proteins.get_complex_chains(pdb_id)

            seq_id = proteins.get_sequence_id(u_ac,pdb_id,chain)
            if seq_id == None:
                continue
            cov = proteins.get_coverage(u_ac,pdb_id,chain)

            qual = templateFiltering.qualityScore(resolution,cov,seq_id)

            for pos in positions:
                aacbase = proteins.get_aac_base(u_ac,pos)
                if not aacbase in sub_infos:
                    continue
                sub_info = sub_infos[aacbase]
                res_nr = sub_info[0]

                if res_nr == None:
                    continue

                if not proteins.contains_residue(pdb_id,chain,res_nr):
                    continue
                res_aa = proteins.get_residue_aa(pdb_id,chain,res_nr)

                identical_aa = res_aa == aacbase[0]

                sld = proteins.get_residue_sld(pdb_id,chain,res_nr)
                scd = proteins.get_residue_scd(pdb_id,chain,res_nr)
                homomer_dists = proteins.get_residue_homomer_dists(pdb_id,chain,res_nr)
                centralities = proteins.get_residue_centralities(pdb_id,chain,res_nr)
                modres = proteins.get_residue_modres(pdb_id,chain,res_nr)
                b_factor = proteins.get_residue_b_factor(pdb_id,chain,res_nr)
                rsa = proteins.get_residue_rsa(pdb_id,chain,res_nr)
                ssa = proteins.get_residue_ssa(pdb_id,chain,res_nr)
                profile = proteins.get_residue_interaction_profile(pdb_id,chain,res_nr)

                phi = proteins.get_residue_phi(pdb_id,chain,res_nr)
                psi = proteins.get_residue_psi(pdb_id,chain,res_nr)

                intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower = proteins.get_residue_link_information(pdb_id,chain,res_nr)

                (inter_chain_median_kd,inter_chain_dist_weighted_kd,
                inter_chain_median_rsa,inter_chain_dist_weighted_rsa,intra_chain_median_kd,
                intra_chain_dist_weighted_kd,intra_chain_median_rsa,intra_chain_dist_weighted_rsa) = proteins.get_residue_milieu(pdb_id,chain,res_nr)
                
                min_hd,min_ld,min_md,min_id,min_cd,min_rd,min_dd,min_lig,min_metal,min_ion,iacs = proteins.structures[(pdb_id,chain)].residues[res_nr].get_shortest_distances(chains)

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

                if min_minimal_distances < 1.2:
                    continue

                if rsa == None:
                    sc = None
                else:
                    if rsa > config.surface_threshold:
                        sc = "Surface"
                    else:
                        sc = "Core"

                raw_rin_class,raw_rin_simple_class = profile.getClass()

                if raw_rin_class == 'No interaction':
                    rin_class = sc
                else:
                    rin_class = raw_rin_class

                if raw_rin_simple_class == 'No interaction':
                    rin_simple_class = sc
                else:
                    rin_simple_class = raw_rin_simple_class

                Class = rin_class
                simpleClass = rin_simple_class

                proteins.add_residue_classification(pdb_id,chain,res_nr,Class,simpleClass)

                mapping = (qual,seq_id,cov,rsa,ssa,min_ld,min_md,min_id,min_cd,min_rd,min_dd,min_hd,profile,centralities,
                            phi,psi,intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
                            inter_chain_median_kd,inter_chain_dist_weighted_kd,
                            inter_chain_median_rsa,inter_chain_dist_weighted_rsa,intra_chain_median_kd,
                            intra_chain_dist_weighted_kd,intra_chain_median_rsa,intra_chain_dist_weighted_rsa,b_factor,modres,
                            Class,simpleClass,identical_aa)

                proteins.add_pos_res_mapping(u_ac,pos,pdb_id,chain,res_nr,mapping)

    for u_ac in proteins.get_protein_u_acs():
        for pos in proteins.get_position_ids(u_ac):
            proteins.classify(u_ac,pos,config)

    return

#@profile
def paraAnnotate(config,manager,lock,proteins, lite = False):
    if config.profiling:
        profile = cProfile.Profile()
        profile.enable()
    annotation_processes = config.annotation_processes
    anno_session_mapping = config.anno_session_mapping
    error_annotations_into_db = config.error_annotations_into_db

    smiles_path = config.smiles_path
    inchi_path = config.inchi_path
    pdb_path = config.pdb_path
    verbose = config.verbose
    errorlog = config.errorlog

    t0 = time.time()
    
    input_queue = manager.Queue()

    #new structure for paralell annotation: give each annotation process a pdb and dict of chain and structure information,
    #all chains without structure information shall go into the database as fully annotated structures, which are not aligned (yet)
    #these are unmapped chains, or protein interaction partners, their presence in the database is necesary for the computation of the structural neighborhood
    #a new structure-structure entry has to be created and has to be inserted into the database before the residue insertion
    #homooligomer information is not present for such cases, this information has to be updated, if the chain is mapped at a later iteration or run

    '''
    sys.path.append("/TL/sin/work/agress/RIN")
    import createRINdb
    createRINdb.calculateRINsFromPdbList(pdb_structure_map.keys(),fromScratch=True,forceCentrality=True)
    '''

    size_map = {}
    max_size = 0

    complex_list = proteins.get_complex_list()

    for pdb_id in complex_list:
        if proteins.is_complex_stored(pdb_id):
            continue
        s = len(proteins.get_complex_chains(pdb_id))
        if not s in size_map:
            size_map[s] = set()
            if s > max_size:
                max_size = s
        size_map[s].add(pdb_id) 

    n_of_structs = 0
    while max_size > 0:
        if max_size in size_map:
            for pdb_id in size_map[max_size]:
                if not lite:
                    input_queue.put(pdb_id)
                else:
                    target_dict = proteins.get_target_dict(pdb_id)
                    input_queue.put((pdb_id,target_dict))

                n_of_structs += 1
        max_size -= 1

    if config.verbosity >= 2:
        t01 = time.time()
        print("Annotation Part 0.1: %s" % (str(t01-t0)))
    

    total_structural_analysis = {}

    interaction_structures = set()

    out_queue = manager.Queue()
    error_queue = manager.Queue()
    err_queue = manager.Queue()
    processes = {}
    if annotation_processes > n_of_structs:
        annotation_processes = n_of_structs
    if config.verbosity >= 1:
        print("Going into Annotation with ",annotation_processes,' processes and ',n_of_structs,'structures')

    if config.profiling:
        profile.disable()
    for i in range(1,annotation_processes + 1):
        p = multiprocessing.Process(target=annotate, args=(config,input_queue,out_queue,error_queue,lock,err_queue,lite))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()
    if config.profiling:
        profile.enable()

    if config.verbosity >= 2:
        t02 = time.time()
        print("Annotation Part 0.2: %s" % (str(t02-t01)))

    err_queue.put(None)
    while True:
        err = err_queue.get()
        if err == None:
            break
        (e,f,g,pdb_id) = err
        errortext = "Annotation Error: %s\n%s\n\n" % (pdb_id,'\n'.join([str(e),str(f),str(g)]))
        errorlog.add_error(errortext)

    if config.verbosity >= 2:
        t1 = time.time()
        print("Annotation Part 1: %s" % (str(t1-t02)))
    
    with lock:
        out_queue.put(None)

    structure_list = proteins.get_structure_list()
    cum_stats = None
    if config.profiling:
        anno_profile = ['u_ac\ttime\tlen(anno_chain_dict)\tpid']

    while True:
        out = out_queue.get()
        if out == None:
            break
        (structural_analysis_dict,pdb_id,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles,proc_id,profile_entry) = out
        if config.profiling:
            anno_profile.append('\t'.join(profile_entry))

            stat_file = 'profile-%s.out' % proc_id
            if os.path.exists(stat_file):
                if cum_stats == None:
                    cum_stats = pstats.Stats(stat_file)
                else:
                    cum_stats.add(stat_file)
                os.remove(stat_file)
        
        for chain in structural_analysis_dict:
            total_structural_analysis[(pdb_id,chain)] = structural_analysis_dict[chain]

            if not (pdb_id,chain) in structure_list:
                interaction_structures.add((pdb_id,chain))

            if lite:
                for res_id in structural_analysis_dict[chain]:
                    residue = structural_analysis_dict[chain][res_id]
                    proteins.add_residue(pdb_id,chain,res_id,residue)

        proteins.set_lig_profile(pdb_id,ligand_profiles)
        proteins.set_ion_profile(pdb_id,ion_profiles)
        proteins.set_metal_profile(pdb_id,metal_profiles)
        proteins.set_chain_chain_profile(pdb_id,chain_chain_profiles)

    del out_queue

    if not lite:
        db,cursor = config.getDB()
        interacting_structure_ids = database.insertInteractingChains(interaction_structures,proteins,db,cursor,smiles_path,inchi_path,pdb_path)

        database.insertComplexes(proteins,db,cursor,smiles_path,inchi_path,pdb_path)

        database.insertResidues(total_structural_analysis,db,cursor,interacting_structure_ids,proteins)

    if config.verbosity >= 2:
        t3 = time.time()
        print("Annotation Part 2: %s" % (str(t3-t1)))

    classification(proteins,config)

    if not lite:

        database.insertClassifications(db,cursor,proteins,config)

        db.close()

    if config.verbosity >= 2:
        t4 = time.time()
        print("Annotation Part 3: %s" % (str(t4-t3)))

    if config.profiling:
        if cum_stats != None:
            cum_stats.add(profile)
            cum_stats.sort_stats('cumulative')
            cum_stats.print_stats()

    if config.profiling:
        logfile = '/'.join(config.errorlog_path.split('/')[:-1]) + '/anno_profile.tsv'
        f = open(logfile,'w')
        f.write('\n'.join(anno_profile))
        f.close()

    return

#@profile
def annotate(config,input_queue,out_queue,error_queue,lock,err_queue,lite):
    if config.profiling:
        pr = cProfile.Profile()
        pr.enable()
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
            if config.profiling:
                pr.disable()
                stat_file = 'profile-%s.out' % multiprocessing.current_process().pid
                pr.dump_stats(stat_file)
            return
        if config.profiling:
            t0 = time.time()

        try:
            if not lite:
                pdb_id = inp

                (annotation_chain_dict,errorlist,ligand_profiles,
                metal_profiles,ion_profiles,chain_chain_profiles) = templateFiltering.structuralAnalysis(pdb_id,pdb_path,dssp_path,rin_db_path,
                    neighborhood_calculation=neighborhood_calculation,dssp=dssp,milieu_threshold=config.milieu_threshold,
                    calculate_interaction_profiles=calculate_interaction_profiles,verbosity = config.verbosity)

            else:
                pdb_id,target_dict = inp
                (annotation_chain_dict,errorlist,ligand_profiles,
                metal_profiles,ion_profiles,chain_chain_profiles) = templateFiltering.structuralAnalysis(pdb_id,pdb_path,dssp_path,rin_db_path,
                    neighborhood_calculation=neighborhood_calculation,dssp=dssp,milieu_threshold=config.milieu_threshold,
                    calculate_interaction_profiles=calculate_interaction_profiles,verbosity = config.verbosity,target_dict = target_dict)

            with lock:
                if config.profiling:
                    t1 = time.time()
                    profile_entry = (pdb_id,str(len(annotation_chain_dict)),str(t1-t0),str(multiprocessing.current_process().pid))
                else:
                    profile_entry = None
                out_queue.put((annotation_chain_dict,pdb_id,ligand_profiles,metal_profiles,
                               ion_profiles,chain_chain_profiles,multiprocessing.current_process().pid,profile_entry))
                if len(errorlist) > 0:
                    for (error,e,f,g) in errorlist:
                        err_queue.put((error,f,g,pdb_id))
        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc()
            with lock:
                err_queue.put((e,f,g,pdb_id))

#@profile
def main(filename,config,output_path,main_file_path):
    n_of_cores = config.proc_n
    mrna_fasta = config.mrna_fasta
    num_of_cores = config.proc_n
    verbose = config.verbose
    search_tool = config.search_tool
    errorlog = config.errorlog
    session = 0 #This can later be used, structman.py could give a specific session id and the pipeline can then expand that session

    background_process_MS = None

    #SyncManager.register('Proteins',sdsc.Proteins)
    #baseManager = SyncManager()
    #baseManager.start()
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
        if config.verbosity >= 1:
            print('Convert vcf file format using Annovar')
        if mrna_fasta != None:
            '... and using mrna file: ',mrna_fasta
        nfname = annovar.annovar_pipeline(filename,config.tax_id,config.annovar_path,config.db_adress,config.db_user_name,config.db_password,anno_db,mrna_fasta,ref_id=config.ref_genome_id)
    else:
        nfname = filename

    mrna_fasta_for_annovar = True #make this a option, for the case of mrna fasta files with non-standard protein-ids (may include further debugging)

    if mrna_fasta_for_annovar:
        mrna_fasta = None

    t01 = time.time()

    if config.verbosity >= 2:
        print("Time for preparation before buildQueue: %s" % (str(t01-t0)))

    junksize = 500

    if config.verbosity >= 1:
        print("Call buildQueue with chunksize: %s and file: %s" % (str(junksize),nfname))
    if config.fasta_input:
        proteins_chunks = parseFasta(config,nfname)
    else:
        proteins_chunks = buildQueue(config,nfname,junksize,mrna_fasta=mrna_fasta)

    t02 = time.time()
    if config.verbosity >= 2:
        print("Time for buildQueue: %s" % (str(t02-t01)))
    if config.verbosity >= 1:
        print("Number of chunks: ",len(proteins_chunks))

    newsession = False
    if session == 0:
        db,cursor = config.getDB()     
        starttime = SQLDateTime()
        session = database.insertSession(starttime,nfname,db,cursor)
        db.close()
        newsession = True

    errorlog.start(nfname,session)

    junk_nr = 1 
    for proteins in proteins_chunks:
        if len(proteins) == 0:
            continue

        #transform the protein map into a Proteins object
        proteins = sdsc.Proteins(proteins)

        if config.verbosity >= 3:
            print('Proteins state after object initialization:')
            proteins.print_protein_state()
        if config.verbosity >= 1:
            print("Chunk %s/%s" % (str(junk_nr),str(len(proteins_chunks))))
        junk_nr+=1
        try:
            os.stat("%s/tmp_structman_pipeline" %(output_path))
        except:
            os.mkdir("%s/tmp_structman_pipeline" %(output_path))  
        os.chdir("%s/tmp_structman_pipeline" %(output_path))
        cwd = "%s/tmp_structman_pipeline" %(output_path)

        try:
            t1 = time.time()
            if config.verbosity >= 1:
                print("Before protCheck")
            #check for already stored genes and make all the necessary database interactions
            #sets the fields database_id and stored in the protein objects as well as the stored and non stored ids in the proteins object
            db,cursor = config.getDB()
            database.protCheck(proteins,session,db,cursor)

            if config.verbosity >= 3:
                print('Proteins state after protCheck:')
                proteins.print_protein_state()

            t2 = time.time()
            if config.verbosity >= 2:
                print("Time for protCheck: %s" % (str(t2-t1)))
            #check for already stored mutations or position twins and make all the necessary database interactions
            #sets the fields database_id and stored in the position objects
            background_process_MS = database.positionCheck(proteins,session,db,cursor,config)

            if config.verbosity >= 3:
                print('Proteins state after positionCheck:')
                proteins.print_protein_state()

            t3 = time.time()
            if config.verbosity >= 2:
                print("Time for positionCheck: %s" % (str(t3-t2)))

            if config.verbosity >= 1:
                print("Before getSequences")
            background_iu_process = getSequences(proteins,config,manager,lock)
            db.close()

            t4 = time.time()
            if config.verbosity >= 2:
                print("Time for getSequences: %s" % (str(t4-t3)))

            if config.verbosity >= 1:
                print("Before autoTemplateSelection")
            autoTemplateSelection(config,manager,lock,proteins,search_tool=search_tool)

            t5 = time.time()
            if config.verbosity >= 2:
                print("Time for Template Selection: %s" % (str(t5-t4)))

            if config.verbosity >= 1:
                print("Before paraAlignment")
            paraAlignment(config,manager,lock,proteins)

            t6 = time.time()
            if config.verbosity >= 2:
                print("Time for Alignment: %s" % (str(t6-t5)))

            if background_iu_process != None:
                background_iu_process.join() #Disorder values have to finished before the classification happens

            if config.verbosity >= 1:
                print("Before paraAnnotate")
            paraAnnotate(config,manager,lock,proteins)

            t7 = time.time()
            if config.verbosity >= 2:
                print("Time for Annotation: %s" % (str(t7-t6)))

            t1 = time.time()
            #join the background inserts
            if background_process_MS != None:
                background_process_MS.join()
            #if background_process_AS != None:
            #    background_process_AS.join()


            t2 = time.time()
            if config.verbosity >= 2:
                print('Resttime for background inserts: ',t2-t1)

        #Error-Handling for a whole input line
        except:
            
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc()
            #print "Pipeline Core Error: ",e,f,g
            errortext = '\n'.join([str(e),str(f),str(g)]) + '\n\n'
            errorlog.add_error(errortext)
            if background_process_MS != None:
                background_process_MS.join()

        os.chdir(output_path)
        #"""
        try:
            shutil.rmtree(cwd)
        except:
            pass
        #"""

    errorlog.stop()

    db,cursor = config.getDB()
    if newsession:
        endtime = SQLDateTime()
        database.updateSession(session,endtime,db,cursor)
    db.close()

    tend = time.time()
    if config.verbosity >= 2:
        print((tend-t0))
    return session
