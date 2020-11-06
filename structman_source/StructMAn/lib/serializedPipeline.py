import pymysql as MySQLdb
import os
import sys
import errno
import getopt
import time
import multiprocessing
import subprocess
import shutil
import psutil
import ray

import traceback
import math
import json

import pdbParser
import templateSelection
import templateFiltering
import globalAlignment
import database
import uniprot
import blast
import annovar
import MMseqs2
import sdsc
import indel_analysis
import rin
import output

import cProfile
import pstats
import io
import gc

try:
    from memory_profiler import profile
except:
    pass

#Taken from https://stackoverflow.com/questions/2023608/check-what-files-are-open-in-python
def list_fds():
    """List process currently open FDs and their target """
    if not sys.platform.startswith('linux'):
        raise NotImplementedError('Unsupported platform: %s' % sys.platform)

    ret = {}
    base = '/proc/self/fd'
    for num in os.listdir(base):
        path = None
        try:
            path = os.readlink(os.path.join(base, num))
        except OSError as err:
            # Last FD is always the "listdir" one (which may be closed)
            if err.errno != errno.ENOENT:
                raise
        ret[int(num)] = path

    return ret

def SQLDateTime():
    return time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())

#For memory profiling, thanks to Fred Cirera
def sizeof_fmt(num, suffix='B'):
    ''' by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified'''
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)

#@profile
def sequenceScan(config,proteins,indels):
    pdb_path = config.pdb_path

    sequenceScanProteins = {}
    sequenceScanPDB = {}

    for u_ac in proteins:
        uni_pos,tags = proteins[u_ac].popNone()
        
        if u_ac.count(':') > 0:
            sequenceScanPDB[u_ac] = tags,uni_pos #PDB inputs are always processed by the sequence scan, but the positions are only added if uni_pos is true
        elif not sdsc.is_mutant_ac(u_ac):
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

            tags,uni_pos = sequenceScanProteins[u_ac]

            for (pos,aa) in enumerate(seq):
                seq_pos = pos+1
                if seq_pos in proteins[u_ac].positions:
                    proteins[u_ac].positions[seq_pos].check(aa,overwrite = config.overwrite_incorrect_wt_aa)
                    proteins[u_ac].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    position = sdsc.Position(pos=seq_pos,wt_aa=aa,tags=tags,checked=True)
                    proteins[u_ac].positions[seq_pos] = position

            proteins[u_ac].set_disorder_scores(disorder_scores)
            proteins[u_ac].set_disorder_regions(disorder_regions_datastruct)

            #sanity check filter
            del_list = []
            for pos in proteins[u_ac].positions:
                if not proteins[u_ac].positions[pos].checked:
                    del_list.append(pos)
            for pos in del_list:
                warn_text = 'Filtered position through sanity check:',u_ac,pos
                if config.verbosity >= 2:
                    print(warn_text)
                config.errorlog.add_warning(warn_text)
                del proteins[u_ac].positions[pos]

    if len(sequenceScanPDB) > 0:
        pdb_sequence_map,pdb_pos_map = pdbParser.getSequences(sequenceScanPDB.keys(),pdb_path)
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
                    position = sdsc.Position(pos=seq_pos,pdb_res_nr=res_id,wt_aa=aa,tags=tags)
                    proteins[u_ac].positions[seq_pos] = position
                    proteins[u_ac].res_id_map[res_id] = position

    for indel in indels:
        indel.mutate_sequence(proteins)

    return proteins,indels

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
def buildQueue(config,filename,chunksize):
    t0 =time.time()

    proteins = {}
    u_ids = set()
    u_acs = set()
    id_map = {}
    ac_map = {}
    np_map = {}

    pdb_map = {}

    if isinstance(filename,str):
        f = open(filename, "r")
        lines = f.read().split('\n')
        f.close()
    else:
        #In case of single line input
        lines = ['\t'.join(filename)]

    tag_map = {}

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
                indel = None

            else:
                aachange = words[1].replace("\n","")
                
                if not sp_id.count(':') == 1: #this means sp_id is not a pdb-id
                    if aachange.count('delins') == 1:
                        flanks,inserted_sequence = aachange.split('delins')
                        left_flank,right_flank = flanks.split('_')

                        indel = sdsc.Substitution(left_flank = left_flank,right_flank = right_flank,inserted_sequence = inserted_sequence,tags = tags)
                    elif aachange.count('del') == 1:
                        flanks = aachange.split('del')[0]
                        left_flank,right_flank = flanks.split('_')

                        indel = sdsc.Deletion(left_flank = left_flank,right_flank = right_flank,tags = tags)

                    elif aachange.count('ins') == 1:
                        flanks,inserted_sequence = aachange.split('ins')
                        if flanks.count('_') == 1:
                            left_flank,right_flank = flanks.split('_')
                        elif flanks[1:] == '1':
                            left_flank = flanks
                            right_flank = None
                        else:
                            left_flank = None
                            right_flank = flanks

                        indel = sdsc.Insertion(left_flank = left_flank,right_flank = right_flank,inserted_sequence = inserted_sequence,tags = tags)

                    else:
                        indel = None
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
                    indel = None
                    if aachange.count('_') == 1: #This means a mutant amino acid type is given
                        aachange = aachange.replace('_','')
                        position = sdsc.Position(pdb_res_nr=aachange[1:-1],wt_aa=aachange[0],mut_aas = set(aachange[-1]),tags=tags)
                    else:
                        position = sdsc.Position(pdb_res_nr=aachange[1:],wt_aa=aachange[0],tags=tags)

        except:
            config.errorlog.add_error("File Format Error: %s" % line)

        if sp_id[2] == "_":
            if not sp_id in np_map:
                np_map[sp_id] = [],[]
            if indel == None:
                np_map[sp_id][0].append(position)
            else:
                np_map[sp_id][1].append(indel)
            
        else:
            if sp_id.count('_') > 0:
                u_ids.add(sp_id)
                if not sp_id in id_map:
                    id_map[sp_id] = [],[]
                if indel == None:
                    id_map[sp_id][0].append(position)
                else:
                    id_map[sp_id][1].append(indel)
            elif len(sp_id) == 6 and sp_id.count(':') == 1:
                pdb_chain_tuple = '%s:%s' % (sp_id[:4].upper(),sp_id[-1]) #enforce uppercase pdb-id 
                if not pdb_chain_tuple in pdb_map:
                    pdb_map[pdb_chain_tuple] = [position]
                else:
                    pdb_map[pdb_chain_tuple].append(position)
            else:
                u_acs.add(sp_id)
                if not sp_id in ac_map:
                    ac_map[sp_id] = [],[]
                if indel == None:
                    ac_map[sp_id][0].append(position)
                else:
                    ac_map[sp_id][1].append(indel)


    t1 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 1: ",str(t1-t0))

    proteins,indels = uniprot.IdMapping(config,ac_map,id_map,np_map,pdb_map)
    
    t2 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 2: ",str(t2-t1))

    proteins,indels = sequenceScan(config,proteins,indels)

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

    amount_of_indel_proteins = 2*len(indels)

    if s > chunksize:

        n_of_indel_batches = amount_of_indel_proteins//chunksize
        if amount_of_indel_proteins%chunksize != 0:
            n_of_indel_batches += 1

        n_of_batches = s//chunksize
        batchsize = s//n_of_batches
        if n_of_indel_batches > 0:
            indel_rest = s%n_of_indel_batches
        outlist = []
        only_indels = s == amount_of_indel_proteins

        new_map = {}
        for i in range(0,n_of_indel_batches):
            new_map = {}
            indel_batch = []
            if i == n_of_indel_batches - 1: # the last batch
                batchsize += indel_rest
            for j in range(0,batchsize/2):
                indel = indels[j+batchsize*i]
                new_map[indel.wt_prot] = proteins[indel.wt_prot]
                del proteins[indel.wt_prot]
                new_map[indel.mut_prot] = proteins[indel.mut_prot]
                del proteins[indel.mut_prot]
                indel_batch.append(indel)
            if only_indels:
                outlist.append((new_map,indel_batch))
            elif not i == n_of_indel_batches - 1:
                outlist.append((new_map,indel_batch))

        mix_map = len(new_map) > 0

        n_of_batches = n_of_batches - len(outlist)
        if s%chunksize != 0:
            n_of_batches += 1
        rest = s%n_of_batches

        for i in range(0,n_of_batches):
            if not mix_map:
                new_map = {}
                mix_map = False
                effective_batchsize = batchsize
                indel_map = []
            else:
                effective_batchsize = batchsize - len(new_map)
                mix_map = False

            if i == n_of_batches - 1: # the last batch
                effective_batchsize += rest
            for j in range(0,effective_batchsize):
                if len(proteins) == 0:
                    continue
                (key,value) = proteins.popitem()
                new_map[key] = value
            outlist.append((new_map,indel_map))
        new_map = {}
        while len(proteins) > 0:
            (key,value) = proteins.popitem()
            new_map[key] = value
        if len(new_map) > 0:
            outlist.append((new_map,indel_mapap))
    else:
        outlist.append((proteins,indels))

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

    if not skip_db:
        database.addProtInfos(proteins,config)

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
            if disorder_scores == None:
                if not mobi_lite:
                    in_queue.put((u_ac,seq))
                    n_disorder += 1
                else:
                    mobi_list.append('>%s\n%s\n' % (u_ac,seq))
            elif disorder_scores != 'Stored':
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

            del in_queue
            del out_queue

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
                    #seq = proteins.get_sequence(u_ac)
                    #disorder_scores = proteins.get_disorder_scores(u_ac)
                    #disorder_regions = proteins.get_disorder_regions(u_ac)
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

#@profile
def paraIupred(config,in_queue,out_queue,lock):
    iupred_path = config.iupred_path
    with lock:
        in_queue.put(None)
    while True:
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
def autoTemplateSelection(config,manager,lock,proteins):

    blast_processes = config.blast_processes
    option_seq_thresh = config.option_seq_thresh
    option_res_thresh = config.option_res_thresh
    cwd = os.getcwd()
    pdb_path = config.pdb_path
    errorlog = config.errorlog

    verbose = config.verbose
    search_tool= config.search_tool

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

        templateSelection.filterRawStructureMap(raw_structure_map,pdb_ids,pdb_path,option_res_thresh,blast_processes,proteins,manager,lock)

        t2 = time.time()
        if config.verbosity >= 2:
            print("Template Selection Part 1: %s" % (str(t1-t0)))
            print("Template Selection Part 2: %s" % (str(t2-t1)))

    return

#@profile
def paraAlignment(config,manager,lock,proteins,skip_db = False):

    alignment_processes = config.alignment_processes
    errorlog = config.errorlog
    t0 = time.time()

    t1 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 1: %s" % (str(t1-t0)))

    if not skip_db:
        database.getAlignments(proteins,config)

    alignment_results = []

    u_acs = proteins.get_protein_u_acs()

    mapping_dump = ray.put(config)

    for u_ac in u_acs:
        annotation_list = proteins.get_protein_annotation_list(u_ac)
        for (pdb_id,chain) in annotation_list:

            if proteins.is_annotation_stored(pdb_id,chain,u_ac):
                continue

            oligo = proteins.get_oligo(pdb_id,chain)
            seq = proteins.get_sequence(u_ac)
            aaclist = proteins.getAACList(u_ac)

            alignment_results.append(align.remote(mapping_dump,u_ac,pdb_id,chain,oligo,seq,aaclist))

    t2 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 2: %s" % (str(t2-t1)))

    mapping_results = []
    u_acs = proteins.get_protein_u_acs()

    for u_ac in u_acs:
        if not proteins.is_protein_stored(u_ac):
            continue

        if proteins.is_completely_stored(u_ac):
            continue

        annotation_list = proteins.get_protein_annotation_list(u_ac)

        for (pdb_id,chain) in annotation_list:
            if not proteins.is_annotation_stored(pdb_id,chain,u_ac):
                continue
            structure_id = proteins.get_structure_db_id(pdb_id,chain)
            target_seq,template_seq = proteins.get_alignment(u_ac,pdb_id,chain)
            aaclist = proteins.getAACList(u_ac)
            prot_id = proteins.get_protein_db_id(u_ac)
            mapping_results.append(paraMap.remote(mapping_dump,u_ac,pdb_id,chain,structure_id,target_seq,template_seq,aaclist,prot_id))

    mappings_outs = ray.get(mapping_results)
    for (u_ac,pdb_id,chain,sub_infos) in mappings_outs:
        proteins.set_sub_infos(u_ac,pdb_id,chain,sub_infos)

    t3 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 3: %s" % (str(t3-t2)))

    align_outs = ray.get(alignment_results)
    alignment_insertion_list = []
    structure_insertion_list = set()
    warn_map = set()

    for out in align_outs:

        if len(out) == 3:
            u_ac,pdb_id,chain = out
            proteins.remove_annotation(u_ac,pdb_id,chain)
            continue

        if len(out) == 4:
            (u_ac,pdb_id,chain,warn_text) = out
            proteins.remove_annotation(u_ac,pdb_id,chain)
            if not u_ac in warn_map:
                errorlog.add_warning(warn_text)
            warn_map.add(u_ac)
            continue

        (u_ac,pdb_id,chain,alignment,seq_id,coverage,interaction_partners,chain_type_map,oligo,sub_infos) = out

        proteins.set_coverage(u_ac,pdb_id,chain,coverage)
        proteins.set_sequence_id(u_ac,pdb_id,chain,seq_id)
        proteins.set_sub_infos(u_ac,pdb_id,chain,sub_infos)
        proteins.set_alignment(u_ac,pdb_id,chain,alignment)

        proteins.set_interaction_partners(pdb_id,interaction_partners)
        proteins.set_chain_type_map(pdb_id,chain_type_map)

        proteins.set_oligo(pdb_id,chain,oligo)

        structure_insertion_list.add((pdb_id,chain))
        prot_id = proteins.get_protein_db_id(u_ac)
        alignment_insertion_list.append((u_ac,prot_id,pdb_id,chain,alignment))

    if config.verbosity >= 2:
        t4 = time.time()
        print("Alignment Part 4: %s" % (str(t4-t3)))
        

    if not skip_db:
        if config.verbosity >= 2:
            t5 = time.time()
            print("Alignment Part 5: %s" % (str(t5-t4)))

        database.insertStructures(structure_insertion_list,proteins,config)

        if config.verbosity >= 2:
            t6 = time.time()
            print("Alignment Part 6: %s" % (str(t6-t5)))

        database.insertAlignments(alignment_insertion_list,proteins,config)

        if config.verbosity >= 2:
            t7 = time.time()
            print("Alignment Part 7: %s" % (str(t7-t6)))

        database.getStoredResidues(proteins,config)

        if config.verbosity >= 2:
            t8 = time.time()
            print("Alignment Part 8: %s" % (str(t8-t7)))

    return

@ray.remote(num_cpus=1)
def paraMap(mapping_dump,u_ac,pdb_id,chain,structure_id,target_seq,template_seq,aaclist,prot_id):
    #hack proposed by the devs of ray to prevent too many processes being spawned
    resources = ray.ray.get_resource_ids() 
    cpus = [v[0] for v in resources['CPU']]
    psutil.Process().cpu_affinity(cpus)

    config = mapping_dump

    pdb_path = config.pdb_path

    template_page = pdbParser.standardParsePDB(pdb_id,pdb_path,obsolete_check=True)

    seq_res_map = globalAlignment.createTemplateFasta(template_page,pdb_id,chain,config,onlySeqResMap = True)

    sub_infos,aaclist = globalAlignment.getSubPos(config,u_ac,target_seq,template_seq,aaclist,seq_res_map)

    return (u_ac,pdb_id,chain,sub_infos)

@ray.remote(num_cpus=1)
def align(align_dump,u_ac,pdb_id,chain,oligo,seq,aaclist):
    #hack proposed by the devs of ray to prevent too many processes being spawned
    resources = ray.ray.get_resource_ids() 
    cpus = [v[0] for v in resources['CPU']]
    psutil.Process().cpu_affinity(cpus) 

    config = align_dump

    pdb_path = config.pdb_path
    option_seq_thresh = config.option_seq_thresh

    try:
        parse_out = pdbParser.getStandardizedPdbFile(pdb_id,pdb_path,oligo=oligo)

        if parse_out == None:
            return (u_ac,pdb_id,chain,'pdbParser failed')

        (template_page,interaction_partners,chain_type_map,oligo) = parse_out

        align_out = globalAlignment.alignBioPython(config,u_ac,seq,pdb_id,template_page,chain,aaclist)

        if isinstance(align_out, str):
            return (u_ac,pdb_id,chain,align_out)

        (coverage,seq_id,sub_infos,alignment_pir,times,aaclist) = align_out

        if sub_infos == None or seq_id == None:
            return (u_ac,pdb_id,chain)

        if 100.0*seq_id >= option_seq_thresh:
            return (u_ac,pdb_id,chain,alignment_pir,seq_id,coverage,interaction_partners,chain_type_map,oligo,sub_infos)
        else:
            return (u_ac,pdb_id,chain)

    except:
        [e,f,g] = sys.exc_info()
        g = traceback.format_exc()
        return (u_ac,pdb_id,chain,'%s,%s\n%s\n%s\n%s' % (pdb_id,chain,e,str(f),g))

@ray.remote(num_cpus=1)
def para_classify_remote_wrapper(classification_dump,package):
    #hack proposed by the devs of ray to prevent too many processes being spawned
    resources = ray.ray.get_resource_ids() 
    cpus = [v[0] for v in resources['CPU']]
    psutil.Process().cpu_affinity(cpus)
    return para_classify(classification_dump,package)

def para_classify(classification_dump,package):
    t0 = time.time()
    config = classification_dump

    outs = []

    for u_ac,classification_inp in package:
        for pos,mappings,disorder_score,disorder_region in classification_inp:

            mappings_obj = sdsc.Mappings()

            for mapping in mappings:

                (pdb_id,chain,res_nr,qual,seq_id,cov,rsa,ssa,profile_or_str,centralities_or_str,
                phi,psi,intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
                inter_chain_median_kd,inter_chain_dist_weighted_kd,
                inter_chain_median_rsa,inter_chain_dist_weighted_rsa,intra_chain_median_kd,
                intra_chain_dist_weighted_kd,intra_chain_median_rsa,intra_chain_dist_weighted_rsa,b_factor,modres,
                identical_aa,lig_dists,chain_distances,homomer_distances,chains,resolution,res_aa) = mapping

                gsd_return = sdsc.get_shortest_distances(chains,lig_dists,chain_distances,homomer_distances)

                if gsd_return == None:
                    continue
                homo_dist,lig_dist,metal_dist,ion_dist,chain_dist,rna_dist,dna_dist,min_lig,min_metal,min_ion,iacs = gsd_return

                if isinstance(profile_or_str,str):
                    profile = rin.Interaction_profile(profile_str = profile_or_str)
                else:
                    profile = profile_or_str

                if isinstance(centralities_or_str,str):
                    centralities = rin.Centrality_scores(code_str = centralities_or_str)
                else:
                    centralities = centralities_or_str

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


                mapping = (qual,seq_id,cov,rsa,ssa,lig_dist,metal_dist,ion_dist,chain_dist,rna_dist,dna_dist,homo_dist,profile,centralities,
                            phi,psi,intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
                            inter_chain_median_kd,inter_chain_dist_weighted_kd,
                            inter_chain_median_rsa,inter_chain_dist_weighted_rsa,intra_chain_median_kd,
                            intra_chain_dist_weighted_kd,intra_chain_median_rsa,intra_chain_dist_weighted_rsa,b_factor,modres,rin_class,rin_simple_class,
                            identical_aa,resolution,res_aa)

                mappings_obj.add_mapping((pdb_id,chain,res_nr),mapping)

            mappings_obj.weight_all(config,disorder_score,disorder_region)

            mapping_results = mappings_obj.get_raw_result()

            outs.append((u_ac,pos,mapping_results))

    t1 = time.time()

    return outs,t1-t0

def classification(proteins,config):

    if config.verbosity >= 2:
        t0 = time.time()

    classification_results = []

    if config.verbosity >= 2:
        t11 = time.time()
        print('Time for classification part 1.1:',t11-t0)

    size_sorted = []
    u_acs = proteins.get_protein_u_acs()
    for u_ac in u_acs:
        if proteins.is_completely_stored(u_ac):
            continue
        annotation_list = proteins.get_protein_annotation_list(u_ac)
        size_sorted.append((u_ac,len(annotation_list)))

    size_sorted = sorted(size_sorted,reverse = True, key= lambda x:x[1])

    packages = []
    package_counter = 0
    counter_direction_forward = True

    for u_ac,size in size_sorted:
        classification_inp = []
        positions = proteins.get_position_ids(u_ac)
        for pos in positions:
            if proteins.is_position_stored(u_ac,pos):
                continue
            mappings = []
            aacbase = proteins.get_aac_base(u_ac,pos)
            disorder_score = proteins.get_disorder_score(u_ac,pos)
            disorder_region = proteins.get_disorder_region(u_ac,pos)
            annotation_list = proteins.get_protein_annotation_list(u_ac)
            for (pdb_id,chain) in annotation_list:

                sub_infos = proteins.get_sub_infos(u_ac,pdb_id,chain)

                if not pos in sub_infos:
                    if config.verbosity >= 3:
                        print('Skipped classification of',u_ac,'due to pos',pos,'was not in sub_infos (len:',len(sub_infos),')')
                    continue

                resolution = proteins.get_resolution(pdb_id)
                chains = proteins.get_complex_chains(pdb_id)

                seq_id = proteins.get_sequence_id(u_ac,pdb_id,chain)
                if seq_id == None:
                    continue
                cov = proteins.get_coverage(u_ac,pdb_id,chain)

                qual = templateFiltering.qualityScore(resolution,cov,seq_id)

                sub_info = sub_infos[pos]
                res_nr = sub_info[0]

                if res_nr == None:
                    continue

                if not proteins.contains_residue(pdb_id,chain,res_nr):
                    if config.verbosity >= 3:
                        print('Skipped classification of',u_ac,pos,'due to res_nr',res_nr,'was not contained in ',pdb_id,chain)
                    continue
                res_aa = proteins.get_residue_aa(pdb_id,chain,res_nr)

                identical_aa = res_aa == aacbase[0]

                centralities_or_str = proteins.get_residue_centralities(pdb_id,chain,res_nr,get_whats_there = True)
                modres = proteins.get_residue_modres(pdb_id,chain,res_nr)
                b_factor = proteins.get_residue_b_factor(pdb_id,chain,res_nr)
                rsa = proteins.get_residue_rsa(pdb_id,chain,res_nr)
                ssa = proteins.get_residue_ssa(pdb_id,chain,res_nr)
                profile_or_str = proteins.get_residue_interaction_profile(pdb_id,chain,res_nr,get_whats_there = True)

                phi = proteins.get_residue_phi(pdb_id,chain,res_nr)
                psi = proteins.get_residue_psi(pdb_id,chain,res_nr)

                intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower = proteins.get_residue_link_information(pdb_id,chain,res_nr)

                (inter_chain_median_kd,inter_chain_dist_weighted_kd,
                inter_chain_median_rsa,inter_chain_dist_weighted_rsa,intra_chain_median_kd,
                intra_chain_dist_weighted_kd,intra_chain_median_rsa,intra_chain_dist_weighted_rsa) = proteins.get_residue_milieu(pdb_id,chain,res_nr)
                
                lig_dists = proteins.get_residue_sld(pdb_id,chain,res_nr)
                chain_distances = proteins.get_residue_scd(pdb_id,chain,res_nr)
                homomer_distances = proteins.get_residue_homomer_dists(pdb_id,chain,res_nr)

                mappings.append((pdb_id,chain,res_nr,qual,seq_id,cov,rsa,ssa,profile_or_str,centralities_or_str,
                            phi,psi,intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
                            inter_chain_median_kd,inter_chain_dist_weighted_kd,
                            inter_chain_median_rsa,inter_chain_dist_weighted_rsa,intra_chain_median_kd,
                            intra_chain_dist_weighted_kd,intra_chain_median_rsa,intra_chain_dist_weighted_rsa,b_factor,modres,
                            identical_aa,lig_dists,chain_distances,homomer_distances,chains,resolution,res_aa))

            classification_inp.append((pos,mappings,disorder_score,disorder_region))


        if len(packages) == package_counter:
            packages.append([])
        packages[package_counter].append((u_ac,classification_inp))

        if package_counter == config.annotation_processes - 1 and counter_direction_forward:
            counter_direction_forward = False
        elif package_counter == config.annotation_processes - 1:
            package_counter -= 1
            if package_counter == -1:
                package_counter = 0
        elif package_counter == 0 and not counter_direction_forward:
            counter_direction_forward = True
        elif package_counter == 0:
            package_counter += 1
        elif counter_direction_forward:
            package_counter += 1
        else:
            package_counter -= 1

    para = len(packages) > 1 #only use ray, when we have more than one package

    if len(packages) == 0:
        return

    if config.verbosity >= 2:
        t11 = time.time()
        print('Time for classification part 1.1:',t11-t0)

    if para:
        classification_dump = ray.put(config)
        for package in packages:
            classification_results.append(para_classify_remote_wrapper.remote(classification_dump,package))
    else:
        classification_results.append(para_classify(config,packages[0]))

    if config.verbosity >= 2:
        t12 = time.time()
        print('Time for classification part 1.2:',t12-t11)

    for u_ac in u_acs:
        proteins.remove_protein_annotations(u_ac)
    proteins.semi_deconstruct()
    gc.collect()

    if config.verbosity >= 2:
        t13 = time.time()
        print('Time for classification part 1.3:',t13-t12)

    if config.verbosity >= 2:
        t1 = time.time()
        print('Time for classification part 1:',t1-t0,para)

    if para:
        max_comp_time = 0
        total_comp_time = 0.
        max_comp_time_pos = None
        amount_of_positions = 0

        while True:
            ready, not_ready = ray.wait(classification_results)

            if len(ready) > 0:
                for outs,comp_time in ray.get(ready):
                    for u_ac,pos,mapping_results in outs:
                        proteins.protein_map[u_ac].positions[pos].mappings = sdsc.Mappings(raw_results = mapping_results)

                    total_comp_time += comp_time
                    amount_of_positions += 1

                    if comp_time > max_comp_time:
                        max_comp_time = comp_time
                        max_comp_time_pos = u_ac
                classification_results = not_ready

            if len(not_ready) == 0:
                break
    else:
        for outs,comp_time in classification_results:
            for u_ac,pos,mapping_results in outs:
                proteins.protein_map[u_ac].positions[pos].mappings = sdsc.Mappings(raw_results = mapping_results)

    if config.verbosity >= 2:
        t2 = time.time()
        print('Time for classification part 2:',t2-t1)
        if para:
            print('Longest computation for:',max_comp_time_pos,'with:',max_comp_time,'In total',amount_of_positions,'proteins','Accumulated time:',total_comp_time)

    return

#@profile
def paraAnnotate(config,manager,lock,proteins, lite = False):

    annotation_processes = config.annotation_processes

    errorlog = config.errorlog

    t0 = time.time()
    
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
    n_of_chain_thresh = 10 #Structures with n_of_chain_thresh or more chains get a nested paralellization

    size_map = {}

    total_structural_analysis = {}
    interaction_structures = set()
    structure_list = proteins.get_structure_list()

    complex_list = proteins.get_complex_list()

    for pdb_id in complex_list:
        if proteins.is_complex_stored(pdb_id):
            continue
        s = len(proteins.get_complex_chains(pdb_id))
        if not s in size_map:
            size_map[s] = set()
        size_map[s].add(pdb_id)

    sorted_sizes = sorted(size_map.keys(),reverse = True)

    anno_result_ids = []

    package_cost = 0
    assigned_costs = {}
    conf_dump = ray.put(config)
    for s in sorted_sizes:
        if s < n_of_chain_thresh:
            cost = 1
        else:
            cost = min([s,config.annotation_processes])
        
        del_list = []
        for pdb_id in size_map[s]:

            if cost + package_cost <= config.annotation_processes:
                if not lite:
                    target_dict = None
                else:
                    target_dict = proteins.get_target_dict(pdb_id)

                anno_result_ids.append(annotate.remote(conf_dump,pdb_id,target_dict))
                package_cost += cost
                del_list.append(pdb_id)
                assigned_costs[pdb_id] = cost
            else:
                break

        for pdb_id in del_list:
            size_map[s].remove(pdb_id)
        if len(size_map[s]) == 0:
            del size_map[s]
    sorted_sizes = sorted(size_map.keys(),reverse = True)


    if config.verbosity >= 2:
        t11 = time.time()
        print("Annotation Part 1.1: %s" % (str(t11-t0)))

    #anno_results = ray.get(anno_result_ids)

    if config.verbosity >= 2:
        t1 = time.time()
        print("Annotation Part 1.2: %s" % (str(t1-t11)))
        print("Annotation Part 1: %s" % (str(t1-t0)))

    max_comp_time = 0.
    total_comp_time = 0.
    max_comp_time_structure = None
    amount_of_structures = 0

    while True:
        ready, not_ready = ray.wait(anno_result_ids)

        if len(ready) > 0:
            for (ret_pdb_id,structural_analysis_dict,errorlist,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles,comp_time) in ray.get(ready):
                new_anno_result_ids = []
                if len(size_map) > 0:
                #Start new jobs regarding the freed resources
                    freed_cost = assigned_costs[ret_pdb_id]
                    package_cost -= freed_cost
                    
                    for s in sorted_sizes:
                        if s < n_of_chain_thresh:
                            cost = 1
                        else:
                            cost = min([s,config.annotation_processes])

                        del_list = []
                        for pdb_id in size_map[s]:
                            if cost + package_cost <= config.annotation_processes:
                                if not lite:
                                    target_dict = None
                                else:
                                    target_dict = proteins.get_target_dict(pdb_id)

                                new_anno_result_ids.append(annotate.remote(conf_dump,pdb_id,target_dict))
                                package_cost += cost
                                del_list.append(pdb_id)
                                assigned_costs[pdb_id] = cost
                            else:
                                break

                        for pdb_id in del_list:
                            size_map[s].remove(pdb_id)
                        if len(size_map[s]) == 0:
                            del size_map[s]
                    sorted_sizes = sorted(size_map.keys(),reverse = True)

                if comp_time > max_comp_time:
                    max_comp_time = comp_time
                    max_comp_time_structure = ret_pdb_id
                amount_of_structures += 1
                total_comp_time += comp_time
                if len(errorlist) > 0:
                    for error_text in errorlist:
                        errorlog.add_warning(error_text)

                for chain in structural_analysis_dict:
                    if lite:
                        for res_id in structural_analysis_dict[chain]:
                            residue = structural_analysis_dict[chain][res_id]
                            proteins.add_residue(ret_pdb_id,chain,res_id,residue)
                    else:
                        total_structural_analysis[(ret_pdb_id,chain)] = structural_analysis_dict[chain]

                        if not (ret_pdb_id,chain) in structure_list:
                            interaction_structures.add((ret_pdb_id,chain))

                proteins.set_lig_profile(ret_pdb_id,ligand_profiles)
                proteins.set_ion_profile(ret_pdb_id,ion_profiles)
                proteins.set_metal_profile(ret_pdb_id,metal_profiles)
                proteins.set_chain_chain_profile(ret_pdb_id,chain_chain_profiles)
            anno_result_ids = new_anno_result_ids + not_ready

        if len(anno_result_ids) == 0 and len(size_map) == 0:
            break

    if config.verbosity >= 2:
        t2 = time.time()
        print("Annotation Part 2: %s" % (str(t2-t1)))
        print('Longest computation for:',max_comp_time_structure,'with:',max_comp_time,'In total',amount_of_structures,'structures','Accumulated time:',total_comp_time)

    if not lite:
        interacting_structure_ids = database.insertInteractingChains(interaction_structures,proteins,config)

        if config.verbosity >= 2:
            t32 = time.time()
            print('Annotation Part 3.1:',t32-t2)

        database.insertComplexes(proteins,config)

        if config.verbosity >= 2:
            t33 = time.time()
            print('Annotation Part 3.2:',t33-t32)

        background_insert_residues_process = database.insertResidues(total_structural_analysis,interacting_structure_ids,proteins,config)

        if config.verbosity >= 2:
            t34 = time.time()
            print('Annotation Part 3.3:',t34-t33)
    else:
        background_insert_residues_process = None

    if config.verbosity >= 2:
        t3 = time.time()
        print("Annotation Part 3: %s" % (str(t3-t2)))

    classification(proteins,config)

    if config.verbosity >= 2:
        t4 = time.time()
        print("Annotation Part 4: %s" % (str(t4-t3)))

    if not lite:
        database.insertClassifications(proteins,config)

    if config.verbosity >= 2:
        t5 = time.time()
        print("Annotation Part 5: %s" % (str(t5-t4)))

    gc.collect()

    if config.verbosity >= 2:
        t6 = time.time()
        print("Time for garbage collection: %s" % (str(t6-t5)))

    return background_insert_residues_process

@ray.remote(num_cpus=1)
def annotate(config,pdb_id,target_dict):
    #hack proposed by the devs of ray to prevent too many processes being spawned
    resources = ray.ray.get_resource_ids() 
    cpus = [v[0] for v in resources['CPU']]
    psutil.Process().cpu_affinity(cpus)

    t0 = time.time()
    (structural_analysis_dict,errorlist,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles) = templateFiltering.structuralAnalysis(pdb_id,config,target_dict = target_dict)
    t1 = time.time()

    if config.verbosity >= 3:
        print('Time for structural analysis of',pdb_id,':',t1-t0)

    return (pdb_id,structural_analysis_dict,errorlist,ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles,t1-t0)


def core(protein_list,indels,config,session,manager,lock,outfolder,session_name,out_objects):
    if len(protein_list) == 0:
        return

    background_insert_residues_process = None
    background_process_MS = None

    #transform the protein map into a Proteins object
    proteins = sdsc.Proteins(protein_list,indels, lite = config.lite)

    if config.verbosity >= 4:
        print('Proteins state after object initialization:')
        proteins.print_protein_state()

    try:
        t1 = time.time()
        if config.verbosity >= 1:
            print("Before protCheck")
        #check for already stored genes and make all the necessary database interactions
        #sets the fields database_id and stored in the protein objects as well as the stored and non stored ids in the proteins object

        if not config.lite:
            database.protCheck(proteins,session,config)

        if config.verbosity >= 4:
            print('Proteins state after protCheck:')
            proteins.print_protein_state()

        t2 = time.time()
        if config.verbosity >= 2:
            print("Time for protCheck: %s" % (str(t2-t1)))
        #check for already stored mutations or position twins and make all the necessary database interactions
        #sets the fields database_id and stored in the position objects

        if not config.lite:
            background_process_MS = database.positionCheck(proteins,session,config)

            database.indelCheck(proteins,session,config)

        if config.verbosity >= 4:
            print('Proteins state after positionCheck:')
            proteins.print_protein_state()

        t3 = time.time()
        if config.verbosity >= 2:
            print("Time for positionCheck: %s" % (str(t3-t2)))

        if config.verbosity >= 1:
            print("Before getSequences")

        background_iu_process = getSequences(proteins,config,manager,lock,skip_db = config.lite)

        t4 = time.time()
        if config.verbosity >= 2:
            print("Time for getSequences: %s" % (str(t4-t3)))

        if config.verbosity >= 1:
            print("Before autoTemplateSelection")
        autoTemplateSelection(config,manager,lock,proteins)

        t5 = time.time()
        if config.verbosity >= 2:
            print("Time for Template Selection: %s" % (str(t5-t4)))

        if config.verbosity >= 1:
            print("Before paraAlignment")
        paraAlignment(config,manager,lock,proteins,skip_db = config.lite)

        t6 = time.time()
        if config.verbosity >= 2:
            print("Time for Alignment: %s" % (str(t6-t5)))

        if background_iu_process != None:
            background_iu_process.join() #Disorder values have to finished before the classification happens

        if config.verbosity >= 1:
            print("Before paraAnnotate")

        background_insert_residues_process = paraAnnotate(config,manager,lock,proteins, lite = config.lite)

        t7 = time.time()
        if config.verbosity >= 2:
            print("Time for Annotation: %s" % (str(t7-t6)))

        if not config.lite:
            indel_analysis.para_indel_analysis(proteins,config,manager,lock)
        else:
            out_objects = output.appendOutput(proteins,outfolder,session_name,out_objects = out_objects)

        t1 = time.time()
        if config.verbosity >= 2:
            print("Time for Indelanalysis: %s" % (str(t1-t7)))

        #join the background inserts
        if background_process_MS != None:
            background_process_MS.join()
        if background_insert_residues_process != None:
            background_insert_residues_process.join()

        t2 = time.time()
        if config.verbosity >= 2:
            print('Resttime for background inserts: ',t2-t1)

        if config.verbosity >= 3:
            for name, size in sorted(((name, sys.getsizeof(value)) for name, value in locals().items()),
                     key= lambda x: -x[1])[:10]:
                print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))
            print(list_fds())

    #Error-Handling for a whole input line
    except:

        [e,f,g] = sys.exc_info()
        g = traceback.format_exc()
        #print "Pipeline Core Error: ",e,f,g
        errortext = '\n'.join([str(e),str(f),str(g)]) + '\n\n'
        config.errorlog.add_error(errortext)
        if background_process_MS != None:
            background_process_MS.join()
        if background_insert_residues_process != None:
            background_insert_residues_process.join()

    return out_objects

#@profile
def main(filename,config,output_path,main_file_path):
    n_of_cores = config.proc_n
    mrna_fasta = config.mrna_fasta
    num_of_cores = config.proc_n
    verbose = config.verbose
    session = 0 #This can later be used, structman.py could give a specific session id and the pipeline can then expand that session

    config.errorlog.start(filename,output_path)

    t0 = time.time()

    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/lib'
    os.environ["PYTHONPATH"] = parent_dir + ":" + os.environ.get("PYTHONPATH", "")

    ray.init(num_cpus = config.annotation_processes, include_dashboard = False, ignore_reinit_error=True)

    '''
    mem_tracked_processes = [] 
    for proc in psutil.process_iter(['pid', 'name', 'username']):
        if proc.info['username'] == 'agr18':
            p_id = proc.info['pid']
            mem_tracked_processes.append(psutil.Process())
    '''

    if mrna_fasta != None:
        if not os.path.exists(mrna_fasta):
            raise NameError("mRNA path not found: %s" % mrna_fasta)

    #annovar-pipeline in case of vcf-file
    if isinstance(filename,str):
        if filename.rsplit(".",1)[1] == "vcf":
            anno_db = "%s_annovar" % db_name.rsplit("_",1)[0]
            if config.verbosity >= 1:
                print('Convert vcf file format using Annovar')
            if mrna_fasta != None:
                '... and using mrna file: ',mrna_fasta
            nfname = annovar.annovar_pipeline(filename,config.tax_id,config.annovar_path,config.db_adress,config.db_user_name,config.db_password,anno_db,mrna_fasta,ref_id=config.ref_genome_id)
        else:
            nfname = filename
    else:
        nfname = 'Single line input'
        single_line_inputs = filename

    t01 = time.time()

    if config.verbosity >= 2:
        print("Time for preparation before buildQueue: %s" % (str(t01-t0)))

    chunksize = config.chunksize

    if config.verbosity >= 1:
        print("Call buildQueue with chunksize: %s and file: %s" % (str(chunksize),nfname))

    if nfname != 'Single line input':
        if config.fasta_input:
            proteins_chunks = parseFasta(config,nfname)
        else:
            proteins_chunks = buildQueue(config,nfname,chunksize)
    else:
        proteins_chunks = buildQueue(config,single_line_inputs,chunksize)
        nfname = '/%s.' % (' '.join(single_line_inputs))

    t02 = time.time()
    if config.verbosity >= 2:
        print("Time for buildQueue: %s" % (str(t02-t01)))
    if config.verbosity >= 1:
        print("Number of chunks: ",len(proteins_chunks))

    newsession = False
    if session == 0 and not config.lite:
        starttime = SQLDateTime()
        session = database.insertSession(starttime,nfname,config)
        newsession = True
    session_name = (nfname.rsplit("/",1)[1]).rsplit(".",1)[0]

    try:
        os.stat("%s/tmp_structman_pipeline" %(output_path))
    except:
        os.mkdir("%s/tmp_structman_pipeline" %(output_path))  
    os.chdir("%s/tmp_structman_pipeline" %(output_path))
    cwd = "%s/tmp_structman_pipeline" %(output_path)

    config.temp_folder = cwd

    manager = multiprocessing.Manager()
    lock = manager.Lock()

    out_objects = None

    chunk_nr = 1 
    for protein_list,indels in proteins_chunks:
        if config.verbosity >= 1:
            print("Chunk %s/%s" % (str(chunk_nr),str(len(proteins_chunks))))
        chunk_nr+=1

        out_objects = core(protein_list,indels,config,session,manager,lock,output_path,session_name,out_objects)

    os.chdir(output_path)

    if config.verbosity >= 2:
        t03 = time.time()

    try:
        shutil.rmtree(cwd)
    except:
        pass

    if config.verbosity >= 2:
        t04 = time.time()
        print('Time for folder cleanup:',t04-t03)

    config.errorlog.stop()

    if newsession:
        endtime = SQLDateTime()
        database.updateSession(session,endtime,config)

    tend = time.time()
    if config.verbosity >= 2:
        print((tend-t0))

    '''
    total_memory_peak = 0
    for p in mem_tracked_processes:
        total_memory_peak += p.memory_info().rss

    total_memory_peak = total_memory_peak/1024./1024./1024.

    print('Accumulated memory peak',total_memory_peak,'Gb')
    '''

    return session
