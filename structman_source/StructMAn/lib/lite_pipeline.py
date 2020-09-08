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

import serializedPipeline
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

def appendOutput(proteins,outfile):
    lines = []

    u_acs = proteins.get_protein_u_acs()

    for u_ac in u_acs:
        u_id = proteins.get_u_id(u_ac)
        refseqs = proteins.get_ref_ids(u_ac)
        aaclist = proteins.getAACList(u_ac)
        gpan = ';'.join(refseqs)
        for aac_base in aaclist:
            aa2 = aaclist[aac_base]
            pos = int(aac_base[1:])
            
            m = (u_ac,pos)
            aa1 = aac_base[0]

            aac = '%s%s' % (aac_base,aa2)
            tag = proteins.get_mut_tags(u_ac,pos,aa2)

            position = proteins.get_position(u_ac,pos)
            if position == None:
                continue
            mappings = position.mappings

            Class = mappings.Class
            conf = mappings.classification_conf
            weighted_sc = mappings.weighted_location
            recommended_structure_str = mappings.get_recommended_res_str(proteins)
            recommended_structure,seq_id,cov,resolution = database.process_recommend_structure_str(recommended_structure_str)
            max_seq_structure_str = mappings.get_max_seq_structure_res_str(proteins)
            max_seq_structure,max_seq_seq_id,max_seq_cov,max_seq_resolution = database.process_recommend_structure_str(max_seq_structure_str)

            amount_of_structures = len(mappings.qualities)
            mv_sec_ass = mappings.weighted_ssa
            simple_class = mappings.simple_class
            interaction_str = str(mappings.interaction_recommendations)

            input_res_id = proteins.get_res_id(u_ac,pos)
            input_pdb_id = ''
            
            if input_res_id != None:
                input_pdb_id = u_ac

            lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (u_ac,u_id,gpan,input_pdb_id,
                        input_res_id,aa1,pos,u_ac,tag,weighted_sc,Class,simple_class,interaction_str,str(conf),mv_sec_ass,
                        recommended_structure,seq_id,cov,resolution,max_seq_structure,max_seq_seq_id,max_seq_cov,
                        max_seq_resolution,str(amount_of_structures)))

    f = open(outfile,'a')
    f.write('\n'.join(lines) + '\n')
    f.close()
    
    return

def main(filename,config,output_path,main_file_path):
    n_of_cores = config.proc_n

    mrna_fasta = config.mrna_fasta
    num_of_cores = config.proc_n
    verbose = config.verbose
    search_tool = config.search_tool
    errorlog = config.errorlog

    manager = multiprocessing.Manager()
    lock = manager.Lock()

    t0 = time.time()

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
    proteins_chunks = serializedPipeline.buildQueue(config,nfname,junksize,mrna_fasta=mrna_fasta)

    t02 = time.time()
    if verbose:
        print("Time for buildQueue: %s" % (str(t02-t01)))
        print("Number of chunks: ",len(proteins_chunks))

    errorlog.start(nfname,session_name)

    junk_nr = 1 
    for proteins in proteins_chunks:
        if len(proteins) == 0:
            continue

        print("Chunk %s/%s" % (str(junk_nr),str(len(proteins_chunks))))
        junk_nr+=1
        try:
            os.stat("%s/tmp_structman_pipeline" %(output_path))
        except:
            os.mkdir("%s/tmp_structman_pipeline" %(output_path))  
        os.chdir("%s/tmp_structman_pipeline" %(output_path))
        cwd = "%s/tmp_structman_pipeline" %(output_path)

        try:
            #transform the protein map into a Proteins object
            proteins = sdsc.Proteins(proteins, lite = True)

            if verbose:
                t3 = time.time()
                print("Before getSequences")

            serializedPipeline.getSequences(proteins,config,manager,lock,skip_db = True)

            if verbose:
                t4 = time.time()
                print("Time for getSequences: %s" % (str(t4-t3)))
                print("Before autoTemplateSelection")

            serializedPipeline.autoTemplateSelection(config,manager,lock,proteins,search_tool=search_tool)


            if verbose:
                t5 = time.time()
                print("Time for Template Selection: %s" % (str(t5-t4)))
                print("Before paraAlignment")

            serializedPipeline.paraAlignment(config,manager,lock,proteins, skip_db = True)

            if verbose:
                t6 = time.time()
                print("Time for Alignment: %s" % (str(t6-t5)))
                print("Before paraAnnotate")

            serializedPipeline.paraAnnotate(config,manager,lock,proteins, lite = True)
            t7 = time.time()
            if verbose:
                print("Time for Annotation: %s" % (str(t7-t6)))

            appendOutput(proteins,outfile)

        #Error-Handling for a whole input line
        except:
            errorlog.add_error('Litepipeline main loop error')


        os.chdir(output_path)
        #"""
        try:
            shutil.rmtree(cwd)
        except:
            pass
        #"""

    errorlog.stop()

    
    tend = time.time()
    print((tend-t0))
    return
