#!/usr/bin/python3
#This is the wrapper script, including and automatizing all functions StructMAn has to offer
import sys
import getopt
import os
import time
import pymysql as MySQLdb
import multiprocessing

main_file_path = os.path.abspath(sys.argv[0])

sys.path.append("%s/lib" % main_file_path.rsplit('/',1)[0])

import output
import database
import serializedPipeline
import lite_pipeline


cwd = os.getcwd()

class Config:
    def __init__(self,config_path,infile,num_of_cores,output_path,overwrite,annotationtable):
        self.db_adress = ""
        self.db_user_name = ""
        self.db_password = ""
        self.db_name = ""
        self.config = ''
        self.infile = ''
        self.outfolder = output_path
        self.overwrite = False
        self.annotationtable = annotationtable
        self.num_of_cores = None
        self.neighborhood_calculation = False
        self.error_annotations_into_db = True
        self.anno_session_mapping = True
        self.calculate_interaction_profiles=False
        self.verbose = False

        self.resources = 'manu'
        self.num_of_cores = None
        self.proc_n = 48
        self.blast_processes = self.proc_n
        self.alignment_processes = self.proc_n
        self.annotation_processes = self.proc_n
        self.number_of_processes = self.proc_n

        self.dssp = True
        self.pdb_input_asymetric_unit = False
        self.search_tool='MMseqs2'

        #active options
        self.option_seq_thresh = 35.0
        self.option_res_thresh = 4.5
        self.option_ral_thresh = 0.5
        self.option_seq_wf = 1.0
        self.option_ral_wf = 0.5
        self.option_res_wf = 0.25
        self.option_rval_wf = 0.1
        self.option_lig_wf = 1.0
        self.option_chain_wf = 1.0
        self.option_number_of_templates = 0
        self.tax_id = None
        self.ref_genome_id = None
        self.mrna_fasta = None
        self.species = None

        self.auto_scale = False
        self.auto_scale_template_threshold = 10
        #auto_scale_template_percentage gets set by buildQueue
        self.auto_scale_template_percentage = None

        #Paths to Helper Programs     
        self.blast_path = ""

        self.mmseqs2_path = ""
        
        self.output_path = ""
        self.pdb_path = ""
        self.annovar_path = ""

        self.dssp_path = ''
        self.rin_db_path = ''
        self.iupred_path = ''
        self.errorlog_path = ""

        self.go = False
        self.anno = False
        self.classification=True
        self.gene = False
        self.path = False
        self.godiff = False
        self.pathdiff = False
        self.do_modelling = False
        self.multi_modelling = False
        self.ligand_file = None
        self.mod_per_mut = 0
        self.mod_per_gene = 0
        self.tanimoto_cutoff = 0.05
        self.distance_threshold = 10.0
        self.ligand_filter = None
        self.proteome = False
        self.intertable_conf=False

        trunk = main_file_path.rsplit('/',1)[0]

        self.base_path = None

        self.blast_db_path = '%s/lib/base/blast_db/pdbba' % trunk
        self.mmseqs2_db_path = '%s/lib/base/blast_db/pdbba_search_db_mmseqs2' % trunk
        self.smiles_path = '%s/lib/base/ligand_bases/Components-smiles-stereo-oe.smi' % trunk
        self.inchi_path = '%s/lib/base/ligand_bases/inchi_base.tsv' % trunk
        self.human_id_mapping_path = '%s/lib/base/id_mapping' % trunk
        self.mmseqs_tmp_folder = '%s/lib/base/blast_db/tmp' % trunk


        self.lite = False

        f = open(config_path,'r')
        lines = f.read().split('\n')
        f.close()

        for line in lines:
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            
            words = line.split('=')
            if len(words) < 2:
                continue
            if words[0] == 'db_adress':
                self.db_adress = words[1]
            if words[0] == 'db_user_name':
                self.db_user_name = words[1]
            if words[0] == 'db_password':
                self.db_password = words[1]
            if words[0] == 'db_name':
                self.db_name = words[1]
            if words[0] == 'do_anno':
                if words[1] == 'True':
                    self.annotationtable=True

            opt = words[0]
            arg = words[1].replace("\n","")
            if opt == 'blast_path':
                self.blast_path = arg
            elif opt == 'pdb':
                self.pdb_path = arg
            elif opt == 'annovar_path':
                self.annovar_path = arg
            elif opt == 'dssp_path':
                self.dssp_path = arg
            elif opt == 'rin_db_path':
                self.rin_db_path = arg
            elif opt == 'iupred_path':
                self.iupred_path = arg
            elif opt == 'mmseqs2_path':
                self.mmseqs2_path = arg

            elif opt == 'base_path':
                self.base_path = arg
            elif opt == 'mmseqs_tmp_folder':
                self.mmseqs_tmp_folder = arg

            elif opt == 'seq_thresh':
                self.option_seq_thresh = float(arg)
                if self.option_seq_thresh <= 1.0:
                    self.option_seq_thresh *= 100.0
            elif opt == 'res_thresh':
                self.option_res_thresh = float(arg)
            elif opt == 'cov_thresh':
                self.option_ral_thresh = float(arg)
                if self.option_ral_thresh > 1.0:
                    self.option_ral_thresh *= 0.01
            elif opt == 'seq_wf':
                self.option_seq_wf = float(arg)
            elif opt == 'cov_wf':
                self.option_ral_wf = float(arg)
            elif opt == 'res_wf':
                self.option_res_wf = float(arg)
            elif opt == 'rval_wf':
                self.option_rval_wf = float(arg)
            elif opt == 'lig_wf':
                self.option_lig_wf = float(arg)
            elif opt == 'chain_wf':
                self.option_chain_wf = float(arg)
            elif opt == 'option_number_of_templates':
                self.option_number_of_templates = int(arg)
            elif opt == 'db_adress':
                db_adress = arg
            elif opt == 'db_user_name':
                self.db_user_name = arg
            elif opt == 'db_password':
                self.db_password = arg
            elif opt == 'db_name':
                self.db_name = arg
            elif opt == 'species':
                self.species = arg
            elif opt == 'tax_id':
                self.tax_id = arg
            elif opt == 'ref':
                self.ref_genome_id = arg
            elif opt == 'resources':
                self.resources = arg
            elif opt == 'auto_scale':
                if arg == 'on':
                    self.auto_scale = True
                else:
                    self.auto_scale = False
            elif opt == 'error_annotations_into_db':
                if arg == 'True':
                    self.error_annotations_into_db = True
                elif arg == 'False':
                    self.error_annotations_into_db = False
            elif opt == 'anno_session_mapping':
                if arg == 'True':
                    self.anno_session_mapping = True
                elif arg == 'False':
                    self.anno_session_mapping = False
            elif opt == 'mrna':
                self.mrna_fasta = arg
            elif opt == 'neighborhood_calculation':
                if arg == 'True':
                    self.neighborhood_calculation = True
            elif opt == 'calculate_interaction_profiles':
                if arg == 'True':
                    self.calculate_interaction_profiles = True
                elif arg == 'False':
                    self.calculate_interaction_profiles = False
            elif opt == 'pdb_input_asymetric_unit':
                if arg == 'True':
                    self.pdb_input_asymetric_unit = True
                elif arg == 'False':
                    self.pdb_input_asymetric_unit = False
            elif opt == 'verbose':
                if arg == 'True':
                    self.verbose = True
                elif arg == 'False':
                    self.verbose = False

            elif opt == 'do_anno':
                if arg == "True":
                    self.anno = True
                else:
                    self.anno = False
            elif opt == 'do_genesort':
                if arg == "True":
                    self.gene = True
                else:
                    self.gene = False
            elif opt == 'do_goterm':
                if arg == "True":
                    self.go = True
                else:
                    self.go = False
            elif opt == 'do_godiff':
                if arg == "True":
                    self.godiff = True
                else:
                    self.godiff = False
            elif opt == 'do_pathway':
                if arg == "True":
                    self.path = True
                else:
                    self.path = False
            elif opt == 'do_pathdiff':
                if arg == "True":
                    self.pathdiff = True
                else:
                    self.pathdiff = False
            elif opt == 'multi_modelling':
                if arg == "True":
                    self.multi_modelling = True
                else:
                    self.multi_modelling = False
            elif opt == 'mod_per_gene':
                self.mod_per_gene = arg
            elif opt == 'mod_per_mut':
                self.mod_per_mut = arg
            elif opt == 'do_modelling':
                if arg == "True":
                    self.do_modelling = True
                else:
                    self.do_modelling = False
            
            elif opt == 'tanimoto_cutoff':
                self.tanimoto_cutoff = float(arg)
            elif opt == 'lig_dist_thresh':
                self.distance_threshold = arg
            elif opt == 'ligand_filter':
                self.ligand_filter = arg
            elif opt == 'proteome':
                if arg == 'True':
                    self.proteome = True
            elif opt == 'intertable':
                if arg == 'True':
                    self.intertable_conf = True
            elif opt == 'lite':
                if arg == 'True':
                    self.lite = True
        
        if self.base_path != None:
            print('Using structman_data from :',self.base_path)
            self.blast_db_path = '%s/base/blast_db/pdbba' % self.base_path
            self.mmseqs2_db_path = '%s/base/blast_db/pdbba_search_db_mmseqs2' % self.base_path
            self.smiles_path = '%s/base/ligand_bases/Components-smiles-stereo-oe.smi' % self.base_path
            self.inchi_path = '%s/base/ligand_bases/inchi_base.tsv' % self.base_path
            self.human_id_mapping_path = '%s/base/id_mapping' % self.base_path

        if self.resources == 'auto' and self.num_of_cores == None:
            self.proc_n = multiprocessing.cpu_count() -1
            if self.proc_n <= 0:
                self.proc_n = 1

            self.blast_processes = self.proc_n
            self.alignment_processes = self.proc_n
            self.annotation_processes = self.proc_n
            self.number_of_processes = self.proc_n

        if self.num_of_cores != None:
            self.proc_n = self.num_of_cores
            
            self.blast_processes = self.proc_n
            self.alignment_processes = self.proc_n
            self.annotation_processes = self.proc_n
            self.number_of_processes = self.proc_n

        print('Using %s core(s)' % str(self.proc_n))
        if not os.path.exists("%s/errorlogs" % self.outfolder):
            os.mkdir("%s/errorlogs" % self.outfolder)
        self.errorlog_path = "%s/errorlogs/errorlog.txt" % self.outfolder



def main(infiles,out_folder,main_file_path,config,overwrite=False,anno=False,intertable=False):
    db_name = config.db_name
    db_adress = config.db_adress
    db_password = config.db_password
    db_user_name = config.db_user_name
    for infile in infiles:
        print('Processing file: ',infile)
        if not config.lite:
            #check if the infile got already processed
            db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
            cursor = db.cursor()
            session_id = database.getSessionId(infile,db,cursor)
            db.close()

        #create the output folder
        [trunk,infilename] = infile.rsplit('/',1)
        if out_folder == None:
            out_folder = "%s/Output" % trunk
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)

        if not config.lite:
            #run the main pipeline
            if session_id == None:
                session_id = serializedPipeline.main(infile,config,out_folder,main_file_path)

            #run the output scripts
            session_name = infilename.rsplit('.',1)[0]
            months = {1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}
            time_struct = time.gmtime()
            year = str(time_struct[0])
            month = months[time_struct[1]]
            day = str(time_struct[2])
            date = "%s%s%s" % (day,month,year)
            outpath = "%s/%s_%s" % (out_folder,session_name,date)
            if not os.path.exists(outpath):
                os.mkdir(outpath)

            output.main(session_id,outpath,config,overwrite=overwrite,anno=anno,intertable=intertable,)
        else:
            lite_pipeline.main(infile,config,out_folder,main_file_path)


if __name__ == "__main__":

    multiprocessing.freeze_support()
    multiprocessing.set_start_method('forkserver')
    #print(multiprocessing.get_start_method())
    disclaimer = 'Usage: structman.py <-i input_file>\n\n\n##### Optional parameter: #####\n\n<-n threads> : Number of cores to be used\n\n<-o output_folder> : Path to the output folder\n\n<-c config_file> : Path to the configuration file\n\n<-l> : lite-mode\n\n<-v> : verbose output'

    argv = sys.argv[1:]
    try:
        opts,args = getopt.getopt(argv,"c:i:n:o:h:lv",['help','overwrite','annotationtable'])
    except getopt.GetoptError:
        print("Illegal Input\n\n",disclaimer)
        sys.exit(2)

    annotationtable = False
    infile = ''
    config_path = ''
    num_of_cores = None
    outfolder = ''
    overwrite = False
    lite = False
    verbose_flag = False

    for opt,arg in opts:
        if opt == '-c':
            config_path = arg
        if opt == '-v':
            verbose_flag = True
        if opt == '-l':
            lite = True
        if opt == '-i':
            infile = arg
        if opt == '-n':
            num_of_cores = int(arg)
        if opt == '-o':
            outfolder = arg
        if opt == '--overwrite':
            overwrite = True
        if opt == '--annotationtable':
            annotationtable = True
        if opt == '-h' or opt == '--help':
            print(disclaimer)
            sys.exit(0)

    if infile == '':
        input_folder = '/structman/input_data/'
        if not os.path.isdir(input_folder):
            print('Did not find the input path\n\n',disclaimer)
            sys.exit(2)
        filelist = os.listdir(input_folder)
        
        infiles = []
        for infile in filelist:
            if infile.split('.')[-1] == 'smlf' or infile.split('.')[-1] == 'vcf':
                infiles.append('%s/%s' % (input_folder,infile))
        if infiles == []:
            print('Did not find any file inside the input folder\n\n',disclaimer)
            sys.exit(2)

    else:
        #make infile to global path
        infile = os.path.abspath(infile)
        infiles = [infile]


    if outfolder == '':
        outfolder = '/structman/results/'
        if not os.path.isdir(outfolder):
            if infile != '':
                outfolder = '%s/Output' % infile.rsplit('/',1)[0]
                if not os.path.exists(outfolder):
                    os.mkdir(outfolder)
            else:
                outfolder = None
        


    if config_path == '':
        if os.path.exists('%s/config.txt' % infile.rsplit('/',1)[0]):
            config_path = '%s/config.txt' % infile.rsplit('/',1)[0]
        elif os.path.exists('%s/config.txt' % cwd):
            config_path = '%s/config.txt' % cwd
        elif os.path.exists('%s/config.txt' % main_file_path.rsplit('/',1)[0]):
            config_path = '%s/config.txt' % main_file_path.rsplit('/',1)[0]
        else:
            print("No config file found, please use -c [Path to config]")
            sys.exit(2)

    print(("Using following config file: %s" % config_path))
    config_path = os.path.abspath(config_path)

    config = Config(config_path,infile,num_of_cores,outfolder,overwrite,annotationtable)

    if verbose_flag:
        config.verbose = True

    if lite:
        config.lite = True

    main(infiles,outfolder,main_file_path,config,overwrite=overwrite,anno=annotationtable)
