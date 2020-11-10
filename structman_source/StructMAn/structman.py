#!/usr/bin/python3
#This is the wrapper script, including and automatizing all functions StructMAn has to offer
import configparser
import sys
import getopt
import os
import time
import pymysql as MySQLdb
import multiprocessing
import subprocess
import traceback

main_file_path = os.path.abspath(sys.argv[0])

sys.path.append("%s/lib" % main_file_path.rsplit('/',1)[0])

import output
import database
import serializedPipeline
import lite_pipeline
import resource
import repairDB


class Config:
    def __init__(self,config_path ,num_of_cores = 1,output_path = '',
                    util_mode = False, output_util = False,external_call = True,profiling = False, verbosity = None,
                    print_all_errors = False,chunksize = 500):
        self.prog_start_time = time.time()
        # read config file, auto add section header so old config files work
        self.config_parser_obj = configparser.ConfigParser()
        try:
            self.config_parser_obj.read(config_path)
        except configparser.MissingSectionHeaderError:
            with open(config_path, 'r') as f:
                self.config_parser_obj.read_string('[user]\n%s' % f.read())
        cfg = self.config_parser_obj['user']

        self.db_adress = cfg.get('db_adress', fallback='')
        self.db_user_name = cfg.get('db_user_name', fallback='')
        self.db_password = cfg.get('db_password', fallback='')
        self.db_name = cfg.get('db_name', fallback='')
        self.mapping_db = cfg.get('mapping_db', fallback=None)

        self.outfolder = output_path
        self.num_of_cores = num_of_cores

        self.neighborhood_calculation = cfg.getboolean('neighborhood_calculation', fallback=True)
        self.error_annotations_into_db = cfg.getboolean('error_annotations_into_db', fallback=True)
        self.anno_session_mapping = cfg.getboolean('anno_session_mapping', fallback=True)
        self.calculate_interaction_profiles = cfg.getboolean('calculate_interaction_profiles', fallback=True)

        self.verbose = cfg.getboolean('verbose', fallback=False)
        self.verbosity = cfg.getint('verbosity', fallback=0)

        self.profiling = profiling
        self.skipref = cfg.getboolean('skipref', fallback=False)

        self.resources = cfg.get('resources', fallback='manu')

        self.proc_n = 48
        self.blast_processes = self.proc_n
        self.alignment_processes = self.proc_n
        self.annotation_processes = self.proc_n
        self.number_of_processes = self.proc_n

        self.dssp = cfg.getboolean('dssp', fallback=True)
        self.pdb_input_asymetric_unit = cfg.getboolean('pdb_input_asymetric_unit', fallback=False)
        self.search_tool = cfg.get('search_tool', fallback='MMseqs2')

       #active options
        self.option_seq_thresh = cfg.getfloat('seq_thresh', fallback=35.0)
        if self.option_seq_thresh <= 1.0:
            self.option_seq_thresh *= 100.0
        self.option_res_thresh = cfg.getfloat('res_thresh', fallback=4.5)
        self.option_ral_thresh = cfg.getfloat('cov_thresh', fallback=0.5)
        if self.option_ral_thresh > 1.0:
            self.option_ral_thresh *= 0.01
        self.option_seq_wf = cfg.getfloat('seq_wf', fallback=1.0)
        self.option_ral_wf = cfg.getfloat('cov_wf', fallback=0.5)
        self.option_res_wf = cfg.getfloat('res_wf', fallback=0.25)
        self.option_rval_wf = cfg.getfloat('rval_wf', fallback=0.1)
        self.option_lig_wf = cfg.getfloat('lig_wf', fallback=1.0)
        self.option_chain_wf = cfg.getfloat('chain_wf', fallback=1.0)
        self.option_number_of_templates = cfg.getint('option_number_of_templates', fallback=0)
        self.tax_id = cfg.get('tax_id', fallback=None)
        self.ref_genome_id = cfg.get('ref', fallback=None)
        self.mrna_fasta = cfg.get('mrna', fallback=None)

        self.surface_threshold = cfg.getfloat('surface_threshold', fallback=0.16)
        self.short_distance_threshold = cfg.getfloat('lig_short_dist_thresh', fallback=5.0)
        self.long_distance_threshold = cfg.getfloat('lig_long_dist_thresh', fallback=8.0)
        self.ligand_interest_sphere = cfg.getfloat('ligand_interest_sphere', fallback=25.0)

        #Paths to Helper Programs
        self.blast_path = cfg.get('blast_path', fallback='')
        self.mmseqs2_path = cfg.get('mmseqs2_path', fallback='')
        self.output_path = cfg.get('output_path', fallback='')
        self.pdb_path = cfg.get('pdb_path', fallback='')

        self.annovar_path = cfg.get('annovar_path', fallback='')
        self.dssp_path = cfg.get('dssp_path', fallback='')
        self.rin_db_path = cfg.get('rin_db_path', fallback='')
        self.iupred_path = cfg.get('iupred_path', fallback='')
        self.errorlog_path = cfg.get('errorlog_path', fallback=None)
        self.base_path = cfg.get('base_path', fallback=None)

        self.go = cfg.getboolean('do_goterm', fallback=False)
        self.anno = cfg.getboolean('do_anno', fallback=False)
        self.classification=cfg.getboolean('do_classification', fallback=True)
        self.gene = cfg.getboolean('do_genesort', fallback=False)
        self.path = cfg.getboolean('do_pathway', fallback=False)
        self.godiff = cfg.getboolean('do_godiff', fallback=False)
        self.pathdiff = cfg.getboolean('do_pathdiff', fallback=False)
        self.do_modelling = cfg.getboolean('do_modelling', fallback=False)
        self.multi_modelling = cfg.getboolean('multi_modelling', fallback=False)
        self.ligand_file = cfg.get('ligand_file', fallback=None)
        self.mod_per_mut = cfg.getint('mod_per_mut', fallback=0)
        self.mod_per_gene = cfg.getint('mod_per_gene', fallback=0)
        self.tanimoto_cutoff = cfg.getfloat('tanimoto_cutoff', fallback=0.05)
        self.milieu_threshold = cfg.getfloat('milieu_threshold', fallback=10.0)
        self.ligand_filter = cfg.get('ligand_filter', fallback=None)
        self.proteome = cfg.getboolean('proteome', fallback=False)
        self.intertable_conf = cfg.getboolean('intertable', fallback=False)
        self.overwrite_incorrect_wt_aa = cfg.getboolean('overwrite_incorrect_wt_aa', fallback=False)

        trunk = os.path.realpath(__file__).rsplit('/',1)[0]

        self.database_source_path = '%s/struct_man_db.sql' % trunk
        self.blast_db_path = '%s/lib/base/blast_db/pdbba' % trunk
        self.mmseqs2_db_path = '%s/lib/base/blast_db/pdbba_search_db_mmseqs2' % trunk
        self.smiles_path = '%s/lib/base/ligand_bases/Components-smiles-stereo-oe.smi' % trunk
        self.inchi_path = '%s/lib/base/ligand_bases/inchi_base.tsv' % trunk
        self.human_id_mapping_path = '%s/lib/base/id_mapping' % trunk
        self.rinerator_base_path = '%s/lib/rinerator' % trunk
        self.rinerator_path = '%s/lib/rinerator/get_chains.py' % trunk
        os.environ["REDUCE_HET_DICT"] = '%s/lib/rinerator/reduce_wwPDB_het_dict.txt' % trunk

        self.pdb_sync_script = '%s/pdb-rsync.sh' % trunk

        self.mmseqs_tmp_folder = cfg.get('mmseqs_tmp_folder', '%s/lib/base/blast_db/tmp' % trunk)

        self.fasta_input = cfg.getboolean('fasta_input', fallback=False)
        self.lite = cfg.getboolean('lite', fallback=False)

        if verbosity != None:
            self.verbosity = verbosity

        self.chunksize = chunksize
        # Checking whether the given paths in config exist or not if it is given and not exist, system gives error message and exits

        if self.search_tool=='MMseqs2':
            if self.mmseqs2_path:
                if self.mmseqs2_path == 'mmseqs':
                    cmd = "mmseqs"
                    ret_val = subprocess.getstatusoutput(cmd)[0]
                    #ret_val = os.system(cmd)
                    if ret_val != 0:
                        print('MMSEQS2 path is not correct, please correct the path')
                        sys.exit()
                else:
                    isExist = os.path.exists(self.mmseqs2_path)
                    if not isExist:
                        print('MMSEQS2 path does not exist, please check the path')
                        sys.exit()
        elif self.search_tool=='Blast':
            isExist = os.path.exists(self.blast_path)
            if not isExist:
                print('Blast path does not exist, please check the path or use mmseqs')
                sys.exit()
        if self.pdb_path:
            isExist = os.path.exists(self.pdb_path)
            if not isExist:
                print('PDB path does not exist, please check the path')
                sys.exit()
        if self.annovar_path:
            isExist = os.path.exists(self.annovar_path)
            if not isExist:
                print('Annovar path does not exist, please check the path')
                sys.exit()
        if self.dssp_path:
            if self.dssp_path == 'mkdssp':
                cmd = "mkdssp"
                ret_val = subprocess.getstatusoutput(cmd)[0]
                #ret_val = os.system(cmd)
                if ret_val != 1:
                    print('DSSP path is not correct, please correct the path')
                    sys.exit()
            else:
                isExist = os.path.exists(self.dssp_path)
                if not isExist:
                    print('DSSP path does not exist, please check the path')
                    sys.exit()
        if self.rin_db_path:
            isExist = os.path.exists(self.rin_db_path)
            if not isExist:
                print('RIN DB path does not exist, please check the path')
                sys.exit()
        if self.iupred_path:
            isExist = os.path.exists(self.iupred_path)
            if not isExist:
                print('IUPred path does not exist, please check the path')
                sys.exit()
        if self.base_path:
            isExist = os.path.exists(self.base_path)
            if not isExist:
                print('Base path does not exist, please check the path')
                sys.exit()
        if self.mmseqs_tmp_folder:
            isExist = os.path.exists(self.mmseqs_tmp_folder)
            if not isExist:
                print('MMSEQS temp folder path does not exist, please check the path')
                sys.exit()
            if not os.access(self.mmseqs_tmp_folder,os.R_OK):
                print('Need writing rights in the MMSEQS temp folder, please check the path')
                sys.exit()

        if self.base_path != None and not util_mode:
            if self.verbosity >= 1:
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

        if self.proc_n > multiprocessing.cpu_count():
            if self.verbosity >= 1:
                print('More processes annotated (',self.proc_n,') than cores registered in system (',multiprocessing.cpu_count(),').')
            self.proc_n = multiprocessing.cpu_count()
            self.blast_processes = self.proc_n
            self.alignment_processes = self.proc_n
            self.annotation_processes = self.proc_n
            self.number_of_processes = self.proc_n

        if not util_mode:
            if not external_call and not os.path.exists(self.outfolder):
                os.mkdir(self.outfolder)
            if self.verbosity >= 1:
                print('Using %s core(s)' % str(self.proc_n))
            if (not external_call) and (not print_all_errors): #no need for errorlogs, when the config is generated not from the main script
                if not os.path.exists("%s/errorlogs" % self.outfolder):
                    os.mkdir("%s/errorlogs" % self.outfolder)
                self.errorlog_path = "%s/errorlogs/errorlog.txt" % self.outfolder

        self.errorlog = Errorlog(path = self.errorlog_path)

        #Determine maximal package size from database
        try:
            db,cursor = self.getDB(server_connection = True)
        except:
            if self.verbosity >= 2:
                print('Database connection failed in config initialization')
            db = None

        if db != None:
            cursor.execute("SHOW VARIABLES WHERE variable_name = 'max_allowed_packet'")
            self.max_package_size = int(cursor.fetchone()[1])*100//99
            db.close()
        else:
            self.max_package_size = None

    def getDB(self,server_connection = False,mapping_db = False):
        if server_connection:
            db = MySQLdb.connect(self.db_adress,self.db_user_name,self.db_password,None)
            cursor = db.cursor()
        elif mapping_db:
            db = MySQLdb.connect(self.db_adress,self.db_user_name,self.db_password,self.mapping_db)
            cursor = db.cursor()
        else:
            db = MySQLdb.connect(self.db_adress,self.db_user_name,self.db_password,self.db_name)
            cursor = db.cursor()
        return db,cursor

class Errorlog:
    def __init__(self,path = None):
        self.path = path
        self.error_counter = 0
        self.warning_counter = 0

    def start(self,nfname,session):
        if self.path == None:
            return
        date = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
        errortext  = "###############################################################################\n%s:\n%s - Session: %s\n" % (date,nfname,str(session))
        f = open(self.path,'a')
        f.write(errortext)
        f.close()
        return

    def add_error(self,error_text):
        self.error_counter += 1
        g = ''.join(traceback.format_list(traceback.extract_stack()))
        error_text = 'Error %s:\n%s\n%s\n' % (str(self.error_counter),error_text,g)
        if self.path == None:
            print(error_text)
            return
        f = open(self.path,'a')
        f.write(error_text)
        f.close()
        return

    def stop(self):
        if self.path == None:
            return
        if self.error_counter == 0:
            errortext  = "Finished without any error\n###############################################################################\n"
            f = open(self.path,'a')
            f.write(errortext)
            f.close()
        else:
            errortext  = "###############################################################################\n"
            f = open(self.path,'a')
            f.write(errortext)
            f.close()
            print("\n\nAt least one error occured, please check the errorlog.\n\n")
        return

    def add_warning(self,warn_text,lock = None):
        self.warning_counter += 1
        g = ''.join(traceback.format_list(traceback.extract_stack()))
        warn_text = 'Warning %s:\n%s\n%s\n' % (str(self.warning_counter),warn_text,g)
        if self.path == None:
            print(warn_text)
            return
        if lock != None:
            with lock:
                f = open(self.path,'a')
                f.write(warn_text)
                f.close()
            return
        f = open(self.path,'a')
        f.write(warn_text)
        f.close()
        return

def main(infiles,out_folder,main_file_path,config,intertable=False):
    db_name = config.db_name
    db_adress = config.db_adress
    db_password = config.db_password
    db_user_name = config.db_user_name
    #resource.setrlimit(resource.RLIMIT_AS,(6442450944,8589934592))
    #resource.setrlimit(resource.RLIMIT_AS,(442450944,589934592))
    for infile in infiles:
        if isinstance(infile,str):
            if config.verbosity >= 1:
                print('Processing file: ',infile)
            if not config.lite:
                #check if the infile got already processed
                session_id = database.getSessionId(infile,config)

            #create the output folder
            [trunk,infilename] = infile.rsplit('/',1)

            if infilename[-6:] == '.fasta':
                config.fasta_input = True

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
            outpath = "%s/%s" % (out_folder,session_name)
            if not os.path.exists(outpath):
                os.mkdir(outpath)

            output.main(session_id,outpath,config,intertable=intertable,)
        else:
            session_id = serializedPipeline.main(infile,config,out_folder,main_file_path)


if __name__ == "__main__":
    main_file_path = os.path.realpath(__file__)

    sys.path.append("%s/lib" % main_file_path.rsplit('/',1)[0])

    import output
    import database
    import serializedPipeline
    import lite_pipeline


    multiprocessing.freeze_support()
    multiprocessing.set_start_method('forkserver')
    #print(multiprocessing.get_start_method())
    disclaimer = ''.join([
                'Usage: structman.py <-i input_file>\n',
                'more functionalities can be used giving a second key word:\n',
                'structman.py database    gives you more info about the database utility functions\n',
                'structman.py config      gives you more info about the config utility functions\n',
                'structman.py update      gives you more info about the different update options\n\n\n',
                '##### Optional parameter: #####\n\n',
                '<-n threads> :           Number of cores to be used\n\n',
                '<-o output_folder> :     Path to the output folder\n\n',
                '<-c config_file> :       Path to the configuration file\n\n',
                '<-d> :                   database mode\n\n',
                '<-l> :                   lite-mode\n\n',
                '<--verbosity> [1-4]:     verbose output\n\n',
                '<--norin>                disable all RIN-based calculation'])

    database_util_disclaimer = ''.join([
                'Usage: structman.py database [command] <-c config_file>\n\n',
                '#### Commands: ####\n\n',
                'reset :     deletes all content of the database\n\n',
                'destroy :   completely removes the database\n\n',
                'create :    creates an empty instance of the database, usually called after database destroy'
                ])

    update_util_disclaimer = ''.join([
                'Usage: structman.py update [commands] <-c config_file> <-p path_to_local_pdb>\n',
                'The update functionalities manage local database, that significantly increases the performance of StructMAn.\n',
                'They need a lot of disk space though. Thus the usage of local databases is only recommended when one wants to process larger inputs with StructMAn.\n'
                'Only needs the path to the local instance of the PDB (-p), if it is not already specified in the config file.\n\n'
                '#### Commands: ####\n\n',

                'pdb :       uses rsync to update all files of the local PDB that could potentially used by StructMAn.\n'
                '            If the given path to the local PDB is empty, a large portion of the PDB will be downloaded, be sure that there is enought disk space available.\n\n',

                'rindb :     calculates and stores the RIN for each PDB structure in the local PDB.\n',
                '            When called for the first time, this can take multiple hours to complete. Also requires a large amount of disk space.'
                ])

    config_util_disclaimer = ''.join([
                'Usage: structman.py config [command] [value] <-c config_file>\n',
                'The config commands enable the expansion of StructMAn by giving it access to additional functionalities, which are too large to be included by default or require an external license.\n\n',
                '#### Commands: ####\n\n',
                'set_local_pdb_path <path_to_local_pdb> :               enables the usage of a local instance of the PDB. If this is not available, StructMAn has to download all structural information from the web.\n\n',
                'set_local_iupred_path <path_to_iupred_executable> :    enables StructMAn to use the predicted disordered regions performed by iupred2a.'
                ])

    argv = sys.argv[1:]

    if len(argv) == 0:
        print(disclaimer)
        sys.exit(1)

    database_util = False
    if argv[0] == 'database':
        database_util = True
        argv = argv[1:]

        if len(argv) == 0 or argv[0] == '-h' or argv[0] == '--help':
            print(database_util_disclaimer)
            sys.exit(1)
        if argv[0] == 'reset':
            db_mode = 'reset'
            argv = argv[1:]
        elif argv[0] == 'out':
            db_mode = 'out'
            argv = argv[1:]
        elif argv[0] == 'create':
            db_mode = 'create'
            argv = argv[1:]
        elif argv[0] == 'destroy':
            db_mode = 'destroy'
            argv = argv[1:]
        else:
            print(database_util_disclaimer)
            sys.exit(1)

    update_util = False
    update_pdb = False
    update_rindb = False
    if argv[0] == 'update':
        update_util = True
        argv = argv[1:]

        if len(argv) == 0 or argv[0] == '-h' or argv[0] == '--help':
            print(update_util_disclaimer)
            sys.exit(1)
        if 'pdb' in argv:
            update_pdb = True
        if 'rindb' in argv:
            update_rindb = True
        if not (update_pdb or update_rindb):
            print(update_util_disclaimer)
            sys.exit(1)

    configure_mode = False
    if argv[0] == 'config':
        configure_mode = True
        argv = argv[1:]

        if len(argv) == 0 or argv[0] == '-h' or argv[0] == '--help':
            print(config_util_disclaimer)
            sys.exit(1)

        conf_update_pdb_path = None
        conf_update_iupred_path = None

        if argv[0] == 'set_local_pdb_path':
            if len(argv) == 1:
                print(config_util_disclaimer)
                sys.exit(1)
            conf_update_pdb_path = argv[1]
            if not os.path.exists(conf_update_pdb_path):
                print('Did not found given path')
                sys.exit(1)

        elif argv[0] == 'set_local_iupred_path':
            if len(argv) == 1:
                print(config_util_disclaimer)
                sys.exit(1)
            conf_update_iupred_path = argv[1]
            if not os.path.exists(conf_update_iupred_path):
                print('Did not found given path')
                sys.exit(1)

        else:
            print(config_util_disclaimer)
            sys.exit(1)

        argv = argv[2:]

    output_util = False
    ppi_output = False
    if argv[0] == 'out':
        output_util = True
        if argv[1] == 'PPI':
            ppi_output = True
        else:
            print(disclaimer)
            sys.exit(1)
        argv = argv[2:]

    util_mode = database_util or configure_mode or update_util

    #Custom single line input preparsing
    insert_flag_pos = None
    remove_args = []
    single_line_inputs = []
    for pos,arg in enumerate(argv):
        if insert_flag_pos != None:
            #If the arg after the -i flag gives us a path, then it is not a single line input
            if pos == insert_flag_pos + 1:
                if os.path.exists(arg):
                    break
                else:
                    single_line_inputs.append(arg)
                    remove_args.append(insert_flag_pos)
                    remove_args.append(pos)
            else:
                if arg[0] == '-':
                    break
                else:
                    single_line_inputs.append(arg)
                    remove_args.append(pos)

        if arg == '-i':
            insert_flag_pos = pos
    #ignore --overwrrite by removing it, there will be a bug if someone combines single line input with overwrite
    for pos,arg in enumerate(argv):
        if arg == '--overwrite':
            remove_args.append(pos)
        if arg == 'pdb':
            remove_args.append(pos)
        if arg == 'rindb':
            remove_args.append(pos)

    remove_args.reverse()
    for pos in remove_args:
        del argv[pos]
    
    try:
        opts,args = getopt.getopt(argv,"c:i:n:o:h:lvdp:",['help','profile','skipref','rlimit=','verbosity=','printerrors','chunksize=','norin'])

    except getopt.GetoptError:
        print("Illegal Input\n\n",disclaimer)
        sys.exit(2)

    infile = ''
    config_path = ''
    num_of_cores = None
    outfolder = ''
    lite = False
    verbose_flag = False
    verbosity = None
    profiling = False
    skipref = False
    print_all_errors = False
    chunksize = 500
    '''
    #mmcif mode flag is added
    mmcif_mode = False
    '''
    norin = False
    minus_p_path = None #different modes can use this way a path given with -p

    for opt,arg in opts:
        if opt == '-c':
            config_path = arg
        if opt == '-v':
            verbose_flag = True
        if opt == '-l':
            lite = True
        if opt == '-d':
            lite = False
        if opt == '-i':
            infile = arg
        if opt == '-n':
            num_of_cores = int(arg)
        if opt == '-o':
            outfolder = arg
        if opt == '-p':
            minus_p_path = arg
            if not os.path.exists(minus_p_path):
                try:
                    os.mkdir(minus_p_path)
                except:
                    print('Did not found given path',minus_p_path)
                    sys.exit(1)

        if opt == '--profile':
            profiling = True
        if opt == '-h' or opt == '--help':
            print(disclaimer)
            sys.exit(0)
        if opt == '--rlimit':
            rlimit = arg
        if opt == '--verbosity':
            verbosity = int(arg)
            vebose_flag = True

        if opt == '--skipref':
            skipref = True

        if opt == '--printerrors':
            print_all_errors = True

        if opt == '--chunksize':
            chunksize = int(arg)
        if opt == '--norin':
            norin = True
        '''
        #mmcif option added to call mmcif mode while calling structman
        if opt == '-mmcif':
            mmcif_mode = True
        '''

    if not output_util and not util_mode:
        if infile == '' and len(single_line_inputs) == 0:
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

        elif infile == '' and len(single_line_inputs) != 0:
            lite = True
            infiles = [single_line_inputs]
        else:
            #make infile to global path
            infile = os.path.abspath(infile)
            infiles = [infile]

    if not util_mode:
        if outfolder == '':
            outfolder = '/structman/results/'
            if not os.path.isdir(outfolder):
                if infile != '':
                    outfolder = '%s/Output' % infile.rsplit('/',1)[0]
                    if not os.path.exists(outfolder):
                        os.mkdir(outfolder)
                elif len(single_line_inputs) != 0:
                    outfolder = os.getcwd()
                else:
                    outfolder = None


    if config_path == '':
        if os.path.exists('%s/config.txt' % infile.rsplit('/',1)[0]):
            config_path = '%s/config.txt' % infile.rsplit('/',1)[0]
        elif os.path.exists('%s/config.txt' % os.getcwd()):
            config_path = '%s/config.txt' % os.getcwd()
        elif os.path.exists('%s/config.txt' % main_file_path.rsplit('/',1)[0]):
            config_path = '%s/config.txt' % main_file_path.rsplit('/',1)[0]
        else:
            print("No config file found, please use -c [Path to config]")
            sys.exit(2)


    config_path = os.path.abspath(config_path)

    config = Config(config_path,num_of_cores = num_of_cores,
                    output_path = outfolder,
                    util_mode = util_mode,output_util = output_util ,external_call = False,profiling = profiling,verbosity = verbosity,
                    print_all_errors = print_all_errors,chunksize = chunksize)

    if verbose_flag:
        config.verbose = True
        if verbosity == 0:
            print('-v cannot combined with --verbosity 0; setting verbosity to 1')
            config.verbosity = 1

    if config.verbosity >= 1:
        print(("Using following config file: %s" % config_path))
        config.verbose = True

    config.calculate_interaction_profiles = not norin

    config.skipref = skipref

    if len(single_line_inputs) == 0:
        infile = os.path.abspath(infile)

    if database_util:
        if db_mode == 'reset':
            db = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,config.db_name)
            cursor = db.cursor()
            database.reset(cursor)
            db.close()
        elif db_mode == 'out':
            db = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,config.db_name)
            cursor = db.cursor()
            print(infile)
            session_id = database.getSessionId(infile,db,cursor)
            [trunk,infilename] = infile.rsplit('/',1)
            session_name = infilename.rsplit('.',1)[0]
            outfolder = "%s/Output" % trunk
            outfile = '%s/%s_ligandlist.tsv' % (outfolder,session_name)

            database.getLigandList(db,cursor,session_id,outfile)
            db.close()
        elif db_mode == 'create':
            repairDB.load(config)
        elif db_mode == 'destroy':
            repairDB.destroy(config)


    elif output_util:
        if ppi_output:
            db = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,config.db_name)
            cursor = db.cursor()
            if config.verbosity >= 1:
                print(infile)
            session_id = database.getSessionId(infile,db,cursor)
            [trunk,infilename] = infile.rsplit('/',1)
            session_name = infilename.rsplit('.',1)[0]
            outfolder = "%s/Output" % trunk
            outfile = '%s/%s_ppi_network.tsv' % (outfolder,session_name)
            output.create_ppi_network(session_id,config,outfile)

            db.close()

    elif update_util:
        import update
        if minus_p_path != None:
            config.pdb_path = minus_p_path
            config.config_parser_obj.set('user','pdb_path',minus_p_path)
            f = open(config_path,'w')
            config.config_parser_obj.write(f)
            f.close()
        update.main(config,skipUpdatePDB = not update_pdb,skip_rindb = not update_rindb)

    elif configure_mode:
        if conf_update_pdb_path != None:
            config.config_parser_obj.set('user','pdb_path',conf_update_pdb_path)
        if conf_update_iupred_path != None:
            config.config_parser_obj.set('user','iupred_path',conf_update_iupred_path)
        f = open(config_path,'w')
        config.config_parser_obj.write(f)
        f.close()

    else:

        config.lite = lite

        main(infiles,outfolder,main_file_path,config)
