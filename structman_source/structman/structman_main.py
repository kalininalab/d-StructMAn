#!/usr/bin/python3
# This is the wrapper script, including and automatizing all functions StructMAn has to offer
import configparser
import getopt
import multiprocessing
import os
import resource
import subprocess
import sys
import time
import traceback

import pymysql as MySQLdb
from psutil import virtual_memory

# if running as script, add local structman package to path
if __name__ == "__main__":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__)))))
import structman
from structman.lib import repairDB, serializedPipeline, searchLargeStructures
from structman.lib.output import output
from structman.lib.database import database
from structman.scripts import update
from structman.base_utils.base_utils import resolve_path

# Check if autopep8 is installed for developers. Disable warning for end users by utilizing conda environment variable `STRUCTMAN_WARN`
try:
    if os.environ.get('STRUCTMAN_WARN', 'yes') != 'no':
        import autopep8
except ImportError:
    print("Autopep8 is not installed. For developers, please set up autopep8 git hooks according to the internal wiki")
    user_input = input("If you are not a developer, please input `x` to disable this message and re-activate your conda environment.\n>>> ")
    if user_input.lower().startswith('x'):
        env_path = os.path.join(os.environ['CONDA_PREFIX'], 'etc/conda/')
        os.makedirs(os.path.join(env_path, 'activate.d'), exist_ok=True)
        os.makedirs(os.path.join(env_path, 'deactivate.d'), exist_ok=True)
        with open(os.path.join(env_path, 'activate.d', 'env_vars.sh'), 'a') as file:
            file.write('export STRUCTMAN_WARN=no\n')
        with open(os.path.join(env_path, 'deactivate.d', 'env_vars.sh'), 'a') as file:
            file.write('unset STRUCTMAN_WARN\n')


class Config:
    def __init__(self, config_path, num_of_cores=1, output_path='', basic_util_mode=False,
                 util_mode=False, output_util=False, external_call=True, profiling=False, verbosity=None,
                 print_all_errors=False, print_all_warns=False, restartlog=False):
        self.prog_start_time = time.time()
        # read config file, auto add section header so old config files work
        self.config_parser_obj = configparser.ConfigParser()
        try:
            self.config_parser_obj.read(config_path)
        except configparser.MissingSectionHeaderError:
            with open(config_path, 'r') as f:
                self.config_parser_obj.read_string('[user]\n%s' % f.read())
        cfg = self.config_parser_obj['user']

        # fix typo `adress` -> `address`, but want old configs to still work, so also check for unfixed name
        old_address_name = cfg.get('db_adress', fallback='')
        self.db_address = cfg.get('db_address', fallback=old_address_name)  # fixed address name takes precedence
        self.db_user_name = cfg.get('db_user_name', fallback='')
        self.db_password = cfg.get('db_password', fallback='')
        self.db_name = cfg.get('db_name', fallback='')
        self.mapping_db = cfg.get('mapping_db', fallback=None)

        self.ignore_local_mapping_db = False

        self.user_mail = cfg.get('mail', fallback='')

        self.outfolder = resolve_path(output_path)
        self.num_of_cores = num_of_cores
        self.indels_given_by_input = True
        self.only_wt = False
        self.rare_residues = set()

        self.neighborhood_calculation = cfg.getboolean('neighborhood_calculation', fallback=True)
        self.error_annotations_into_db = cfg.getboolean('error_annotations_into_db', fallback=True)
        self.anno_session_mapping = cfg.getboolean('anno_session_mapping', fallback=True)
        self.calculate_interaction_profiles = cfg.getboolean('calculate_interaction_profiles', fallback=True)
        self.model_indel_structures = cfg.getboolean('model_indel_structures', fallback=False)

        self.verbose = cfg.getboolean('verbose', fallback=False)
        self.verbosity = cfg.getint('verbosity', fallback=1)

        self.profiling = profiling
        self.skipref = cfg.getboolean('skipref', fallback=False)
        self.test_low_mem_system = cfg.getboolean('test_low_mem_system', fallback=False)

        self.resources = cfg.get('resources', fallback='manu')
        self.container_version = cfg.getboolean('container_version', fallback = False)

        self.proc_n = 48
        self.blast_processes = self.proc_n
        self.alignment_processes = self.proc_n
        self.annotation_processes = self.proc_n
        self.number_of_processes = self.proc_n

        self.dssp = cfg.getboolean('dssp', fallback=True)
        self.pdb_input_asymetric_unit = cfg.getboolean('pdb_input_asymetric_unit', fallback=False)
        self.search_tool = cfg.get('search_tool', fallback='MMseqs2')

        # active options
        self.n_of_chain_thresh = cfg.getint('n_of_chain_thresh', fallback=12)  # Structures with n_of_chain_thresh or more chains get a nested paralellization
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
        self.buried_threshold = cfg.getfloat('buried_threshold', fallback=0.05)
        self.short_distance_threshold = cfg.getfloat('lig_short_dist_thresh', fallback=5.0)
        self.long_distance_threshold = cfg.getfloat('lig_long_dist_thresh', fallback=8.0)
        self.ligand_interest_sphere = cfg.getfloat('ligand_interest_sphere', fallback=25.0)

        # Paths to Helper Programs
        self.blast_path = cfg.get('blast_path', fallback='')
        self.mmseqs2_path = cfg.get('mmseqs2_path', fallback='')
        self.output_path = cfg.get('output_path', fallback='')
        self.pdb_path = cfg.get('pdb_path', fallback='')

        self.annovar_path = cfg.get('annovar_path', fallback='')
        self.dssp_path = cfg.get('dssp_path', fallback = '')
        self.rin_db_path = cfg.get('rin_db_path', fallback='')
        self.iupred_path = cfg.get('iupred_path', fallback='')

        self.errorlog_path = cfg.get('errorlog_path', fallback=None)

        self.go = cfg.getboolean('do_goterm', fallback=False)
        self.anno = cfg.getboolean('do_anno', fallback=False)
        self.classification = cfg.getboolean('do_classification', fallback=True)
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
        self.intra_milieu_threshold = cfg.getfloat('intra_milieu_threshold', fallback=15.0)
        self.ligand_filter = cfg.get('ligand_filter', fallback=None)
        self.proteome = cfg.getboolean('proteome', fallback=False)
        self.intertable_conf = cfg.getboolean('intertable', fallback=False)
        self.overwrite_incorrect_wt_aa = cfg.getboolean('overwrite_incorrect_wt_aa', fallback=False)

        self.pdb_sync_script = cfg.get('pdb_sync_script', fallback=structman.settings.PDB_SYNC_SCRIPT)

        self.database_source_path = cfg.get('database_source_path', fallback=structman.settings.STRUCTMAN_DB_SQL)
        self.mapping_db_source_path = cfg.get('mapping_db_source_path', fallback=structman.settings.UNIPROT_DB_SQL)

        self.base_path = cfg.get('base_path', fallback='')

        self.rinerator_base_path = structman.settings.RINERATOR_DIR
        self.rinerator_path = structman.settings.RINERATOR_PATH
        os.environ["REDUCE_HET_DICT"] = structman.settings.REDUCE_HET_DICT

        self.mmseqs_tmp_folder = cfg.get('mmseqs_tmp_folder')
        self.tmp_folder = self.mmseqs_tmp_folder

        self.fasta_input = cfg.getboolean('fasta_input', fallback=False)
        self.lite = cfg.getboolean('lite', fallback=False)

        if verbosity is not None:
            self.verbosity = verbosity

        # Checking whether the given paths in config exist or not if it is given and not exist, system gives error message and exits

        if self.search_tool == 'MMseqs2':
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
        elif self.search_tool == 'Blast':
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
            if not os.access(self.mmseqs_tmp_folder, os.R_OK):
                print('Need writing rights in the MMSEQS temp folder, please check the path')
                sys.exit()

        if self.base_path is not None and not basic_util_mode:
            if self.verbosity >= 1:
                print('Using structman_data from :', self.base_path)
            self.mmseqs2_db_path = '%s/base/blast_db/pdbba_search_db_mmseqs2' % self.base_path
            self.smiles_path = '%s/base/ligand_bases/Components-smiles-stereo-oe.smi' % self.base_path
            self.inchi_path = '%s/base/ligand_bases/inchi_base.tsv' % self.base_path
            self.human_id_mapping_path = '%s/base/id_mapping' % self.base_path

        if self.resources == 'auto' and self.num_of_cores is None:
            self.proc_n = multiprocessing.cpu_count() - 1
            if self.proc_n <= 0:
                self.proc_n = 1

            self.blast_processes = self.proc_n
            self.alignment_processes = self.proc_n
            self.annotation_processes = self.proc_n
            self.number_of_processes = self.proc_n

        if self.num_of_cores is not None:
            self.proc_n = self.num_of_cores

            self.blast_processes = self.proc_n
            self.alignment_processes = self.proc_n
            self.annotation_processes = self.proc_n
            self.number_of_processes = self.proc_n

        if self.proc_n > multiprocessing.cpu_count():
            if self.verbosity >= 1:
                print('More processes annotated (', self.proc_n, ') than cores registered in system (', multiprocessing.cpu_count(), ').')
            self.proc_n = multiprocessing.cpu_count()
            self.blast_processes = self.proc_n
            self.alignment_processes = self.proc_n
            self.annotation_processes = self.proc_n
            self.number_of_processes = self.proc_n

        mem = virtual_memory()
        self.gigs_of_ram = mem.total / 1024 / 1024 / 1024
        self.low_mem_system = self.gigs_of_ram < 40 or self.test_low_mem_system  # Less than 20Gb is a low memory system
        if self.test_low_mem_system:
            self.gigs_of_ram = 8
        if self.low_mem_system:
            self.chunksize = int(max([((self.gigs_of_ram * 90) // self.proc_n) - 120, 60 // self.proc_n, 1]))
        else:
            self.chunksize = int(max([((self.gigs_of_ram * 150) // self.proc_n) - 60, 120 // self.proc_n, 1]))

        if not util_mode:
            if not external_call and not os.path.exists(self.outfolder):
                os.makedirs(self.outfolder)
            if self.verbosity >= 1:
                print('Using %s core(s)' % str(self.proc_n))
            if (not external_call) and (not print_all_warns):  # no need for errorlogs, when the config is generated not from the main script
                self.errorlog_path = os.path.join(self.outfolder, 'errorlogs', 'errorlog.txt')
                errorlog_dir = os.path.dirname(self.errorlog_path)
                if not os.path.exists(errorlog_dir):
                    os.makedirs(errorlog_dir)
                if restartlog and os.path.isfile(self.errorlog_path):
                    os.remove(self.errorlog_path)

        self.errorlog = Errorlog(path=self.errorlog_path, print_all_errors=print_all_errors, print_all_warns=print_all_warns)

        # Determine maximal package size from database
        try:
            db, cursor = self.getDB(server_connection=True)
        except:
            if self.verbosity >= 2:
                print('Database connection failed in config initialization')
            db = None

        if db is not None:
            cursor.execute("SHOW VARIABLES WHERE variable_name = 'max_allowed_packet'")
            self.max_package_size = int(cursor.fetchone()[1]) * 100 // 99
            db.close()
            self.main_db_is_set = True
        else:
            self.max_package_size = None
            self.main_db_is_set = False

        # Check for mapping DB instance
        self.check_mapping_db()

    def check_mapping_db(self):
        # Check for mapping DB instance
        try:
            db, cursor = self.getDB(mapping_db = True)
            db.close()
            self.mapping_db_is_set = True
        except:
            self.mapping_db_is_set = False

    def getDB(self, server_connection=False, mapping_db=False):
        if server_connection:
            db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password)
            cursor = db.cursor()
        elif mapping_db:
            if self.ignore_local_mapping_db:
                return None, None

            if self.mapping_db is None:
                return None, None
            db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password, database=self.mapping_db)
            cursor = db.cursor()
        else:
            db = MySQLdb.connect(host=self.db_address, user=self.db_user_name, password=self.db_password, database=self.db_name)
            cursor = db.cursor()
        return db, cursor


class Errorlog:
    def __init__(self, path=None, print_all_errors=False, print_all_warns=False):
        self.path = path
        self.error_counter = 0
        self.warning_counter = 0
        self.print_all_errors = print_all_errors
        self.print_all_warns = print_all_warns

    def start(self, nfname, session):
        if self.path is None:
            return
        date = time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())
        errortext = "###############################################################################\n%s:\n%s - Session: %s\n" % (date, nfname, str(session))
        f = open(self.path, 'a')
        f.write(errortext)
        f.close()

    def add_error(self, error_text):
        self.error_counter += 1
        g = ''.join(traceback.format_list(traceback.extract_stack()))
        error_text = 'Error %s:\n%s\n%s\n' % (str(self.error_counter), error_text, g)
        if self.path is None or self.print_all_errors or self.print_all_warns:
            print(error_text)
            return
        f = open(self.path, 'a')
        f.write(error_text)
        f.close()

    def stop(self):
        if self.path is None:
            return
        if self.error_counter == 0:
            errortext = "Finished without any error\n###############################################################################\n"
            f = open(self.path, 'a')
            f.write(errortext)
            f.close()
        else:
            errortext = "###############################################################################\n"
            f = open(self.path, 'a')
            f.write(errortext)
            f.close()
            print("\n\nAt least one error occured, please check the errorlog.\n\n")

    def add_warning(self, warn_text, lock=None):
        self.warning_counter += 1
        g = ''.join(traceback.format_list(traceback.extract_stack()))
        warn_text = 'Warning %s:\n%s\n%s\n' % (str(self.warning_counter), warn_text, g)
        if self.path is None or self.print_all_warns:
            print(warn_text)
            return
        if lock is not None:
            with lock:
                f = open(self.path, 'a')
                f.write(warn_text)
                f.close()
            return
        f = open(self.path, 'a')
        f.write(warn_text)
        f.close()

def limit_memory(soft_limit):
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (soft_limit, hard))

def main(infiles, out_folder, config, intertable=False):
    db_name = config.db_name
    db_address = config.db_address
    db_password = config.db_password
    db_user_name = config.db_user_name
    # resource.setrlimit(resource.RLIMIT_AS,(6442450944,8589934592))
    # resource.setrlimit(resource.RLIMIT_AS,(442450944,589934592))
    for infile in infiles:
        if isinstance(infile, str):
            if config.verbosity >= 1:
                print('Processing file: ', infile)
            if not config.lite:
                # check if the infile got already processed
                session_id = database.getSessionId(infile, config)

            # create the output folder
            trunk, infilename = os.path.split(infile)

            if infilename[-6:] == '.fasta':
                config.fasta_input = True

            if out_folder is None:
                out_folder = resolve_path(os.path.join(trunk, "Output"))

        if not os.path.exists(out_folder):
            os.makedirs(out_folder)

        if not config.lite:
            # run the main pipeline
            if session_id is None:
                session_id = serializedPipeline.main(infile, config)

            # run the output scripts
            session_name = infilename.rsplit('.', 1)[0]
            months = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun', 7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
            time_struct = time.gmtime()
            year = str(time_struct[0])
            month = months[time_struct[1]]
            day = str(time_struct[2])
            date = "%s%s%s" % (day, month, year)
            outpath = os.path.join(out_folder, session_name)
            if not os.path.exists(outpath):
                os.makedirs(outpath)

            output.main(session_id, outpath, config, intertable=intertable,)
        else:
            session_id = serializedPipeline.main(infile, config)


# needed for entry point in setup.py
def structman_cli():
    main_file_path = resolve_path(__file__)

    disclaimer = ''.join([
        'Usage: structman.py <-i input_file>\n',
        'more functionalities can be used giving a second key word:\n',
        'structman.py database    gives you more info about the database utility functions\n',
        'structman.py config      gives you more info about the config utility functions\n',
        'structman.py update      gives you more info about the different update options\n\n\n',
        '| Optional parameter  | Default value           | Description                               |\n',
        '|---------------------|-------------------------|-------------------------------------------|\n',
        '| <-n threads>        | All available -1        | Number of cores to be used                |\n',
        '|---------------------|-------------------------|-------------------------------------------|\n',
        '| <-o output_folder>  | standard installation:  | Path to the output folder                 |\n',
        '|                     | Output/infilename/      |                                           |\n',
        '|                     | container installation: |                                           |\n',
        '|                     | structman/results/      |                                           |\n',
        '|---------------------|-------------------------|-------------------------------------------|\n',
        '| <-c config_file>    | standard installation:  | Path to the configuration file            |\n',
        '|                     | config.txt              |                                           |\n',
        '|                     | container installation: | create custom config.txt inside mounted   |\n',
        '|                     | not mounted             | folder, or use <structman config>         |\n',
        '|---------------------|-------------------------|-------------------------------------------|\n',
        '| <-l>                | False                   | lite mode (deactivates the database)      |\n',
        '|---------------------|-------------------------|-------------------------------------------|\n',
        '| <--verbosity> [0-5] | 1                       | verbosity of command line messages        |\n',
        '|                     |                         | 0: silence; 1: minimal progress reports   |\n',
        '|                     |                         | 2: more detailed progress reports         |\n',
        '|                     |                         | 3: details for debugging, 4: more details |\n',
        '|                     |                         | 5: too much details                       |\n',
        '|                     |                         | pipe output to a file for verbosity >= 3  |\n',
        '|---------------------|-------------------------|-------------------------------------------|\n',
        '| <--printerrors>     | False                   | prints error messages in the console      |\n',
        '|                     |                         | instead of logging them. Warning messages |\n',
        '|                     |                         | still go into the log                     |\n',
        '|---------------------|-------------------------|-------------------------------------------|\n',
        '| <--printwarnings>   | False                   | prints error and warning messages in      |\n',
        '|                     |                         | the console instead of logging them       |\n',
        '|---------------------|-------------------------|-------------------------------------------|\n',
        '| <--restartlog>      | False                   | wipes the log before starting the session |\n',
        '|---------------------|-------------------------|-------------------------------------------|\n',
        '| <--norin>           | False                   | disable all RIN-based calculation         |\n',
        '|---------------------|-------------------------|-------------------------------------------|\n'
    ])

    argv = sys.argv[1:]

    if len(argv) == 0:
        print(disclaimer)
        sys.exit(1)

    database_util_disclaimer = ''.join([
        'Usage: structman.py database [command] <-c config_file>\n\n',
        '#### Commands: ####\n\n',
        'reset :\n  deletes all content of the database\n\n',
        'clear :\n  deletes all content of the database, but keeps the results of all stored structures\n\n',
        'export <-p path_to_a_folder> :\n  exports the database to a .sql.gz format file\n\n',
        'destroy :\n  completely removes the database (database name is taken from config file)\n\n',
        'create  :\n  creates an instance of the database (database name is taken from config file)'
    ])

    database_util = False
    if argv[0] == 'database':
        database_util = True
        argv = argv[1:]

        possible_key_words = set(['reset', 'out', 'create', 'destroy', 'clear', 'export', 'reduce', 'remove_sessions'])

        if len(argv) == 0 or argv[0] == '-h' or argv[0] == '--help':
            print(database_util_disclaimer)
            sys.exit(1)
        if argv[0] in possible_key_words:
            db_mode = argv[0]
            argv = argv[1:]
        else:
            print(database_util_disclaimer)
            sys.exit(1)

    update_util_disclaimer = ''.join([
        'Usage: structman.py update [commands] <-c config_file> <-p path_to_local_pdb>\n',
        'The update functionalities manage local database, that significantly increases the performance of StructMAn.\n',
        'They need a lot of disk space though. Thus the usage of local databases is only recommended when one wants to process larger inputs with StructMAn.\n'
        'Only needs the path to the local instance of the PDB (-p), if it is not already specified in the config file.\n\n'
        '#### Commands: ####\n\n',

        'pdb                       uses rsync to update all files of the local PDB that could potentially used by StructMAn.\n',
        '                          If the given path to the local PDB is empty, a large portion of the PDB will be downloaded,\n',
        '                          be sure that there is enought disk space available.\n\n',

        'rindb                     calculates and stores the RIN for each PDB structure in the local PDB.\n',
        '                          When called for the first time, this can take multiple hours to complete.\n',
        '                          Also requires a large amount of disk space.\n\n',

        'rindb_from_scratch        same as rindb, but overwrites all RINs currently in the RINdb.\n\n',

        'mapping_db                downloads protein identifier mapping information and sequences from UNIPROT.\n',
        '                          Builds a SQL database and removes the downloaded files again.\n',
        '                          Requires a configured SQL database server connection.\n',
        '                          Requires around 100Gb of temporal disk space and 35Gb of space on the SQL server.\n\n',

        'mapping_db_from_scratch   same as mapping_db but removes an existing instance of the mapping DB.\n'
    ])

    update_util = False
    update_pdb = False
    update_rindb = False
    update_rindb_from_scratch = False
    update_mapping_db = False
    update_mapping_db_from_scratch = False
    update_mapping_db_keep_raw_files = False

    if len(argv) > 0:
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
            if 'rindb_from_scratch' in argv:
                update_rindb = True
                update_rindb_from_scratch = True
            if 'mapping_db' in argv:
                update_mapping_db = True
            if 'mapping_db_from_scratch' in argv:
                update_mapping_db = True
                update_mapping_db_from_scratch = True
            if 'mapping_db_keep_raw_files' in argv:
                update_mapping_db = True
                update_mapping_db_keep_raw_files = True
            if not (update_pdb or update_rindb or update_mapping_db):
                print(update_util_disclaimer)
                sys.exit(1)

    config_util_disclaimer = ''.join([
        'Usage: structman.py config [command] [value] <-c config_file>\n',
        'The config commands enable the expansion of StructMAn by giving it access to additional functionalities, which are too large to be included by default or require an external license.\n\n',
        '#### Commands: ####\n\n',
        'set_local_pdb_path <path_to_local_pdb> :               enables the usage of a local instance of the PDB. If this is not available, StructMAn has to download all structural information from the web.\n\n',
        'set_local_iupred_path <path_to_iupred_executable> :    enables StructMAn to use the predicted disordered regions performed by iupred3.\n\n',
        '<any config variable name> <any value> :               modifies any variable in the config file. Warning: does not check if given input makes sense.',
    ])

    configure_mode = False
    if len(argv) > 0:
        if argv[0] == 'config':
            configure_mode = True
            argv = argv[1:]

            if len(argv) == 0 or argv[0] == '-h' or argv[0] == '--help':
                print(config_util_disclaimer)
                sys.exit(1)

            conf_update_pdb_path = None
            conf_update_iupred_path = None

            config_variable_name = None

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
                if len(argv) == 1:
                    print(config_util_disclaimer)
                    sys.exit(1)
                config_variable_name = argv[0]
                config_variable_value = argv[1]

            argv = argv[2:]

    output_util = False

    if len(argv) > 0:
        if argv[0] == 'out':
            argv = argv[1:]
            possible_key_words = set(['suggest', 'PPI', 'LSS', 'FAT', 'MM', 'RAS'])
            out_util_mode = None
            output_util = True
            if argv[0] in possible_key_words:
                out_util_mode = argv[0]
                argv = argv[1:]
            else:
                print(disclaimer)
                sys.exit(1)

    basic_util_mode = database_util or configure_mode
    util_mode = database_util or configure_mode or update_util or output_util

    # Custom single line input preparsing
    insert_flag_pos = None
    remove_args = []
    single_line_inputs = []
    for pos, arg in enumerate(argv):
        if insert_flag_pos is not None:
            # If the arg after the -i flag gives us a path, then it is not a single line input
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
    # ignore --overwrrite by removing it, there will be a bug if someone combines single line input with overwrite
    for pos, arg in enumerate(argv):
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
        long_paras = [
            'help', 'profile', 'skipref', 'rlimit=', 'verbosity=',
            'printerrors', 'printwarnings', 'chunksize=', 'norin',
            'dbname=', 'restartlog', 'only_snvs', 'skip_indel_analysis', 'only_wt', 'mem_limit=',
            'model_indel_structures', 'ignore_local_pdb', 'ignore_local_mapping_db'
        ]
        opts, args = getopt.getopt(argv, "c:i:n:o:h:lvdp:", long_paras)

    except getopt.GetoptError:
        print("Illegal Input\n\n", disclaimer)
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
    print_all_warns = False
    chunksize = None
    dbname = None
    restartlog = False
    only_snvs = False
    skip_indel_analysis = False
    only_wt = False
    model_indel_structures = None
    ignore_local_pdb = False
    ignore_local_rindb = False
    ignore_local_mapping_db = False
    '''
    #mmcif mode flag is added
    mmcif_mode = False
    '''
    norin = False
    minus_p_path = None  # different modes can use this way a path given with -p
    mem_limit = None

    for opt, arg in opts:
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
                    os.makedirs(minus_p_path)
                except:
                    print('Did not found given path', minus_p_path)
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

        if opt == '--skipref':
            skipref = True

        if opt == '--printerrors':
            print_all_errors = True

        if opt == '--printwarnings':
            print_all_warns = True

        if opt == '--chunksize':
            chunksize = int(arg)
        if opt == '--norin':
            norin = True
        if opt == '--dbname':
            dbname = arg
        if opt == '--restartlog':
            restartlog = True
        '''
        #mmcif option added to call mmcif mode while calling structman
        if opt == '-mmcif':
            mmcif_mode = True
        '''
        if opt == '--only_snvs':
            only_snvs = True

        if opt == '--skip_indel_analysis':
            skip_indel_analysis = True

        if opt == '--only_wt':
            only_wt = True
        
        if opt == '--mem_limit':
            mem_limit = float(arg)

        if opt == '--model_indel_structures':
            model_indel_structures = True

        if opt == '--ignore_local_pdb':
            ignore_local_pdb = True

        if opt == '--ignore_local_rindb':
            ignore_local_rindb = True

        if opt == '--ignore_local_mapping_db':
            ignore_local_mapping_db = True

    if not output_util and not util_mode:
        if infile == '' and len(single_line_inputs) == 0:
            input_folder = '/structman/input_data/'
            if not os.path.isdir(input_folder):
                print('Did not find the input path\n\n', disclaimer)
                sys.exit(2)
            filelist = os.listdir(input_folder)

            infiles = []
            for infile in filelist:
                if infile.split('.')[-1] == 'smlf' or infile.split('.')[-1] == 'vcf':
                    infiles.append(os.path.join(input_folder, infile))
            if infiles == []:
                print('Did not find any file inside the input folder\n\n', disclaimer)
                sys.exit(2)

        elif infile == '' and len(single_line_inputs) != 0:
            lite = True
            infiles = [single_line_inputs]
        else:
            # make infile to global path
            infile = resolve_path(infile)
            infiles = [infile]

    if not util_mode:
        if outfolder == '':
            outfolder = '/structman/results/'
            if not os.path.isdir(outfolder):
                if infile != '':
                    outfolder = resolve_path(os.path.join(os.path.dirname(infile), 'Output'))
                    if not os.path.exists(outfolder):
                        os.makedirs(outfolder)
                elif len(single_line_inputs) != 0:
                    outfolder = os.getcwd()
                else:
                    outfolder = None

    if config_path == '':
        # search highest to lowest priority
        search_paths = [
            os.path.dirname(infile),  # same directory as input file
            os.getcwd(),  # directory where structman command was triggered
            os.path.dirname(os.path.dirname(main_file_path)),  # top level of structman repo
            os.path.dirname(main_file_path),  # inside structman package (same level as structman_main.py)
            os.path.expanduser('~'),  # user's home folder
            '/structman/resources/', # container source folder
        ]
        search_names = ['config.txt', 'structman_config.txt']  # possible config filenames
        # search all possible paths with all possible filenames
        search_files = [os.path.join(path, name) for path in search_paths for name in search_names]
        for cfg_file in search_files:
            if os.path.exists(cfg_file):
                config_path = resolve_path(cfg_file)
                break

        if not config_path:
            print("No config file found, please use -c [Path to config]")
            sys.exit(2)

    config = Config(config_path, num_of_cores=num_of_cores,
                    output_path=outfolder, util_mode=util_mode,
                    basic_util_mode=basic_util_mode, output_util=output_util, external_call=False, profiling=profiling, verbosity=verbosity,
                    print_all_errors=print_all_errors, print_all_warns=print_all_warns, restartlog=restartlog)

    config.only_snvs = only_snvs
    config.skip_indel_analysis = skip_indel_analysis
    config.only_wt = only_wt
    if model_indel_structures is not None:
        config.model_indel_structures = model_indel_structures

    if ignore_local_pdb:
        config.pdb_path = ''

    if ignore_local_rindb:
        config.rin_db_path = ''

    config.ignore_local_mapping_db = ignore_local_mapping_db
    if ignore_local_mapping_db:
        config.mapping_db_is_set = False

    if mem_limit is not None:
        mem = virtual_memory()
        lim = mem.total * mem_limit
        limit_memory(lim)
        print('Limit the memory consumption to %s%% (%s GB of RAM)' % (str(100 * mem_limit), str(lim / 1024. / 1024. / 1024.)))

    if chunksize is not None:
        config.chunksize = chunksize

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
        infile = resolve_path(infile)

    if database_util:
        if db_mode == 'reset':
            repairDB.empty(config)
        elif db_mode == 'clear':
            repairDB.clear(config)
        elif db_mode == 'out':
            db = MySQLdb.connect(config.db_address, config.db_user_name, config.db_password, config.db_name)
            cursor = db.cursor()
            # print(infile)
            session_id = database.getSessionId(infile, db, cursor)
            trunk, infilename = os.path.split(infile)
            session_name = infilename.rsplit('.', 1)[0]
            outfolder = resolve_path(os.path.join(trunk, 'Output'))
            outfile = os.path.join(outfolder, '%s_ligandlist.tsv' % session_name)

            database.getLigandList(db, cursor, session_id, outfile)
            db.close()
        elif db_mode == 'create':
            if infile is not None and os.path.isfile(infile):
                config.database_source_path = infile

            if dbname is not None:
                if dbname[:len(config.db_user_name)] != config.db_user_name:
                    dbname = '%s_%s' % (config.db_user_name, dbname)
                config.config_parser_obj.set('user', 'db_name', dbname)
                config.db_name = dbname
                f = open(config_path, 'w')
                config.config_parser_obj.write(f)
                f.close()

            repairDB.load(config)
        elif db_mode == 'destroy':
            repairDB.destroy(config)
        elif db_mode == 'export':
            if minus_p_path is None:
                print('database export needs a target_path given with -p [path to the target folder]')
                sys.exit(1)

            repairDB.export(config, minus_p_path)
        elif db_mode == 'reduce':
            if infile is None:
                sys.exit(1)
            repairDB.reduceToStructures(config, infile)
        elif db_mode == 'remove_sessions':
            repairDB.remove_sessions(config)

    elif output_util:
        if config.verbosity >= 1:
            print(infile)
        session_id = database.getSessionId(infile, config)
        trunk, infilename = os.path.split(infile)
        session_name = infilename.rsplit('.', 1)[0]
        outfolder = os.path.join(trunk, 'Output')

        if out_util_mode == 'suggest':
            outfile = os.path.join(outfolder, f'{session_name}_template_suggestions.tsv')
            anno_table = os.path.join(outfolder, '%s_full_annotation_table.tsv' % session_name)
            if not os.path.exists(anno_table):
                output.create_annotation_table(session_id, config, anno_table)
            output.suggest(session_id, config, anno_table, outfile, ranking='consensus')
        elif out_util_mode == 'PPI':
            outfile = os.path.join(outfolder, '%s_ppi_network.tsv' % session_name)
            output.create_ppi_network(session_id, config, outfile)
        elif out_util_mode == 'LSS':
            searchLargeStructures.search(config, chunksize, infile=infile)
        elif out_util_mode == 'FAT':
            outfile = os.path.join(outfolder, '%s_full_annotation_table.tsv' % session_name)
            output.create_annotation_table(session_id, config, outfile)
        elif out_util_mode == 'MM':
            outfolder = os.join(outfolder, session_name)
            output.mass_model(session_id, config, outfolder)
        elif out_util_mode == 'RAS':
            outfile = os.path.join(outfolder, '%s_all_annotated_structures.smlf' % session_name)
            output.retrieve_annotated_structures(session_id, config, outfile)

    elif update_util:
        if minus_p_path is not None:
            config.pdb_path = minus_p_path
            config.config_parser_obj.set('user', 'pdb_path', minus_p_path)
            f = open(config_path, 'w')
            config.config_parser_obj.write(f)
            f.close()
        elif config.container_version:
            if update_pdb:
                if not os.path.exists('/structman/resources/pdb'):
                    os.mkdir('/structman/resources/pdb')
                config.pdb_path = '/structman/resources/pdb'
                config.config_parser_obj.set('user', 'pdb_path', '/structman/resources/pdb')

            if update_rindb:
                if not os.path.exists('/structman/resources/rindb'):
                    os.mkdir('/structman/resources/rindb')
                config.rin_db_path = '/structman/resources/rindb'
                config.config_parser_obj.set('user', 'rin_db_path', '/structman/resources/rindb')

            if update_mapping_db:
                config.mapping_db = 'struct_man_db_mapping'
                config.config_parser_obj.set('user', 'mapping_db', 'struct_man_db_mapping')
                config.check_mapping_db()

            f = open(config_path, 'w')
            config.config_parser_obj.write(f)
            f.close()
        update.main(config, skipUpdatePDB=not update_pdb,
                            skip_rindb=not update_rindb,
                            rin_fromScratch=update_rindb_from_scratch,
                            update_mapping_db = update_mapping_db,
                            mapping_db_from_scratch = update_mapping_db_from_scratch,
                            update_mapping_db_keep_raw_files = update_mapping_db_keep_raw_files
                    )

    elif configure_mode:
        if conf_update_pdb_path is not None:
            config.config_parser_obj.set('user', 'pdb_path', conf_update_pdb_path)
        if conf_update_iupred_path is not None:
            config.config_parser_obj.set('user', 'iupred_path', conf_update_iupred_path)

        if config_variable_name is not None:
            config.config_parser_obj.set('user', config_variable_name, config_variable_value)

        f = open(config_path, 'w')
        config.config_parser_obj.write(f)
        f.close()

    else:
        config.lite = lite

        main(infiles, outfolder, config)


if __name__ == "__main__":
    structman_cli()
