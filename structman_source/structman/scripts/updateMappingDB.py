import os
import subprocess
import sys
import traceback
import gzip

from structman.lib.database import database
from structman.utils import pack

def check_instance(config, fromScratch = False):
    fresh_instance = False
    if config.mapping_db_is_set and fromScratch:
        db, cursor = config.getDB(mapping_db = True)
        sql = 'DROP DATABASE %s' % config.mapping_db
        cursor.execute(sql)
        db.close()
        fresh_instance = True
    if not config.mapping_db_is_set or fresh_instance:
        db, cursor = config.getDB(server_connection=True)

        sql = 'CREATE DATABASE %s' % config.mapping_db

        try:
            cursor.execute(sql)
            db.close()
        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            print('\n'.join([str(e), str(f), str(g)]))
            db.close()
            sys.exit()

        if config.mapping_db_source_path[-3:] == '.gz':
            f = gzip.open(config.mapping_db_source_path, 'rb')
            binary = True
        else:
            f = open(config.mapping_db_source_path, 'r')
            binary = False
        lines = f.readlines()
        f.close()

        new_lines = []
        for line in lines:
            if binary:
                line = line.decode('ascii')
            if line[:4] == 'USE ':
                new_lines.append('USE `%s`;\n' % config.mapping_db)
            else:
                new_lines.append(line)

        db_file = 'tmp_mapping_database_file.sql'
        f = open(db_file, 'w')
        f.write(''.join(new_lines))
        f.close()

        cmds = ' '.join(['mysql', '-u', config.db_user_name, '-h', config.db_address, '--password=%s' % config.db_password, '<', db_file])

        p = subprocess.Popen(cmds, shell=True)
        p.wait()

        os.remove(db_file)
        
    return fresh_instance

def retrieve_raw_data(config, fromScratch = False):
    mapping_file_path = retrieve_data_from_uniprot(config, 'idmapping', 'idmapping.dat.gz', fromScratch = fromScratch)
    sequence_file_paths = []
    for seq_file_name in ['uniprot_sprot_varsplic.fasta.gz','uniprot_sprot.fasta.gz','uniprot_trembl.fasta.gz']:
        sequence_file_paths.append(retrieve_data_from_uniprot(config, 'complete', seq_file_name, fromScratch = fromScratch))
    return mapping_file_path, sequence_file_paths

def retrieve_data_from_uniprot(config, uniprot_sub_folder, uniprot_file_name, fromScratch = False):
    file_path = '%s/%s' % (config.tmp_folder, uniprot_file_name)
    if os.path.isfile(file_path) and fromScratch:
        os.remove(file_path)
    if not os.path.isfile(file_path):
        cmds = ' '.join(['wget', 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/%s/%s' % (uniprot_sub_folder, uniprot_file_name)])
        p = subprocess.Popen(cmds, shell=True, cwd = config.tmp_folder)
        p.wait()
    return file_path

def put_seqs_to_database(seq_map, config):
    values = []
    for u_ac in seq_map:
        values.append((u_ac, pack(seq_map[u_ac])))

    database.update(config, 'UNIPROT', ['Uniprot_Ac','Sequence'], values, mapping_db = True)

def main(config, fromScratch = False):
    #Step 1: Check if mapping SQL DB instance is there, create if not. Recreate for fromScratch mode
    fresh_instance = check_instance(config, fromScratch = fromScratch)

    print('Mapping DB instance checked, fresh_instance:', fresh_instance)

    #Step 2: Check for raw files and download if necessary
    mapping_file_path, seq_file_paths = retrieve_raw_data(config, fromScratch = fromScratch)

    print('\nDownloading all raw data files done.\n')

    #Step 3: Update the database

    ac_id_values = []
    ac_ref_values = []
    ac_ref_nt_values = []

    max_values_at_a_time = 1000000 * config.gigs_of_ram

    with gzip.open(mapping_file_path, 'rb') as f:
        for line in f:
            words = line.decode('ascii').split()
            if len(words) == 0:
                continue
            u_ac = words[0]
            id_name = words[1]
            id_value = words[2]

            if u_ac[-2:] == '-1':
                u_ac = u_ac[:-2]

            if id_name == 'UniProtKB-ID':
                ac_id_values.append((u_ac, id_value))

                if len(ac_id_values) == max_values_at_a_time:
                    database.update(config, 'UNIPROT', ['Uniprot_Ac', 'Uniprot_Id'], ac_id_values, mapping_db = True)
                    ac_id_values = []

            elif id_name == 'RefSeq':
                ac_ref_values.append((u_ac, id_value))

                if len(ac_ref_values) == max_values_at_a_time:
                    database.update(config, 'UNIPROT', ['Uniprot_Ac', 'RefSeq'], ac_ref_values, mapping_db = True)
                    ac_ref_values = []

            elif id_name == 'RefSeq_NT':
                ac_ref_nt_values.append((u_ac, id_value))

                if len(ac_ref_nt_values) == max_values_at_a_time:
                    database.update(config, 'UNIPROT', ['Uniprot_Ac', 'RefSeq_NT'], ac_ref_nt_values, mapping_db = True)
                    ac_ref_nt_values = []

    print('\nParsing of mapping data file done.\n')

    if len(ac_id_values) > 0:
        database.update(config, 'UNIPROT', ['Uniprot_Ac', 'Uniprot_Id'], ac_id_values, mapping_db = True)
    print('\nDatabase update of Uniprot IDs done.\n')

    if len(ac_ref_values) > 0:
        database.update(config, 'UNIPROT', ['Uniprot_Ac', 'RefSeq'], ac_ref_values, mapping_db = True)
    print('\nDatabase update of RefSeq IDs done.\n')

    if len(ac_ref_nt_values) > 0:
        database.update(config, 'UNIPROT', ['Uniprot_Ac', 'RefSeq_NT'], ac_ref_nt_values, mapping_db = True)
    print('\nDatabase update of RefSeq NTs done.\n')

    max_seqs_at_a_time = 100000 * config.gigs_of_ram

    for seq_file in seq_file_paths:
        with gzip.open(seq_file, 'rb') as f:
            seq_map = {}
            for line in f:
                line = line.decode('ascii')
                if len(line) == 0:
                    continue
                line = line[:-1]
                if line[0] == '>':

                    if len(seq_map) == max_seqs_at_a_time:
                        put_seqs_to_database(seq_map, config)
                        seq_map = {}

                    u_ac = line.split('|')[1]
                    seq_map[u_ac] = ''
                else:
                    seq_map[u_ac] += (line)
            if len(seq_map) > 0:
                put_seqs_to_database(seq_map, config)
            
            print('\nDatabase update of sequences done.\n')

    #Step 4: Remove the raw files
    os.remove(mapping_file_path)
    for seq_file in seq_file_paths:
        os.remove(seq_file)

    print('\nRemoving raw data files done.\n')
