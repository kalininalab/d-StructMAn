#!/usr/bin/python3
import gzip
import os
import shutil
import subprocess
import sys

from structman import Config, settings


if __name__ == "__main__":
    target_folder = sys.argv[1]
    config_path = sys.argv[2]

    if not os.path.isfile(config_path):
        print('ERROR: Need path to config file as second argument.')
        sys.exit(1)

    structman_target_folder = f'{target_folder}/StructMAn'
    lib_target_folder = f'{structman_target_folder}/lib'
    rinerator_target_folder = f'{lib_target_folder}/rinerator'
    database_target_file = f'{target_folder}/StructMAn_db/struct_man_db.sql.gz'

    f = open(settings.STRUCTMAN_DB_SQL, 'r')

    lines = f.readlines()
    f.close()
    new_lines = []
    for pos, line in enumerate(lines):
        #line = line.decode('ascii')
        if line[:13] == '-- Datenbank:' or line[:4] == 'USE ':

            new_lines.append(b'--\n')
            new_lines.append(b'-- Database: `struct_man_db_1`\n')
            new_lines.append(b'--\n')
            new_lines.append(b'CREATE DATABASE IF NOT EXISTS `struct_man_db_1` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci;\n')
            new_lines.append(b'USE `struct_man_db_1`;\n')
        else:
            new_lines.append(line.encode())

    f = gzip.open(database_target_file, 'wb')
    f.write(b''.join(new_lines))
    f.close()

    #p = subprocess.Popen(['split', '-b', '45M', 'struct_man_db.sql.gz', 'db_split'], cwd='%s/StructMAn_db/' % target_folder)
    #p.wait()

    #if os.path.isfile(database_target_file):
    #    os.remove(database_target_file)

    shutil.copy(f'{settings.ROOT_DIR}/structman_main.py', structman_target_folder)
    shutil.copytree(settings.SCRIPTS_DIR, structman_target_folder, dirs_exist_ok = True)

    utils_lib_path = f'{settings.ROOT_DIR}/structman/utils'
    utils_target_path = f'{target_folder}/StructMAn/utils'
    for utils_file in os.listdir(utils_lib_path):
        source_path = f'{utils_lib_path}/{utils_file}'
        target_path = f'{utils_target_path}/{utils_file}'
        if os.path.isfile(source_path):
            shutil.copy(source_path, target_path)

    for libfile in os.listdir(settings.LIB_DIR):
        source_path = f'{settings.LIB_DIR}/{libfile}'
        target_path = f'{lib_target_folder}/{libfile}'
        if os.path.isfile(source_path):
            shutil.copy(source_path, target_path)

    for rinfile in os.listdir(settings.RINERATOR_DIR):
        source_path = f'{settings.RINERATOR_DIR}/{rinfile}'
        target_path = f'{rinerator_target_folder}/{rinfile}'
        if os.path.isfile(source_path):
            shutil.copy(source_path, target_path)

    lib_subfolders = ['database', 'output', 'sdsc']
    for lib_subfolder in lib_subfolders:
        subfolder_path = f'{settings.LIB_DIR}/{lib_subfolder}'
        target_subfolder_path = f'{lib_target_folder}/{lib_subfolder}'
        for sub_file in os.listdir(subfolder_path):
            source_path = f'{subfolder_path}/{sub_file}'
            target_path = f'{target_subfolder_path}/{sub_file}'
            if os.path.isfile(source_path):
                shutil.copy(source_path, target_path)

    consts_path = f'{settings.LIB_DIR}/sdsc/consts'
    target_consts_path = f'{lib_target_folder}/sdsc/consts'
    for const_file in os.listdir(consts_path):
        source_path = f'{consts_path}/{const_file}'
        target_path = f'{target_consts_path}/{const_file}'
        if os.path.isfile(source_path):
            shutil.copy(source_path, target_path)

    config = Config(config_path, external_call=True)

    mmseqs2_db_path = config.mmseqs2_db_path
    search_db_base_path = mmseqs2_db_path.rsplit('/', 1)[0]
    pdb_fasta_name = 'pdbba_mmseqs2'

    source_mmseqs_db_path = '%s/%s' % (search_db_base_path, pdb_fasta_name)
    database_target_folder = '%s/StructMAn/lib/base/blast_db' % target_folder
    target_mmseqs_db_path = '%s/%s' % (database_target_folder, pdb_fasta_name)
    if os.path.isfile(target_mmseqs_db_path):
        os.remove(target_mmseqs_db_path)
    shutil.copy(source_mmseqs_db_path, target_mmseqs_db_path)
    print('Source for mmseqs db:', source_mmseqs_db_path, 'Target for mmseqs db:', target_mmseqs_db_path)

    p = subprocess.Popen(['mmseqs', 'createdb', pdb_fasta_name, 'pdbba_search_db_mmseqs2'], cwd=database_target_folder)
    p.wait()

    if os.path.isfile(target_mmseqs_db_path):
        os.remove(target_mmseqs_db_path)
