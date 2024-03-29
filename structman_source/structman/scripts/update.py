#!/usr/bin/python3
import os
import subprocess
import sys

# if running as script, add local structman package to path
if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__))))))
import structman

from structman.scripts import createPdbBaDb, updateMappingDB, updateAlphafoldModelDB


def main(config, skipUpdatePDB=False, skip_rindb=False, rin_fromScratch=False, update_mapping_db = False, mapping_db_from_scratch = False, update_mapping_db_keep_raw_files = False, update_alphafold_db = False):
    mmseqs_fromScratch = False
    skipStructureDBs = False

    verbose = True

    pdb_path = config.pdb_path
    pdb_update_script = config.pdb_sync_script
    print(pdb_update_script)

    rinerator_base_path = config.rinerator_base_path
    rin_db_path = config.rin_db_path

    mmseqs2_db_path = config.mmseqs2_db_path
    search_db_base_path = mmseqs2_db_path.rsplit('/', 1)[0]
    mmseqs2_tmp = config.mmseqs_tmp_folder
    recently_modified_structures = set()
    if not skipStructureDBs:
        if not skipUpdatePDB:
            if not os.path.exists(pdb_path):
                print('ERROR: Did not found local pdb')
                return
            # Set the BASE_DIR variable in the sync script
            f = open(pdb_update_script, 'r')
            lines = f.readlines()
            f.close()
            newlines = []
            for line in lines:
                if line.count('BASE_DIR=') == 1:
                    line = 'BASE_DIR="%s"\n' % pdb_path
                newlines.append(line)
            f = open(pdb_update_script, 'w')
            f.write(''.join(newlines))
            f.close()

            # update local pdb
            if verbose:
                p = subprocess.Popen([pdb_update_script])
            else:
                FNULL = open(os.devnull, 'w')
                p = subprocess.Popen([pdb_update_script], stderr=FNULL, stdout=FNULL)
            p.wait()

            p = subprocess.Popen(['chmod', '777', '-R', '.'], cwd=pdb_path)
            p.wait()

            print('Update PDB done')

        if update_alphafold_db:
            updateAlphafoldModelDB.main(config)

        if not skip_rindb:
            # update rin db
            recently_modified_structures = structman.lib.createRINdb.main(fromScratch=rin_fromScratch, pdb_p=pdb_path, rin_db_path=rin_db_path, n_proc=config.proc_n, rinerator_base_path=rinerator_base_path, process_model_db = config.model_db_active, config = config)

            p = subprocess.Popen(['chmod', '777', '-R', '.'], cwd=rin_db_path)
            p.wait()

            print('Update RIN db done')

            if recently_modified_structures is not None:
                print('Recently modified structures: ', len(recently_modified_structures), recently_modified_structures)

    if not skipUpdatePDB:
        pdb_fasta_name = 'pdbba_mmseqs2'
        config.pdb_fasta_path = f'{search_db_base_path}/{pdb_fasta_name}'
        model_db_fasta_name = 'model_db_mmseqs2'
        config.model_db_fasta_path = f'{search_db_base_path}/{model_db_fasta_name}'

        # update pdbba for mmseqs2
        createPdbBaDb.main(config)

        print("Update search database fasta for MMseqs2 done")

        '''
        #rerun makeblastdb
        p = subprocess.Popen(['makeblastdb','-in','pdbba','-dbtype','prot'],cwd=search_db_base_path)
        p.wait()

        print("Search database for blast created!")
        '''

        # rerun mmseqs2 createdb
        p = subprocess.Popen(['mmseqs', 'createdb', pdb_fasta_name, 'pdbba_search_db_mmseqs2'], cwd=search_db_base_path)
        p.wait()

        p = subprocess.Popen(['rm', '-R', mmseqs2_tmp], cwd=search_db_base_path)
        p.wait()

        if not os.path.isdir(mmseqs2_tmp):
            os.mkdir(mmseqs2_tmp)

        p = subprocess.Popen(['chmod', '777', '-R', mmseqs2_tmp], cwd=search_db_base_path)
        p.wait()

        p = subprocess.Popen(['mmseqs', 'createindex', 'pdbba_search_db_mmseqs2', mmseqs2_tmp, '-s', '7.5'], cwd=search_db_base_path)
        p.wait()

        print("Search database for MMseqs2 created!")

        if config.model_db_fasta_created:
            p = subprocess.Popen(['mmseqs', 'createdb', model_db_fasta_name, 'model_db_search_db_mmseqs2'], cwd=search_db_base_path)
            p.wait()

            p = subprocess.Popen(['rm', '-R', mmseqs2_tmp], cwd=search_db_base_path)
            p.wait()

            if not os.path.isdir(mmseqs2_tmp):
                os.mkdir(mmseqs2_tmp)

            p = subprocess.Popen(['chmod', '777', '-R', mmseqs2_tmp], cwd=search_db_base_path)
            p.wait()

            p = subprocess.Popen(['mmseqs', 'createindex', 'model_db_search_db_mmseqs2', mmseqs2_tmp, '-s', '7.5'], cwd=search_db_base_path)
            p.wait()

            print("Model search database for MMseqs2 created!")


    # update the mapping database
    if update_mapping_db:
        updateMappingDB.main(config, fromScratch = mapping_db_from_scratch, update_mapping_db_keep_raw_files = update_mapping_db_keep_raw_files)

    # update the human proteome mmseqs db, TODO if we want a simple mutation-calling for fasta inputs.

    print('Hurray!')


if __name__ == "__main__":
    config_path = sys.argv[1]
    if not os.path.isfile(config_path):
        print('ERROR: Need path to config file as second argument.')
        sys.exit(1)
    config = structman.Config(config_path, external_call=True)
    main(config)
