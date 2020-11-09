#!/usr/bin/python3
import subprocess
import os
import sys
import structman

def main(config,skipUpdatePDB = False,skip_rindb = False):
    main_file_path = (os.path.abspath(sys.argv[0])).rsplit('/',1)[0]
    rin_fromScratch = False
    forceCentrality = False
    mmseqs_fromScratch = False
    skipStructureDBs = False
    
    verbose = True

    n_proc = config.proc_n

    pdb_path = config.pdb_path
    pdb_update_script = config.pdb_sync_script

    rinerator_base_path = config.rinerator_base_path
    rin_db_path = config.rin_db_path

    search_db_base_path = config.blast_db_path
    mmseqs2_db_path = config.mmseqs2_db_path
    mmseqs2_tmp = config.mmseqs_tmp_folder

    if not skipStructureDBs:
        if not skipUpdatePDB:
            if not os.path.exists(pdb_path):
                print('ERROR: Did not found local pdb')
                return
            #Set the BASE_DIR variable in the sync script
            f = open(pdb_update_script,'r')
            lines = f.readlines()
            f.close()
            newlines = []
            for line in lines:
                if line.count('BASE_DIR=') == 1:
                    line = 'BASE_DIR="%s"\n' % pdb_path
                newlines.append(line)
            f = open(pdb_update_script,'w')
            f.write(''.join(newlines))
            f.close()

            #update local pdb
            if verbose:
                p = subprocess.Popen([pdb_update_script])
            else:
                FNULL = open(os.devnull, 'w')
                p = subprocess.Popen([pdb_update_script],stderr=FNULL,stdout=FNULL)
            p.wait()
            print('Update PDB done')
        if not skip_rindb:
            #update rin db
            sys.path.append(rinerator_base_path)
            import createRINdb
            recently_modified_structures = createRINdb.main(fromScratch=rin_fromScratch,forceCentrality=forceCentrality,update_days=30.,pdb_p=pdb_path,rin_db_path=rin_db_path,n_proc=n_proc,rinerator_base_path = rinerator_base_path)

        print('Update RIN db done')
    else:
        recently_modified_structures = set()

    print('Recently modified structures: ',len(recently_modified_structures),recently_modified_structures)


    #update pdbba for mmseqs2
    import createPdbBaDb
    createPdbBaDb.main(mmseqs2_db_path,recently_modified_structures,fromScratch=mmseqs_fromScratch,pdb_p = pdb_path)

    print("Update search database fasta for MMseqs2 done")

    '''
    #rerun makeblastdb
    p = subprocess.Popen(['makeblastdb','-in','pdbba','-dbtype','prot'],cwd=search_db_base_path)
    p.wait()

    print("Search database for blast created!")
    '''

    #rerun mmseqs2 createdb
    p = subprocess.Popen(['mmseqs','createdb','pdbba_mmseqs2','pdbba_search_db_mmseqs2'],cwd=search_db_base_path)
    p.wait()

    p = subprocess.Popen(['rm','-R',mmseqs2_tmp],cwd=search_db_base_path)
    p.wait()

    if not os.path.isdir(mmseqs2_tmp):
        os.mkdir(mmseqs2_tmp)

    p = subprocess.Popen(['chmod','777','-R',mmseqs2_tmp],cwd=search_db_base_path)
    p.wait()

    p = subprocess.Popen(['mmseqs','createindex','pdbba_search_db_mmseqs2',mmseqs2_tmp,'-s','7.5'],cwd=search_db_base_path)
    p.wait()

    print("Search database for MMseqs2 created!")

    #update the mapping database, TODO

    #update the human proteome mmseqs db, TODO if we want a simple mutation-calling for fasta inputs.

    print('Hurray!')

if __name__ == "__main__":
    config_path = sys.argv[1]
    if not os.path.isfile(config_path):
        print('ERROR: Need path to config file as second argument.')
        sys.exit(1)
    config = structman.Config(config_path,external_call = True)
    main(config)

