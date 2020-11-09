#!/usr/bin/python3
import os
import gzip
import subprocess
import sys
import traceback
import gzip
import zlib
import resource
import centrality
import time
from multiprocessing import Process, Queue, Manager, Value, Lock

#forbidden symbols: : _ ( ) . -

#lowercase_map = {'a':'!','b':'@','c':'#'}#,'d':'$','e':'%','f':'^','g':'&','h':'*','i':'+','j':'=','k':'{','l':'}','m':'[','n':']','o':';','p':'<','q':'>','r':',','s':'?','t':'`','u':'~','v':'|','w':'/','x':'"','y':"'"}

#lowercase_order = {y:x for x,y in lowercase_map.iteritems()}
import pdbParser



def parsePDB(page):

    page = page.decode('ascii')


    lines = page.split('\n')

    original_chains = set([])

    ligands = set()
    modres_map = {}

    new_lines = []

    #changed_chains = []
    firstAltLoc = None
    for line in lines:
        if len(line) > 26:
            record_name = line[0:6].replace(" ","")
            atom_nr = line[6:11].replace(" ","")
            atom_name = line[12:16].replace(" ","")
            res_name = line[17:20].replace(" ","")
            chain_id = line[21]
            
            #if chain_id in lowercase_map:
            #    chain_id = lowercase_map[chain_id]
            #    changed_chains.append(chain_id)

            line = '%s%s%s' % (line[:21],chain_id,line[22:])

            res_nr = line[22:27].replace(" ","")
            if record_name == 'ENDMDL':
                break
            if record_name == 'MODRES':
                chain_id = line[16]
                res_nr = line[18:23].replace(" ","")
                res_name = line[24:27].replace(" ","")
                if not (chain_id,res_nr) in modres_map:
                    modres_map[(chain_id,res_nr)] = res_name
            if record_name == "ATOM" or record_name == "HETATM":
                altLoc = line[16]
                if firstAltLoc == None and altLoc != ' ':
                    firstAltLoc = altLoc #The first found alternative Location ID is set as the major alternative location ID or firstAltLoc
                if altLoc != ' ' and altLoc != firstAltLoc: #Whenever an alternative Location ID is found, which is not firstAltLoc, then skip the line
                    continue
                if len(res_name) < 3:
                    ligands.add((res_name,chain_id,res_nr))
                else:
                    if record_name == "HETATM":
                        if not (chain_id,res_nr) in modres_map and not res_name in pdbParser.threeToOne:
                            ligands.add((res_name,chain_id,res_nr))
                        else:
                            original_chains.add(chain_id)
                    else:
                        original_chains.add(chain_id)
        new_lines.append(line)

    page = '\n'.join(new_lines)
    #return original_chains,ligands,page,changed_chains
    return original_chains,ligands,page

'''
def changeBackChains(changed_chains,n_sif_file,n_intsc_file,n_nrint_file,n_res_file):
    for fn in [n_sif_file,n_intsc_file,n_nrint_file,n_res_file]:
        f = open(fn,'r')
        page = f.read()
        f.close()

        for special_chain_id in changed_chains:
            original_chain = lowercase_order[special_chain_id]
            page.replace(special_chain_id,original_chain)

        f = open(fn,'w')
        f.write(page)
        f.close()

    return
'''

def calcRIN(page,out_path,pdb_id,rinerator_path,remove_tmp_files,verbosity):
    original_chains,ligands,page = parsePDB(page)

    if len(original_chains) == 0:
        return
    tmp_pdb = "%s/%s.pdb" % (out_path,pdb_id)
    tmp_chain = "%s/chains_%s.txt" % (out_path,pdb_id)
    tmp_ligands = "%s/ligands_%s.txt" % (out_path,pdb_id)

    f = open(tmp_pdb,'w')
    f.write(page)
    f.close()

    f = open(tmp_chain,'w')
    f.write(",".join(original_chains))
    f.close()

    f = open(tmp_ligands,'w')
    f.write("\n".join([','.join(x) for x in ligands]))
    f.close()
    rinerator_base_path = rinerator_path.rsplit('/',1)[0]
    sys.path.append(rinerator_base_path)
    import get_chains

    reduce_cmd = '%s/reduce' % rinerator_base_path
    probe_cmd = '%s/probe' % rinerator_base_path


    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')
    get_chains.main(tmp_pdb,out_path,tmp_chain,tmp_ligands,True,reduce_cmd,probe_cmd)
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

    '''
    cmd = ['python3',rinerator_path,'-d',tmp_pdb,out_path,tmp_chain,tmp_ligands]

    FNULL = open(os.devnull, 'w')
    if verbosity >= 2:
        p = subprocess.Popen(cmd,stdout=FNULL)
        outs,errs = p.communicate()
        print('Rinerator error output for ',pdb_id,'\n',errs)
    else:
        p = subprocess.Popen(cmd,stderr=FNULL,stdout=FNULL)
        p.communicate()
    '''

    if remove_tmp_files:
        os.remove(tmp_pdb)
        os.remove(tmp_chain)
        os.remove(tmp_ligands)
    else:
        print(tmp_pdb)
        print(tmp_chain)
        print(tmp_ligands)

    reduce_file = "%s/%s_h.ent" % (out_path,pdb_id)
    probe_file = "%s/%s_h.probe" % (out_path,pdb_id)
    sif_file = "%s/%s_h.sif" % (out_path,pdb_id)
    n_sif_file = "%s/%s.sif" % (out_path,pdb_id)

    intsc_file = "%s/%s_h_intsc.ea" % (out_path,pdb_id)
    n_intsc_file = "%s/%s_intsc.ea" % (out_path,pdb_id)
    nrint_file = "%s/%s_h_nrint.ea" % (out_path,pdb_id)
    n_nrint_file = "%s/%s_nrint.ea" % (out_path,pdb_id)
    res_file = "%s/%s_h_res.txt" % (out_path,pdb_id)
    n_res_file = "%s/%s_res.txt" % (out_path,pdb_id)

    '''
    if p.returncode != 0:
        try:
            os.remove(reduce_file)
            os.remove(probe_file)
            os.remove(sif_file)
            os.remove(intsc_file)
            os.remove(nrint_file)
            os.remove(res_file)
        except:
            pass
        raise IOError("RINerator failed")
    '''

    os.remove(reduce_file)
    os.remove(probe_file)

    os.rename(sif_file,n_sif_file)
    os.rename(intsc_file,n_intsc_file)
    os.rename(nrint_file,n_nrint_file)
    os.rename(res_file,n_res_file)

    #if changed_chains != []:
    #    changeBackChains(changed_chains,n_sif_file,n_intsc_file,n_nrint_file,n_res_file)

    if os.path.isfile("%s.gz" % n_sif_file):
        os.remove("%s.gz" % n_sif_file)
    if os.path.isfile("%s.gz" % n_intsc_file):
        os.remove("%s.gz" % n_intsc_file)
    if os.path.isfile("%s.gz" % n_nrint_file):
        os.remove("%s.gz" % n_nrint_file)
    if os.path.isfile("%s.gz" % n_res_file):
        os.remove("%s.gz" % n_res_file)

    centrality.main(n_sif_file)

    os.system("gzip %s" % n_sif_file)
    os.system("gzip %s" % n_intsc_file)
    os.system("gzip %s" % n_nrint_file)
    os.system("gzip %s" % n_res_file)
    return

def createRinProc(in_queue,lock,fromScratch,i,forceCentrality,remove_tmp_files,base_path,rinerator_path,errorlog):
    with lock:
        in_queue.put(None)
    while True:
        with lock:
            in_t = in_queue.get()
            if in_t == None:
                #print('Terminate Rinerator Process: ',i)
                return
            (pdbgz_path,pdb_id,recently_modified) = in_t
        try:
            out_path = "%s/%s/%s" % (base_path,pdb_id[1:-1],pdb_id)
            if not os.path.exists(out_path):
                os.mkdir(out_path)

            n_sif_file = "%s/%s.sif" % (out_path,pdb_id)
            #skip, if network is already there
            if os.path.isfile("%s.gz" % n_sif_file) and not fromScratch and not recently_modified:
                if forceCentrality:
                    centrality.main(n_sif_file)
                continue

            f = gzip.open(pdbgz_path, 'rb')
            page = f.read()
            f.close()
            #original_chains,ligands,page,changed_chains = parsePDB(page)

            calcRIN(page,out_path,pdb_id,rinerator_path,remove_tmp_files,0)

        except:
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            print("Error: ",e,f,g)
            error = '\n'.join([pdb_id,pdbgz_path,str(e),str(f),str(g)])
            errortext  = "###############################################################################\n%s\n###############################################################################\n" % (error)

            with lock:
                f = open(errorlog,'a')
                f.write(errortext)
                f.close()

def test_single_file(pdb_id):
    calculateRINsFromPdbList([pdb_id],remove_tmp_files=False)

def calculateRINsFromPdbList(pdbs,fromScratch=True,forceCentrality=True,remove_tmp_files=True,n_proc=32):

    pdbs = set([x.lower() for x in pdbs])

    lim = 100*1024*1024*1024

    resource.setrlimit(resource.RLIMIT_AS, (lim, lim))

    if not os.path.isfile(het_dict_path):
        print("%s not found" % het_dict_path)
        sys.exit(1)

    os.environ["REDUCE_HET_DICT"] = het_dict_path

    num_of_proc = n_proc

    manager = Manager()
    lock = manager.Lock()

    in_queue = manager.Queue()

    bio_pdbs = set()

    total_structures = 0

    subfolders = os.listdir(bio_assembly_path)
    for subfolder in subfolders:
        sub_path = "%s/%s" % (bio_assembly_path,subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists("%s/%s" % (base_path,subfolder)):
            os.mkdir("%s/%s" % (base_path,subfolder))

        for fn in files:
            if fn.count('.pdb1.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path,fn)
                if os.path.getsize(pdbgz_path) > 50*1024*1024:
                    continue
                pdb_id = fn.replace('.pdb1.gz','')
                if not pdb_id in pdbs:
                    continue

                bio_pdbs.add(pdb_id)
                in_queue.put((pdbgz_path,pdb_id))
                total_structures += 1

    subfolders = os.listdir(AU_path)
    for subfolder in subfolders:

        sub_path = "%s/%s" % (AU_path,subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists("%s/%s" % (base_path,subfolder)):
            os.mkdir("%s/%s" % (base_path,subfolder))

        for fn in files:
            if fn.count('.ent.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path,fn)
                if os.path.getsize(pdbgz_path) > 50*1024*1024:
                    continue
                pdb_id = fn[3:7]
                if not '%s_au' % pdb_id in pdbs:
                    continue
                if pdb_id in bio_pdbs:
                    continue
                in_queue.put((pdbgz_path,pdb_id))
                total_structures += 1

    print('Amount of structures for RINerator: ',total_structures)

    processes = {}
    for i in range(1,num_of_proc + 1):
        p = Process(target=createRinProc, args=(in_queue,lock,fromScratch,i,forceCentrality,remove_tmp_files,base_path,rinerator_path,errorlog))
        processes[i] = p
        print('Start RINerator Process: ',i)
        p.start()
    for i in processes:
        processes[i].join()
    return


def main(fromScratch=False,forceCentrality=False,update_days=2.,pdb_p='',rin_db_path='',n_proc=32,rinerator_base_path = ''):

    pdb_path=pdb_p
    bio_assembly_path = '%s/data/biounit/PDB/divided' % pdb_path
    AU_path = "%s/data/structures/divided/pdb" % pdb_path

    rinerator_path = "%s/get_chains.py" % rinerator_base_path
    het_dict_path = "%s/reduce_wwPDB_het_dict.txt" % rinerator_base_path

    base_path = rin_db_path

    errorlog = "%s/createRINdb_errorlog.txt" % rinerator_base_path

    starttime = time.time()
    update_window = 60.*60.*24.*update_days

    lim = 100*1024*1024*1024

    resource.setrlimit(resource.RLIMIT_AS, (lim, lim))

    if not os.path.isfile(het_dict_path):
        print("%s not found" % het_dict_path)
        sys.exit(1)

    os.environ["REDUCE_HET_DICT"] = het_dict_path

    num_of_proc = n_proc

    manager = Manager()
    lock = manager.Lock()

    in_queue = manager.Queue()

    bio_pdbs = set()

    subfolders = os.listdir(bio_assembly_path)

    recently_modified_structures = set()

    N = 0

    for subfolder in subfolders:

        sub_path = "%s/%s" % (bio_assembly_path,subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists("%s/%s" % (base_path,subfolder)):
            os.mkdir("%s/%s" % (base_path,subfolder))

        for fn in files:
            if fn.count('.pdb1.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path,fn)
                if os.path.getsize(pdbgz_path) > 50*1024*1024:
                    continue

                modification_time = os.path.getmtime(pdbgz_path)

                pdb_id = fn.replace('.pdb1.gz','')

                if starttime - modification_time < update_window:
                    recently_modified = True
                    recently_modified_structures.add(pdb_id)
                else:
                    recently_modified = False

                bio_pdbs.add(pdb_id)
                in_queue.put((pdbgz_path,pdb_id,recently_modified))
                N += 1

    subfolders = os.listdir(AU_path)
    for subfolder in subfolders:


        sub_path = "%s/%s" % (AU_path,subfolder)
        files = os.listdir(sub_path)
        if not os.path.exists("%s/%s" % (base_path,subfolder)):
            os.mkdir("%s/%s" % (base_path,subfolder))

        for fn in files:
            if fn.count('.ent.gz') == 1:
                pdbgz_path = "%s/%s" % (sub_path,fn)
                if os.path.getsize(pdbgz_path) > 50*1024*1024:
                    continue
                pdb_id = fn[3:7]
                if pdb_id in bio_pdbs:
                    continue

                modification_time = os.path.getmtime(pdbgz_path)
                if starttime - modification_time < update_window:
                    recently_modified = True
                    recently_modified_structures.add(pdb_id)
                else:
                    recently_modified = False

                in_queue.put((pdbgz_path,pdb_id,recently_modified))
                N += 1

    print('Creating RINs for',N,'structures.')

    processes = {}
    for i in range(1,min([num_of_proc,N]) + 1):
        p = Process(target=createRinProc, args=(in_queue,lock,fromScratch,i,forceCentrality,True,base_path,rinerator_path,errorlog))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()
    return recently_modified_structures
if __name__ == "__main__":
    #pdb_id = '1abr'
    #test_single_file(pdb_id)
    main(fromScratch=True)

