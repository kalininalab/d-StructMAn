import random
import string
import os
import subprocess
import sys
import uniprot
import time
import shutil

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

oneToThree = {'C':'CYS', 
'D':'ASP',
'S':'SER',
'V':'VAL',
'Q':'GLN',
'K':'LYS',
'P':'PRO',
'T':'THR',
'F':'PHE',
'A':'ALA',
'H':'HIS',
'G':'GLY',
'I':'ILE',
'L':'LEU',
'R':'ARG',
'W':'TRP',
'N':'ASN',
'Y':'TYR',
'M':'MET',
'E':'GLU'}
#'X':'UNK'}

def parseFasta(path,new_file,lines=None):
    if lines == None:
        f = open(path,'r')
        lines = f.read().split('\n')
        f.close()

    seq_map = {}
    n = 0

    new_lines = []

    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == '>':
            entry_id = line[1:].split()[0]
            seq_map[entry_id] = ''
            n += 1
            new_lines.append(line)
        else:
            
            seq_map[entry_id] += line
            new_lines.append(line)

    f = open(new_file,'w')
    f.write('\n'.join(new_lines))
    f.close()

    return seq_map
    
    
#def parseFastaWithChainId():



def geneSeqMapToFasta(proteins,outfile,config,filtering_db=None):
    lines = []
    n = 0
    m = 0

    u_acs = proteins.get_not_stored_acs()

    if len(u_acs) == 0:
        return False

    for u_ac in u_acs:
        seq = proteins.get_sequence(u_ac)

        if seq == 0 or seq == 1 or seq == '':
            continue
        if filtering_db != None:
            folder_key = gene.split('-')[0][-2:]
            filename = '%s/%s/%s_ref50_gpw.fasta.gz' % (filtering_db,folder_key,u_ac)
            if os.path.isfile(filename):
                n += 1
                continue

        lines.append('>%s' % u_ac)
        lines.append(seq)
        m += 1

    if len(lines) > 0:
        if config.verbosity >= 2:
            print('Filtered ',n,' Proteins before mmseqs2')
            print(m,' sequences going into mmseqs2')
        f = open(outfile,'w')
        f.write('\n'.join(lines))
        f.close()
        return True
    else: 
        return 'Empty fasta file'

##Ask about that
def getSequence(u_ac):
    db_adress = 'bioinfodb'
    db_user_name = 'agress'
    db_password = '3GccJS8c'
    try:
        db = MySQLdb.connect(db_adress,db_user_name,db_password,'struct_man_db_uniprot')
        cursor = db.cursor()
    except:
        db = None
        cursor = None
    gene_sequence_map = uniprot.getSequencesPlain([u_ac],db,cursor)
    if db != None:
        db.close()
    return gene_sequence_map

def parseHits(temp_outfile,option_seq_thresh,small_genes):
    f = open(temp_outfile,'r')
    lines = f.read().split('\n')
    f.close()

    entries = {}

    if len(lines) == 0:
        return(entries)

    pdb_ids = set()

    for line in lines:
        if line == '':
            continue
        words = line.split()
        gene = words[0]
        hitlist = words[1].split(',')
        seq_id = 100.0*float(words[2])
        aln_length = int(words[3])
        
        if aln_length < 50:
            if not gene in small_genes:
                continue
            elif seq_id < (option_seq_thresh*2):
                continue
        #if seq_id < option_seq_thresh:
        #    print seq_id
        #    print hitlist
        #    continue

        if not gene in entries:
            entries[gene] = {}

        hits = {}

        for hit in hitlist:
            pdb,chain = hit.split('-')
            if pdb not in hits:
                hits[pdb] = [chain,set([chain])]
            else:
                hits[pdb][1].add(chain)

        for hit in hits:
            pdb_id = hit
            pdb_ids.add(pdb_id)
            chain = hits[hit][0]
            oligos = hits[hit][1]
            if not len(chain) > 1:
                if not (pdb_id,chain) in entries[gene]:
                    entries[gene][(pdb_id,chain)] = {"Seq_Id":seq_id,"Coverage":'-',"Oligo":oligos,"Length":aln_length}
                else:
                    if aln_legnth > entries[gene][(pdb_id,chain)]["Length"]:
                        entries[gene][pdb_id] = {"Seq_Id":seq_id,"Coverage":'-',"Oligo":oligos,"Length":aln_length}
                    entries[gene][(pdb_id,chain)]["Oligo"].update(oligos)

    return entries,pdb_ids

#called by serializePipeline
def search(proteins,config):
    mmseqs_tmp_folder = '%s/mmseqs_tmp' % config.mmseqs_tmp_folder
    if not os.path.exists(mmseqs_tmp_folder):
        os.mkdir(mmseqs_tmp_folder)

    mmseqs2_path = config.mmseqs2_path
    search_db = config.mmseqs2_db_path
    option_seq_thresh = config.option_seq_thresh

    small_proteins = set()

    u_acs = proteins.get_protein_u_acs()
    for u_ac in u_acs:
        try:
            if len(proteins.get_sequence(u_ac)) < 100:
                small_proteins.add(u_ac)
        except:
            continue

    t0 = time.time()

    temp_fasta = '%s/tmp_%s.fasta' % (mmseqs_tmp_folder,randomString())
    to_fasta_out = geneSeqMapToFasta(proteins,temp_fasta,config)

    if isinstance(to_fasta_out,str):
        config.errorlog.add_warning('%s , mmseqs2 skipped, %s' % (to_fasta_out,str(list(u_acs)[:10])))
        return {},set()
    if not to_fasta_out:
        #All proteins are stored, no need for mmseqs
        return {},set()

    temp_outfile = '%s/tmp_outfile_%s.fasta' % (mmseqs_tmp_folder,randomString())

    t1 = time.time()

    if config.verbosity >= 1:
        print(mmseqs2_path,'easy-search',temp_fasta,search_db,temp_outfile,mmseqs_tmp_folder)

    if config.verbosity >= 3:
        if len(small_proteins) == 0:
            p = subprocess.Popen([mmseqs2_path,'easy-search',temp_fasta,search_db,temp_outfile,mmseqs_tmp_folder,'--max-seqs','999999','--min-aln-len','50','-s','7.5'])
        else:
            p = subprocess.Popen([mmseqs2_path,'easy-search',temp_fasta,search_db,temp_outfile,mmseqs_tmp_folder,'--max-seqs','999999','-s','7.5'])
        p.wait()
    else:
        FNULL = open(os.devnull, 'w')
        if len(small_proteins) == 0:
            p = subprocess.Popen([mmseqs2_path,'easy-search',temp_fasta,search_db,temp_outfile,mmseqs_tmp_folder,'--max-seqs','999999','--min-aln-len','50','-s','7.5'],stdout=FNULL)
        else:
            p = subprocess.Popen([mmseqs2_path,'easy-search',temp_fasta,search_db,temp_outfile,mmseqs_tmp_folder,'--max-seqs','999999','-s','7.5'],stdout=FNULL)
        p.wait()

    os.remove(temp_fasta)

    t2 = time.time()

    hits,pdb_ids = parseHits(temp_outfile,option_seq_thresh,small_proteins)

    for u_ac in small_proteins:
        if u_ac.count(':') == 1:
            if not u_ac in hits:
                pdb_id,chain = u_ac.split(':')
                hits[u_ac] = {(pdb_id,chain):{'Seq_Id':100.0,'Length':len(proteins.get_sequence(u_ac)),'Coverage':'-','Oligo':[chain]}}
                pdb_ids.add(pdb_id)

    os.remove(temp_outfile)

    for fn in os.listdir(mmseqs_tmp_folder):
        subfolder_path = '%s/%s' % (mmseqs_tmp_folder,fn)
        if os.path.exists(subfolder_path):
            if os.path.getmtime(subfolder_path) > config.prog_start_time:
                try:
                    shutil.rmtree(subfolder_path)
                except:
                    config.errorlog.add_warning('Tmp folder wipe failed for: %s' % subfolder_path)

    t3 = time.time()

    if config.verbosity >= 2:
        print("MMseqs2 Part 1: %s" % (str(t1-t0)))
        print("MMseqs2 Part 2: %s" % (str(t2-t1)))
        print("MMseqs2 Part 3: %s" % (str(t3-t2)))

    return hits,pdb_ids


