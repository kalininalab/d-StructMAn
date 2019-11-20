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
    print(n)

    f = open(new_file,'w')
    f.write('\n'.join(new_lines))
    f.close()

    return seq_map

def geneSeqMapToFasta(gene_seq_map,outfile,filtering_db=None):
    lines = []
    n = 0
    m = 0

    for gene in gene_seq_map:
        seq = gene_seq_map[gene][0]
        if seq == 0 or seq == 1 or seq == '':
            continue
        if filtering_db != None:
            folder_key = gene.split('-')[0][-2:]
            filename = '%s/%s/%s_ref50_gpw.fasta.gz' % (filtering_db,folder_key,gene)
            if os.path.isfile(filename):
                n += 1
                continue

        lines.append('>%s' % gene)
        lines.append(seq)
        m += 1

    if len(lines) > 0:
        print('Filtered ',n,' Proteins before mmseqs2')
        print(m,' sequences going into mmseqs2')
        f = open(outfile,'w')
        f.write('\n'.join(lines))
        f.close()
        return None
    else:
        return 'Empty fasta file'

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

def search(gene_seq_map,search_db,mmseqs2_path,mmseqs_tmp_folder,option_seq_thresh,verbose=False):

    #print gene_seq_map

    small_genes = set()
    for gene in gene_seq_map:
        try:
            if len(gene_seq_map[gene]) < 100:
                small_genes.add(gene)
        except:
            continue

    t0 = time.time()

    temp_fasta = '%s/tmp_%s.fasta' % (mmseqs_tmp_folder,randomString())
    error_message = geneSeqMapToFasta(gene_seq_map,temp_fasta)

    if error_message != None:
        print(error_message,', mmseqs2 skiped')
        return {},set()

    temp_outfile = '%s/tmp_outifle_%s.fasta' % (mmseqs_tmp_folder,randomString())

    t1 = time.time()

    FNULL = open(os.devnull, 'w')
    if len(small_genes) == 0:
        p = subprocess.Popen([mmseqs2_path,'easy-search',temp_fasta,search_db,temp_outfile,mmseqs_tmp_folder,'--max-seqs','999999','--min-aln-len','50'],stdout=FNULL)
    else:
        p = subprocess.Popen([mmseqs2_path,'easy-search',temp_fasta,search_db,temp_outfile,mmseqs_tmp_folder,'--max-seqs','999999'],stdout=FNULL)
    p.wait()

    os.remove(temp_fasta)

    t2 = time.time()

    hits,pdb_ids = parseHits(temp_outfile,option_seq_thresh,small_genes)

    for gene in small_genes:
        if gene.count(':') == 1:
            if not gene in hits:
                pdb_id,chain = gene.split(':')
                hits[gene] = {(pdb_id,chain):{'Seq_Id':100.0,'Length':len(gene_seq_map[gene]),'Coverage':'-','Oligo':set([chain])}}
                pdb_ids.add(pdb_id)

    os.remove(temp_outfile)

    shutil.rmtree(mmseqs_tmp_folder)
    os.mkdir(mmseqs_tmp_folder)

    t3 = time.time()

    if verbose:
        print("MMseqs2 Part 1: %s" % (str(t1-t0)))
        print("MMseqs2 Part 2: %s" % (str(t2-t1)))
        print("MMseqs2 Part 3: %s" % (str(t3-t2)))

    #print hits

    return hits,pdb_ids

