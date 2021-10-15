import os
import random
import shutil
import string
import subprocess
import time
import pymysql as MySQLdb

from structman.lib import uniprot


def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))


def geneSeqMapToFasta(proteins, outfile, config, filtering_db=None):
    lines = []
    n = 0
    m = 0

    prot_ids = proteins.get_protein_ids()

    n_not_stored = 0

    for prot_id in prot_ids:

        if proteins[prot_id].stored and prot_id not in proteins.indels:
            continue
        n_not_stored += 1

        seq = proteins.get_sequence(prot_id)

        if seq == 0 or seq == 1 or seq == '':
            continue
        if filtering_db is not None:
            folder_key = prot_id.split('-')[0][-2:]
            filename = '%s/%s/%s_ref50_gpw.fasta.gz' % (filtering_db, folder_key, prot_id)
            if os.path.isfile(filename):
                n += 1
                continue

        lines.append('>%s' % prot_id)
        lines.append(seq)
        m += 1

    if n_not_stored == 0:
        return False

    if len(lines) > 0:
        if config.verbosity >= 2:
            print('Filtered ', n, ' Proteins before mmseqs2')
            print(m, ' sequences going into mmseqs2')
        f = open(outfile, 'w')
        f.write('\n'.join(lines))
        f.close()
        return True
    else:
        return 'Empty fasta file'


def parseHits(temp_outfile, option_seq_thresh, small_genes):
    f = open(temp_outfile, 'r')
    lines = f.read().split('\n')
    f.close()

    entries = {}

    if len(lines) == 0:
        return (entries)

    pdb_ids = set()

    for line in lines:
        if line == '':
            continue
        words = line.split()
        gene = words[0]
        hitlist = words[1].split(',')
        seq_id = 100.0 * float(words[2])
        aln_length = int(words[3])
        target_len = int(words[4])
        coverage = float(words[5])

        if aln_length < 50:
            if gene not in small_genes:
                continue
            elif seq_id < (option_seq_thresh * 2):
                continue

        if gene not in entries:
            entries[gene] = {}

        hits = {}

        for hit in hitlist:
            pdb, chain = hit.split('-')
            if pdb not in hits:
                hits[pdb] = [chain, set([chain])]
            else:
                hits[pdb][1].add(chain)

        for hit in hits:
            pdb_id = hit
            pdb_ids.add(pdb_id)
            chain = hits[hit][0]
            oligos = hits[hit][1]
            if not len(chain) > 1:
                if not (pdb_id, chain) in entries[gene]:
                    entries[gene][(pdb_id, chain)] = [seq_id, coverage, oligos, aln_length, target_len]
                else:
                    if aln_length > entries[gene][(pdb_id, chain)][3]:
                        entries[gene][pdb_id] = [seq_id, coverage, oligos, aln_length, target_len]
                    entries[gene][(pdb_id, chain)][2].update(oligos)

    return entries, pdb_ids


# called by serializePipeline
def search(proteins, config):
    mmseqs_tmp_folder = '%s/mmseqs_tmp' % config.mmseqs_tmp_folder
    if not os.path.exists(mmseqs_tmp_folder):
        os.mkdir(mmseqs_tmp_folder)

    mmseqs2_path = config.mmseqs2_path
    search_db = config.mmseqs2_db_path
    option_seq_thresh = config.option_seq_thresh

    small_proteins = set()

    u_acs = proteins.get_protein_ids()
    for u_ac in u_acs:
        try:
            if len(proteins.get_sequence(u_ac)) < 100:
                small_proteins.add(u_ac)
        except:
            continue

    t0 = time.time()

    temp_fasta = '%s/tmp_%s.fasta' % (mmseqs_tmp_folder, randomString())
    to_fasta_out = geneSeqMapToFasta(proteins, temp_fasta, config)

    if isinstance(to_fasta_out, str):
        config.errorlog.add_warning('%s , mmseqs2 skipped, %s' % (to_fasta_out, str(list(u_acs)[:10])))
        return {}, set()
    if not to_fasta_out:
        # All proteins are stored, no need for mmseqs
        if config.verbosity >= 3:
            print('No need for sequence similarity search, all proteins are stored already.')
        return {}, set()

    temp_outfile = '%s/tmp_outfile_%s.fasta' % (mmseqs_tmp_folder, randomString())

    t1 = time.time()

    if config.verbosity >= 2:
        print(mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder)

    out_format_str = 'query,target,fident,alnlen,tlen,qcov'

    if config.verbosity >= 3:
        if len(small_proteins) == 0:
            cmds = [mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder, '--format-output', out_format_str,
                    '--max-seqs', '999999', '--min-aln-len', '50', '-s', '7.5', '--split-memory-limit', '%1.0fG' % (0.5 * config.gigs_of_ram)]
            p = subprocess.Popen(cmds)
        else:
            cmds = [mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder, '--format-output', out_format_str,
                    '--max-seqs', '999999', '-s', '7.5', '--split-memory-limit', '%1.0fG' % (0.5 * config.gigs_of_ram)]
            p = subprocess.Popen(cmds)
        p.wait()
    else:
        FNULL = open(os.devnull, 'w')
        if len(small_proteins) == 0:
            cmds = [mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder, '--format-output', out_format_str,
                    '--max-seqs', '999999', '--min-aln-len', '50', '-s', '7.5', '--split-memory-limit', '%1.0fG' % (0.5 * config.gigs_of_ram)]
            p = subprocess.Popen(cmds, stdout=FNULL)
        else:
            cmds = [mmseqs2_path, 'easy-search', temp_fasta, search_db, temp_outfile, mmseqs_tmp_folder, '--format-output', out_format_str,
                    '--max-seqs', '999999', '-s', '7.5', '--split-memory-limit', '%1.0fG' % (0.5 * config.gigs_of_ram)]
            p = subprocess.Popen(cmds, stdout=FNULL)
        p.wait()

    os.remove(temp_fasta)

    t2 = time.time()

    hits, pdb_ids = parseHits(temp_outfile, option_seq_thresh, small_proteins)

    for u_ac in small_proteins:
        if u_ac.count(':') == 1:
            if u_ac not in hits:
                pdb_id, chain = u_ac.split(':')
                hits[u_ac] = {(pdb_id, chain): [100.0, len(proteins.get_sequence(u_ac)), 1.0, [chain]]}
                pdb_ids.add(pdb_id)

    os.remove(temp_outfile)

    for fn in os.listdir(mmseqs_tmp_folder):
        subfolder_path = '%s/%s' % (mmseqs_tmp_folder, fn)
        if os.path.exists(subfolder_path):
            if os.path.getmtime(subfolder_path) > config.prog_start_time:
                try:
                    shutil.rmtree(subfolder_path)
                except:
                    if config.verbosity >= 3:
                        config.errorlog.add_warning('Tmp folder wipe failed for: %s' % subfolder_path)

    t3 = time.time()

    if config.verbosity >= 2:
        print("MMseqs2 Part 1: %s" % (str(t1 - t0)))
        print("MMseqs2 Part 2: %s" % (str(t2 - t1)))
        print("MMseqs2 Part 3: %s" % (str(t3 - t2)))

    return hits, pdb_ids
