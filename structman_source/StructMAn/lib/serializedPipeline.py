import errno
import gc
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
import psutil

import ray

from structman.lib import (
    annovar,
    blast,
    database,
    globalAlignment,
    indel_analysis,
    MMseqs2,
    output,
    pdbParser,
    rin,
    sdsc,
    templateFiltering,
    templateSelection,
    uniprot
)
from structman.utils import ray_init, ray_hack, monitor_ray_store, calculate_chunksizes, pack, unpack

try:
    from memory_profiler import profile
except:
    pass

# Taken from https://stackoverflow.com/questions/2023608/check-what-files-are-open-in-python


def list_fds():
    """List process currently open FDs and their target """
    if not sys.platform.startswith('linux'):
        raise NotImplementedError('Unsupported platform: %s' % sys.platform)

    ret = {}
    base = '/proc/self/fd'
    for num in os.listdir(base):
        path = None
        try:
            path = os.readlink(os.path.join(base, num))
        except OSError as err:
            # Last FD is always the "listdir" one (which may be closed)
            if err.errno != errno.ENOENT:
                raise
        ret[int(num)] = path

    return ret


def SQLDateTime():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())


# For memory profiling, thanks to Fred Cirera
def sizeof_fmt(num, suffix='B'):
    ''' by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified'''
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


def remove_sanity_filtered(config, proteins, indels, primary_protein_id):
    del_list = []
    for pos in proteins[primary_protein_id].positions:
        if not proteins[primary_protein_id].positions[pos].checked:
            del_list.append((pos, proteins[primary_protein_id].positions[pos].wt_aa))
    for pos, wt_aa in del_list:
        remove_position(config, proteins, indels, pos, wt_aa, primary_protein_id)


def remove_position(config, proteins, indels, pos, wt_aa, primary_protein_id):
    warn_text = 'Filtered position through sanity check:%s,%s' % (primary_protein_id, pos)
    if config.verbosity >= 2:
        print(warn_text)
    config.errorlog.add_warning(warn_text)

    if proteins[primary_protein_id].multi_mutations is not None:
        mm_dels = []
        for mm_nr, multi_mutation in enumerate(proteins[primary_protein_id].multi_mutations):
            for mut in multi_mutation:
                if isinstance(mut, tuple):
                    position, aa2 = mut
                    if position.pos == pos and wt_aa == position.wt_aa:
                        mm_dels.append(mm_nr)
        for mm_nr in sorted(mm_dels, reverse=True):
            del proteins[primary_protein_id].multi_mutations[mm_nr]
        if len(proteins[primary_protein_id].multi_mutations) == 0:
            proteins[primary_protein_id].multi_mutations = None

    del proteins[primary_protein_id].positions[pos]

    if len(proteins[primary_protein_id].positions) == 0:
        if primary_protein_id in indels:
            if len(indels[primary_protein_id]) == 0:
                del proteins[primary_protein_id]
        else:
            del proteins[primary_protein_id]

# @profile


def sequenceScan(config, proteins, indels):
    def promote_uac_to_primary_id(u_ac, np_ref):
        primary_protein_id = u_ac
        proteins[primary_protein_id] = proteins[np_ref]

        if primary_protein_id in indels and np_ref in indels:
            indels[primary_protein_id] += indels[np_ref]
        elif np_ref in indels:
            indels[primary_protein_id] = indels[np_ref]

        if primary_protein_id in indels:
            for indel in indels[primary_protein_id]:
                indel.wt_prot = primary_protein_id

        if np_ref in indels:
            del indels[np_ref]
        del proteins[np_ref]
        proteins[primary_protein_id].primary_protein_id = primary_protein_id
        proteins[primary_protein_id].u_ac = u_ac
        return u_ac

    pdb_path = config.pdb_path

    sequenceScanProteins = {}
    sequenceScanPDB = {}
    sequenceScanNM = {}
    sequenceScanNP = {}
    fasta_inputs = {}

    for prot_id in proteins:
        uni_pos, tags = proteins[prot_id].popNone()
        if proteins[prot_id].is_sequence_set():
            if uni_pos:
                fasta_inputs[prot_id] = tags
            continue
        if prot_id.count(':') > 0:
            sequenceScanPDB[prot_id] = tags, uni_pos  # PDB inputs are always processed by the sequence scan, but the positions are only added if uni_pos is true
        elif not sdsc.is_mutant_ac(prot_id):
            if proteins[prot_id].input_id[:2] == 'NM':
                if proteins[prot_id].input_id in indels:  # Indel proteins need complete classification
                    uni_pos = True
                sequenceScanNM[proteins[prot_id].input_id] = tags, uni_pos, proteins[prot_id].u_ac
            elif proteins[prot_id].input_id[:2] == 'NP' or proteins[prot_id].input_id[-2] == '.':
                if proteins[prot_id].input_id in indels:  # Indel proteins need complete classification
                    uni_pos = True
                sequenceScanNP[proteins[prot_id].input_id] = tags, uni_pos, proteins[prot_id].u_ac
            else:
                if prot_id in indels:  # Indel proteins need complete classification
                    uni_pos = True
                sequenceScanProteins[prot_id] = tags, uni_pos  # New: process everything to filter input by sanity checks

    if len(sequenceScanNP) > 0:
        if config.verbosity >= 1:
            print("Amount of proteins going into NP sequenceScan: ", len(sequenceScanNP))
        gene_sequence_map = uniprot.get_refseq_sequences(','.join(sequenceScanNP.keys()), config, seq_type='protein')
        stems = set()
        ref_stem_map = {}
        for np_ref in sequenceScanNP:
            (tags, uni_pos, u_ac) = sequenceScanNP[np_ref]
            if u_ac is None:
                continue
            if u_ac[:3] == 'NP_':
                continue
            u_ac_stem = u_ac.split('_')[0].split('-')[0]
            stems.add(u_ac_stem)
            ref_stem_map[np_ref] = u_ac_stem
        isoform_specific_id_map = uniprot.u_ac_isoform_search(gene_sequence_map, stems, ref_stem_map, config)
        for np_ref in sequenceScanNP:
            (tags, uni_pos, u_ac) = sequenceScanNP[np_ref]
            seq = gene_sequence_map[np_ref]
            if np_ref in isoform_specific_id_map:
                primary_protein_id = promote_uac_to_primary_id(isoform_specific_id_map[np_ref], np_ref)
            else:
                primary_protein_id = np_ref

            proteins[primary_protein_id].sequence = seq

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if seq_pos in proteins[primary_protein_id].positions:
                    proteins[primary_protein_id].positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                    proteins[primary_protein_id].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    position = sdsc.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True)
                    proteins[primary_protein_id].positions[seq_pos] = position

            # sanity check filter
            remove_sanity_filtered(config, proteins, indels, primary_protein_id)

    if len(sequenceScanNM) > 0:
        if config.verbosity >= 1:
            print("Amount of proteins going into NM sequenceScan: ", len(sequenceScanNM))
        gene_sequence_map = uniprot.get_refseq_sequences(','.join(sequenceScanNM.keys()), config)
        ref_stem_map = {}
        stems = set()
        for nm_ref in sequenceScanNM:
            (tags, uni_pos, u_ac) = sequenceScanNM[nm_ref]
            if u_ac is None:
                continue
            if u_ac[:3] == 'NM_':
                continue
            u_ac_stem = u_ac.split('_')[0].split('-')[0]
            stems.add(u_ac_stem)
            ref_stem_map[nm_ref] = u_ac_stem

        isoform_specific_id_map = uniprot.u_ac_isoform_search(gene_sequence_map, stems, ref_stem_map, config)
        for nm_ref in sequenceScanNM:
            (tags, uni_pos, u_ac) = sequenceScanNM[nm_ref]
            seq = gene_sequence_map[nm_ref]
            if nm_ref in isoform_specific_id_map:
                primary_protein_id = promote_uac_to_primary_id(isoform_specific_id_map[nm_ref], nm_ref)
            else:
                primary_protein_id = nm_ref

            proteins[primary_protein_id].sequence = seq

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if seq_pos in proteins[primary_protein_id].positions:
                    proteins[primary_protein_id].positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                    proteins[primary_protein_id].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    position = sdsc.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True)
                    proteins[primary_protein_id].positions[seq_pos] = position

            # sanity check filter
            remove_sanity_filtered(config, proteins, indels, primary_protein_id)

    if len(sequenceScanProteins) > 0:
        if config.verbosity >= 1:
            print("Amount of proteins going into Uniprot sequenceScan: ", len(sequenceScanProteins))

        gene_sequence_map = uniprot.getSequencesPlain(sequenceScanProteins.keys(), config)
        for u_ac in gene_sequence_map:
            if gene_sequence_map[u_ac][0] == 1 or gene_sequence_map[u_ac][0] == 0:
                config.errorlog.add_warning("Error in sequenceScan with gene: %s" % u_ac)
                continue
            seq, disorder_scores, disorder_regions_datastruct = gene_sequence_map[u_ac]
            proteins[u_ac].sequence = seq

            tags, uni_pos = sequenceScanProteins[u_ac]

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                if seq_pos in proteins[u_ac].positions:
                    proteins[u_ac].positions[seq_pos].check(aa, overwrite=config.overwrite_incorrect_wt_aa)
                    proteins[u_ac].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    position = sdsc.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True)
                    proteins[u_ac].positions[seq_pos] = position

            proteins[u_ac].set_disorder_scores(disorder_scores)
            proteins[u_ac].set_disorder_regions(disorder_regions_datastruct)

            # sanity check filter
            remove_sanity_filtered(config, proteins, indels, u_ac)

    if len(sequenceScanPDB) > 0:
        pdb_sequence_map, pdb_pos_map = pdbParser.getSequences(sequenceScanPDB.keys(), pdb_path)
        for u_ac in pdb_sequence_map:
            proteins[u_ac].sequence = pdb_sequence_map[u_ac][0]
            tags, uni_pos = sequenceScanPDB[u_ac]

            residue_id_backmap = {}
            res_pos_map = pdb_pos_map[u_ac]

            for res_id in res_pos_map:
                residue_id_backmap[res_pos_map[res_id]] = res_id
            for (pos, aa) in enumerate(pdb_sequence_map[u_ac][0]):
                seq_pos = pos + 1
                if seq_pos not in residue_id_backmap:
                    continue
                res_id = residue_id_backmap[seq_pos]
                if res_id in proteins[u_ac].res_id_map:
                    proteins[u_ac].res_id_map[res_id].pos = seq_pos
                    proteins[u_ac].positions[seq_pos] = proteins[u_ac].res_id_map[res_id]
                    if tags is None:
                        continue
                    proteins[u_ac].positions[seq_pos].add_tags(tags)
                elif uni_pos:
                    position = sdsc.Position(pos=seq_pos, pdb_res_nr=res_id, wt_aa=aa, tags=tags)
                    proteins[u_ac].positions[seq_pos] = position
                    proteins[u_ac].res_id_map[res_id] = position

    if len(fasta_inputs) > 0:
        for prot_id in fasta_inputs:
            seq = proteins[prot_id].sequence

            tags = fasta_inputs[prot_id]

            for (pos, aa) in enumerate(seq):
                seq_pos = pos + 1
                position = sdsc.Position(pos=seq_pos, wt_aa=aa, tags=tags, checked=True)
                proteins[prot_id].positions[seq_pos] = position

    if config.verbosity >= 3:
        print('Before indel mutation with:', len(indels), 'number of indels')

    for primary_protein_id in indels:
        for indel in indels[primary_protein_id]:
            indel.mutate_sequence(proteins)

    multi_mutation_objects = []

    if not config.only_wt:
        if config.verbosity >= 3:
            print('Before multi mutations mutation')
        for primary_protein_id in list(proteins.keys()).copy():  # iterate over snapshot, since the mutant proteins get added to proteins in the loop
            if proteins[primary_protein_id].multi_mutations is not None:
                # creates mutation protein object
                multi_mutation_objects += proteins[primary_protein_id].create_multi_mutations(proteins, config)

    return proteins, indels, multi_mutation_objects


def processAAChange(aachange, pdb_style=False):
    if ord(aachange[0]) > 47 and ord(aachange[0]) < 58:  # if first char is a number
        aa1 = 'X'
        aachange = "X%s" % aachange
    else:
        aa1 = aachange[0]

    if ord(aachange[-1]) > 47 and ord(aachange[-1]) < 58:  # if last char is a number
        pos = int(aachange[1:])
        aa2 = None
    elif pdb_style and aachange.count('_') == 1:
        aachange = aachange.replace('_', '')
        aa2 = aachange[-1]
        pdb_res_nr = aachange[1:-1]
        return aachange, aa1, aa2, pdb_res_nr
    elif pdb_style:
        aa2 = None
        pdb_res_nr = aachange[1:]
        return aachange, aa1, aa2, pdb_res_nr
    else:
        aa2 = aachange[-1]
        pos = int(aachange[1:-1])

    return aachange, aa1, aa2, pos


def parseFasta(config, nfname):
    f = open(nfname, 'r')
    lines = f.readlines()
    f.close()

    seq_map = {}

    pos_map = {}

    for line in lines:
        line = line[:-1]
        if len(line) == 0:
            continue
        if line[0] == '>':
            words = line[1:].split()
            entry_id = words[0]
            if entry_id.count('|') > 1:
                entry_id = entry_id.split('|')[1]
            seq_map[entry_id] = ''
            if len(words) > 1:
                aacs = words[1]
                aac_tag_tuples = aacs.split(';')
                for aac_tag_tuple in aac_tag_tuples:
                    if aac_tag_tuple.count('<') == 1:
                        aac, tags = aac_tag_tuple.split('<')
                    else:
                        aac = aac_tag_tuple
                        tags = set()

                    positions, multi_mutations = process_mutations_str(aac, tags)
                    add_to_prot_map(pos_map, entry_id, positions, multi_mutations, config)
            else:
                position = sdsc.Position()
                multi_mutations = []
                positions = [position]
                add_to_prot_map(pos_map, entry_id, positions, multi_mutations, config)
        else:
            seq_map[entry_id] += line.replace('\n', '').upper()

    proteins = {}
    indels = {}
    for prot_id in seq_map:
        seq = seq_map[prot_id]

        uniprot.integrate_protein(config, proteins, indels, prot_id, prot_id, pos_map)
        proteins[prot_id].sequence = seq

    outlist = input_chunking(config, proteins, indels)

    return outlist


def add_to_prot_map(prot_map, sp_id, positions, multi_mutations, config):
    if sp_id not in prot_map:
        prot_map[sp_id] = [], [], []
    for indel_or_snv in multi_mutations:
        if not isinstance(indel_or_snv, tuple):
            prot_map[sp_id][1].append(indel_or_snv)
        else:
            prot_map[sp_id][0].append(indel_or_snv[0])

    for position in positions:
        prot_map[sp_id][0].append(position)

    if len(multi_mutations) > 1 and not config.only_snvs:
        prot_map[sp_id][2].append(multi_mutations)
    return prot_map


def process_mutations_str(mutation_str, tags, pdb_style=False):
    if isinstance(tags, str):
        tags = set(tags.split(','))
    multi_mutations = []
    positions = []
    aachanges = mutation_str.split(',')
    for aachange in aachanges:
        if aachange == '':
            continue
        if aachange.count('delins') == 1:
            indel = sdsc.Substitution(raw_str=aachange, tags=tags)
            multi_mutations.append(indel)

        elif aachange.count('del') == 1:
            indel = sdsc.Deletion(raw_str=aachange, tags=tags)
            multi_mutations.append(indel)

        elif aachange.count('ins') == 1:
            indel = sdsc.Insertion(raw_str=aachange, tags=tags)
            multi_mutations.append(indel)

        else:
            indel = None
            aachange, aa1, aa2, pos = processAAChange(aachange, pdb_style=pdb_style)  # if pdb_style is True, then pos is actually pdb res_id

            if pdb_style:
                if aa2 is None:
                    position = sdsc.Position(pdb_res_nr=pos, wt_aa=aa1, mut_aas=set(aa2), tags=tags, mut_tags_map={aa2: tags})
                    positions.append(position)
                else:
                    position = sdsc.Position(pdb_res_nr=pos, wt_aa=aa1, tags=tags)
                    multi_mutations.append((position, aa2))
            else:
                if aa2 is None:
                    position = sdsc.Position(pos=pos, wt_aa=aa1, tags=tags)
                    positions.append(position)
                else:
                    position = sdsc.Position(pos=pos, wt_aa=aa1, mut_aas=set(aa2), tags=tags, mut_tags_map={aa2: tags})
                    multi_mutations.append((position, aa2))
    return positions, multi_mutations


def input_chunking(config, proteins, indels):
    outlist = []

    s = len(proteins)
    if config.verbosity >= 1:
        print("Total proteins: ", s)
    if s > config.chunksize:
        n_of_batches = s // config.chunksize
        if s % config.chunksize != 0:
            n_of_batches += 1
        batchsize = s // n_of_batches
        if s % n_of_batches != 0:
            batchsize += 1

        outlist = []

        rest = s % n_of_batches

        no_more_indels = len(indels) == 0

        for i in range(0, n_of_batches):
            new_map = {}
            indel_map = {}
            prot_count = 0
            for j in range(0, batchsize):
                if prot_count >= batchsize:
                    break
                if len(proteins) == 0:
                    continue
                if no_more_indels:
                    (key, value) = proteins.popitem()
                    new_map[key] = value
                    prot_count += 1
                else:
                    key, indellist = indels.popitem()
                    new_map[key] = proteins[key]
                    del proteins[key]
                    prot_count += 1
                    for indel in indellist:
                        new_map[indel.mut_prot] = proteins[indel.mut_prot]
                        del proteins[indel.mut_prot]
                        prot_count += 1
                    indel_map[key] = indellist
                    no_more_indels = len(indels) == 0
            outlist.append((new_map, indel_map))
        new_map = {}
        while len(proteins) > 0:
            (key, value) = proteins.popitem()
            new_map[key] = value
            if key in indels:
                for indel in indels[key]:
                    new_map[indel.mut_prot] = proteins[indel.mut_prot]
                    del proteins[indel.mut_prot]
                indel_map[key] = indels[key]
        if len(new_map) > 0:
            outlist.append((new_map, indel_map))
    else:
        outlist.append((proteins, indels))
    return outlist


# @profile
def buildQueue(config, filename, already_split=False):
    t0 = time.time()

    proteins = {}
    u_ids = set()
    u_acs = set()
    id_map = {}
    ac_map = {}
    np_map = {}
    nm_map = {}
    hgnc_map = {}

    pdb_map = {}

    if isinstance(filename, str):
        f = open(filename, "r")
        lines = f.read().split('\n')
        f.close()
    else:
        # In case of single line input
        lines = ['\t'.join(filename)]

    if config.low_mem_system and not already_split:
        prot_lines = {}
        prots = []
        for line in lines:
            # skip blank lines and remove c++ style `//` comments
            line = line.split('//')[0].strip()
            if line == '':
                continue
            if len(line) < 3:
                if config.verbosity >= 1:
                    print("Skipped input line:\n%s\nToo short.\n" % line)
                continue
            line = line.replace(' ', '\t')
            words = line.split("\t")
            if len(words) < 1:
                if config.verbosity >= 1:
                    print("Skipped input line:\n%s\nToo few words.\n" % line)
                continue
            sp_id = words[0]  # .replace("'","\\'")
            if sp_id not in prot_lines:
                prot_lines[sp_id] = [line]
                prots.append(sp_id)
            else:
                prot_lines[sp_id].append(line)

        total_num_of_raw_ids = len(prot_lines)
        if total_num_of_raw_ids > (10 * config.chunksize):
            temp_infiles = []
            num_of_infiles = total_num_of_raw_ids // (10 * config.chunksize)
            if total_num_of_raw_ids % (10 * config.chunksize) != 0:
                num_of_infiles += 1

            num_of_prots_per_file = total_num_of_raw_ids // num_of_infiles
            if total_num_of_raw_ids % num_of_infiles != 0:
                num_of_prots_per_file += 1

            for i in range(num_of_infiles):
                temp_file_lines = []
                for j in range(i * num_of_prots_per_file, (i + 1) * num_of_prots_per_file):
                    if j < len(prots):
                        temp_file_lines += prot_lines[prots[j]]
                if len(temp_file_lines) == 0:
                    continue
                temp_file_path = '%s/infile_split_%s.smlf' % (config.temp_folder, str(i))
                f = open(temp_file_path, 'w')
                f.write('\n'.join(temp_file_lines))
                f.close()
                temp_infiles.append(temp_file_path)
            return [], temp_infiles

    tag_map = {}

    for line in lines:
        # skip blank lines and remove c++ style `//` comments
        line = line.split('//')[0].strip()
        if line == '':
            continue
        if len(line) < 3:
            if config.verbosity >= 1:
                print("Skipped input line:\n%s\nToo short.\n" % line)
            continue
        line = line.replace(' ', '\t')
        words = line.split("\t")
        if len(words) < 1:
            if config.verbosity >= 1:
                print("Skipped input line:\n%s\nToo few words.\n" % line)
            continue
        sp_id = words[0]  # .replace("'","\\'")

        if len(sp_id) < 2:
            if config.verbosity >= 1:
                print("Skipped input line:\n%s\nID too short.\n" % line)
            continue

        tags = set()

        if len(words) > 2:
            tags = set(words[2].split(','))

        try:
            if len(words) == 1 or words[1] == '':
                position = sdsc.Position(tags=tags)
                indel = None
                multi_mutations = []
                positions = [position]

            else:
                mutation_str = words[1].replace("\n", "")
                multi_mutations = []
                positions = []
                if (not sp_id.count(':') == 1) or sp_id[0:5] == 'HGNC:':  # this means sp_id is not a pdb-id
                    positions, multi_mutations = process_mutations_str(mutation_str, tags)

                else:  # this means sp_id is a pdb-id
                    positions, multi_mutations = process_mutations_str(mutation_str, tags, pdb_style=True)

        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            config.errorlog.add_error("File Format Error: %s\n%s\n%s\n%s" % (line, str(e), str(f), str(g)))

        if sp_id[2] == "_" or sp_id[-2] == '.':  # RefSeq Ids
            if sp_id[:2] == 'NP':
                np_map = add_to_prot_map(np_map, sp_id, positions, multi_mutations, config)
            elif sp_id[:2] == 'NM':
                nm_map = add_to_prot_map(nm_map, sp_id, positions, multi_mutations, config)
            else:
                #config.errorlog.add_warning('Unsupported Refseq Id:',line)
                # Try to put them all into the NP map
                np_map = add_to_prot_map(np_map, sp_id, positions, multi_mutations, config)
        else:
            if sp_id.count('_') > 0:
                u_ids.add(sp_id)
                id_map = add_to_prot_map(id_map, sp_id, positions, multi_mutations, config)

            elif sp_id[:5] == 'HGNC:':
                hgnc_map = add_to_prot_map(hgnc_map, sp_id, positions, multi_mutations, config)

            elif len(sp_id) == 6 and sp_id.count(':') == 1:
                pdb_chain_tuple = '%s:%s' % (sp_id[:4].upper(), sp_id[-1])  # enforce uppercase pdb-id
                pdb_map = add_to_prot_map(pdb_map, pdb_chain_tuple, positions, multi_mutations, config)

            else:
                u_acs.add(sp_id)
                ac_map = add_to_prot_map(ac_map, sp_id, positions, multi_mutations, config)

    t1 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 1: ", str(t1 - t0), len(ac_map))

    proteins, indels = uniprot.IdMapping(config, ac_map, id_map, np_map, pdb_map, hgnc_map, nm_map)

    t2 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 2: ", str(t2 - t1))

    outlist = input_chunking(config, proteins, indels)

    t3 = time.time()
    if config.verbosity >= 2:
        print("buildQueue Part 3: ", str(t3 - t2))

    return outlist, [None]


def nToAA(seq):
    i = 1
    triple = ''
    aa_seq = ''
    for char in list(seq):
        triple += char
        if i % 3 == 0:
            if triple in sdsc.CODONS:
                aa = sdsc.CODONS[triple] if triple not in sdsc.STOP_CODONS else ''
            aa_seq += aa
            triple = ''
        i += 1
    return aa_seq


# @profile
def getSequences(proteins, config, skip_db=False):

    number_of_processes = config.proc_n
    pdb_path = config.pdb_path
    pdb_input_asymetric_unit = config.pdb_input_asymetric_unit
    blast_processes = config.blast_processes
    cwd = os.getcwd()
    verbose = config.verbose

    t0 = time.time()

    uniprot.getSequences(proteins, config)

    t1 = time.time()
    if config.verbosity >= 2:
        print("Time for getSequences Part 1: %s" % str(t1 - t0))

    if not skip_db:
        database.addProtInfos(proteins, config)

    t2 = time.time()
    if config.verbosity >= 2:
        print("Time for getSequences Part 2: %s" % str(t2 - t1))

    u_acs = proteins.get_protein_ids()

    background_iu_process = None
    if config.iupred_path != '':
        iupred_path = config.iupred_path
        if iupred_path.count('mobidb-lite') > 0:
            mobi_lite = True
        else:
            mobi_lite = False

        t0 = time.time()

        if not mobi_lite:
            iupred_results = []
            iupred_dump = ray.put(config.iupred_path)
        else:
            mobi_list = []

        n_disorder = 0
        stored_disorder_ids = []
        for u_ac in u_acs:
            seq = proteins.get_sequence(u_ac)
            if seq == 0 or seq == 1 or seq is None:
                continue
            disorder_scores = proteins.get_disorder_scores(u_ac)
            if disorder_scores is None:
                if not mobi_lite:
                    iupred_results.append(para_iupred.remote(u_ac, seq, iupred_dump))
                    n_disorder += 1
                else:
                    mobi_list.append('>%s\n%s\n' % (u_ac, seq))
            elif disorder_scores != 'Stored':
                proteins.set_disorder_tool(u_ac, 'MobiDB3.0')

        if not mobi_lite:
            iupred_out = ray.get(iupred_results)

            for iupred_parts in iupred_out:
                proteins.set_disorder_scores(iupred_parts[0], iupred_parts[2])
                proteins.set_disorder_regions(iupred_parts[0], iupred_parts[1])
                proteins.set_disorder_tool(iupred_parts[0], 'IUpred')

        else:
            mobi_tmp_file = 'mobi_tmp_file.fasta'
            f = open(mobi_tmp_file, 'w')
            f.write(''.join(mobi_list))
            f.close()
            mobi_bin_path = '%s/binx/' % iupred_path.rsplit('/', 1)[0]
            mobi_threads = min([7, config.proc_n])
            p = subprocess.Popen([iupred_path, mobi_tmp_file, '-t', str(mobi_threads), '-bin', mobi_bin_path, '-l'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            out, err = p.communicate()
            os.remove(mobi_tmp_file)
            if err != '':
                config.errorlog.add_warning('Warning: mobidb-lite threw an error: %s' % err)
            else:
                entries = []
                for line in out.split('}'):
                    if line.count('{') == 0:
                        continue
                    line = line + '}'
                    mobi_results = json.loads(line)
                    u_ac = mobi_results['acc']
                    raw_scores = mobi_results['p']
                    raw_regions = mobi_results['regions']
                    regions = []
                    for a, b in raw_regions:
                        regions.append([a, b, 'disorder'])
                    scores = {}
                    #seq = proteins.get_sequence(u_ac)
                    #disorder_scores = proteins.get_disorder_scores(u_ac)
                    #disorder_regions = proteins.get_disorder_regions(u_ac)
                    for pos, score in enumerate(raw_scores):
                        scores[pos + 1] = score
                    proteins.set_disorder_scores(u_ac, scores)
                    proteins.set_disorder_regions(u_ac, regions)
                    proteins.set_disorder_tool(u_ac, 'mobidb-lite')

        t1 = time.time()
        if config.verbosity >= 2:
            print("Time for addIupred Part 1: %s" % str(t1 - t0))

        if not skip_db:
            background_iu_process = database.addIupred(proteins, config)

            t2 = time.time()
            if config.verbosity >= 2:
                print("Time for addIupred Part 2: %s" % str(t2 - t1))
        else:
            background_iu_process = None

    return background_iu_process


@ray.remote(max_calls=1)
def para_iupred(u_ac, seq, iupred_dump):
    ray_hack()

    iupred_path = iupred_dump

    p = subprocess.Popen(['python3', '%s/iupred2a.py' % iupred_path, seq, 'glob'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    out, err = p.communicate()

    after_glob = False
    iupred_parts = [u_ac, [], {}]
    for line in out.split('\n'):
        if len(line) < 3:
            continue
        if line[0] == '#':
            continue
        if line[0] == '1':
            after_glob = True
        words = line.split()
        if not after_glob:
            if len(words) < 4:
                continue
            if words[0] == 'globular' and words[1] == 'domain':
                lower_bound, upper_bound = words[3].split('-')
                iupred_parts[1].append((int(lower_bound), int(upper_bound), 'globular'))
        if after_glob:
            if len(words) < 3:
                continue
            iupred_parts[2][int(words[0])] = words[2]

    return iupred_parts


def paraBlast(config, process_queue_blast, out_queue, lock, i, err_queue):
    blast_path = config.blast_path
    blast_db_path = config.blast_db_path
    option_number_of_templates = config.option_number_of_templates
    option_seq_thresh = config.option_seq_thresh
    option_ral_thresh = config.option_ral_thresh
    option_res_thresh = config.option_res_thresh
    pdb_path = config.pdb_path
    cwd = os.getcwd()
    with lock:
        process_queue_blast.put(None)
    while True:
        inp = process_queue_blast.get()
        if inp is None:
            return
        try:
            (gene, seq) = inp
            target_name = gene.replace("/", "")
            structures = blast.blast(seq, target_name, blast_path, blast_db_path, nr=option_number_of_templates, seq_thresh=option_seq_thresh, cov_thresh=option_ral_thresh, cwd="%s/%d" % (cwd, i))

            structures = templateSelection.selectTemplates(structures, pdb_path)
            structures = templateFiltering.filterTemplates(structures, option_seq_thresh, option_res_thresh, option_ral_thresh)

            with lock:
                out_queue.put((gene, structures))
        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            with lock:
                err_queue.put((e, f, g, gene))


# @profile
def autoTemplateSelection(config, proteins):
    blast_processes = config.blast_processes
    option_seq_thresh = config.option_seq_thresh
    option_res_thresh = config.option_res_thresh
    cwd = os.getcwd()
    pdb_path = config.pdb_path
    errorlog = config.errorlog

    verbose = config.verbose
    search_tool = config.search_tool

    if config.verbosity >= 1:
        print('Sequence search with: ', search_tool)

    if search_tool == 'MMseqs2':
        t0 = time.time()

        raw_structure_map, pdb_ids = MMseqs2.search(proteins, config)

        t1 = time.time()

        info_map = {}
        filtering_results = []
        filtering_dump = ray.put(config.pdb_path)

        small_chunksize, big_chunksize, n_of_small_chunks, n_of_big_chunks = calculate_chunksizes(config.proc_n, len(pdb_ids))

        n = 0
        chunk = []
        for pdb_id in pdb_ids:
            chunk.append(pdb_id)
            if n < n_of_big_chunks:
                if len(chunk) == big_chunksize:
                    filtering_results.append(filter_structures.remote(chunk, filtering_dump))
                    chunk = []
                    n += 1
            elif n <= n_of_big_chunks + n_of_small_chunks:
                if len(chunk) == small_chunksize:
                    filtering_results.append(filter_structures.remote(chunk, filtering_dump))
                    chunk = []
                    n += 1
            else:
                print('Warning: too few chunks')

        if len(chunk) != 0:
            filtering_results.append(filter_structures.remote(chunk, filtering_dump))

        filtering_out = ray.get(filtering_results)

        for out_chunk in filtering_out:
            for (pdb_id, resolution, homomer_dict) in out_chunk:
                info_map[pdb_id] = (resolution, homomer_dict)

        u_acs = proteins.get_protein_ids()
        structure_list = proteins.get_structure_list()
        complex_list = proteins.get_complex_list()

        for u_ac in raw_structure_map:
            for pdb_id, chain in raw_structure_map[u_ac]:
                if pdb_id not in info_map:
                    continue
                resolution, homomer_dict = info_map[pdb_id]
                if resolution is None:
                    continue
                if resolution > option_res_thresh:
                    continue
                oligo = raw_structure_map[u_ac][(pdb_id, chain)]['Oligo']
                struct_anno = sdsc.StructureAnnotation(u_ac, pdb_id, chain)
                proteins.add_annotation(u_ac, pdb_id, chain, struct_anno)

                if not (pdb_id, chain) in structure_list:
                    struct = sdsc.Structure(pdb_id, chain, oligo=oligo, mapped_proteins=[u_ac])
                    proteins.add_structure(pdb_id, chain, struct)
                else:
                    proteins.add_mapping_to_structure(pdb_id, chain, u_ac)

                if pdb_id not in complex_list:
                    compl = sdsc.Complex(pdb_id, resolution, homomers=homomer_dict)
                    proteins.add_complex(pdb_id, compl)

        t2 = time.time()
        if config.verbosity >= 2:
            print("Template Selection Part 1: %s" % (str(t1 - t0)))
            print("Template Selection Part 2: %s" % (str(t2 - t1)))


@ray.remote(max_calls = 1)
def filter_structures(chunk, filtering_dump):
    ray_hack()

    pdb_path = filtering_dump

    outs = []
    for pdb_id in chunk:
        resolution, homomer_dict = pdbParser.getInfo(pdb_id, pdb_path)
        outs.append((pdb_id, resolution, homomer_dict))

    return outs

# @profile


def paraAlignment(config, proteins, skip_db=False, skip_inserts=False, indel_analysis_follow_up=False, get_all_alignments=False):
    indel_analysis_follow_up = indel_analysis_follow_up or get_all_alignments

    alignment_processes = config.alignment_processes
    t0 = time.time()

    if not skip_db:
        database.getAlignments(proteins, config, get_all_alignments=get_all_alignments)

    t1 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 1: %s" % (str(t1 - t0)))

    prot_ids = list(proteins.get_protein_ids())

    mapping_dump = ray.put(config)

    mapping_results = []

    sus_complexes = set()
    sus_structures = set()

    for prot_id in prot_ids:
        if not proteins.is_protein_stored(prot_id):
            continue
        if prot_id not in proteins.indels:
            if proteins.is_completely_stored(prot_id) and not get_all_alignments:
                continue

        annotation_list = proteins.get_protein_annotation_list(prot_id)

        for (pdb_id, chain) in annotation_list:
            if prot_id not in proteins.indels:
                if not proteins.is_annotation_stored(pdb_id, chain, prot_id) and not get_all_alignments:
                    continue
            structure_id = proteins.get_structure_db_id(pdb_id, chain)
            if indel_analysis_follow_up:
                get_alignment_out = proteins.get_alignment(prot_id, pdb_id, chain)
                if isinstance(get_alignment_out, str):
                    config.errorlog.add_error('Get alignment error: %s. %s %s %s %s' % (get_alignment_out, prot_id, pdb_id, chain, str(get_all_alignments)))
                    continue
                else:
                    target_seq, template_seq = get_alignment_out
            else:
                target_seq, template_seq = proteins.pop_alignment(prot_id, pdb_id, chain)
            aaclist = proteins.getAACList(prot_id)
            prot_db_id = proteins.get_protein_db_id(prot_id)
            mapping_results.append(paraMap.remote(mapping_dump, prot_id, pdb_id, chain, structure_id, target_seq, template_seq, aaclist, prot_db_id))

    gc.collect()

    t2 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 2: %s" % (str(t2 - t1)))

    mappings_outs = ray.get(mapping_results)
    for (prot_id, pdb_id, chain, sub_infos, atom_count, last_residue, first_residue) in mappings_outs:
        proteins.set_sub_infos(prot_id, pdb_id, chain, sub_infos)
        proteins.set_atom_count(pdb_id, atom_count)
        proteins.set_last_residue(pdb_id, chain, last_residue)
        proteins.set_first_residue(pdb_id, chain, first_residue)

    t3 = time.time()
    if config.verbosity >= 2:
        print("Alignment Part 3: %s" % (str(t3 - t2)))

    total_cost = 0
    cost_map = {}

    protein_packages = {}

    for prot_id in prot_ids:
        cost_map[prot_id] = 0
        annotation_list = proteins.get_protein_annotation_list(prot_id)
        seq = proteins.get_sequence(prot_id)
        aaclist = proteins.getAACList(prot_id)
        prot_specific_mapping_dump = ray.put((prot_id, seq, aaclist))
        structure_infos = []
        for (pdb_id, chain) in annotation_list:
            if proteins.is_annotation_stored(pdb_id, chain, prot_id):
                continue
            cost_map[prot_id] += (len(seq))**2
            oligo = proteins.get_oligo(pdb_id, chain)
            structure_infos.append((pdb_id, chain, oligo))

        if cost_map[prot_id] == 0:
            print('Warning: alignment cost zero for', prot_id)
            continue

        protein_packages[prot_id] = (prot_specific_mapping_dump, structure_infos)

        total_cost += cost_map[prot_id]

    optimal_chunk_cost = (total_cost// (4 * config.proc_n)) + 1

    if config.verbosity >= 3:
        print('Total cost for alignment:', total_cost, optimal_chunk_cost)

    alignment_results = []

    chunk = []
    chunk_cost = 0

    prots_todo = set(protein_packages.keys())
    started_processes = 0

    while started_processes < config.proc_n and len(prots_todo) > 0:
        done = []
        for prot_id in prots_todo:
            package_cost = cost_map[prot_id]
            if package_cost >= optimal_chunk_cost:  #packages that are greater than the optimal_chunk_cost get split
                structure_infos_split_a = []
                structure_infos_split_b = []
                cost_unit = (len(proteins.get_sequence(prot_id)))**2
                (prot_specific_mapping_dump, structure_infos) = protein_packages[prot_id]
                split_cost = 0
                for structure_info in structure_infos:
                    if split_cost < optimal_chunk_cost:
                        structure_infos_split_a.append(structure_info)
                        split_cost += cost_unit
                    else:
                        structure_infos_split_b.append(structure_info)

                if config.verbosity >= 3:
                    print('Splitting protein alignment:', prot_id, 'split_cost:', split_cost)

                alignment_results.append(align.remote(mapping_dump, [(prot_specific_mapping_dump, structure_infos_split_a)]))
                started_processes += 1
                if len(structure_infos_split_b) > 0:
                    protein_packages[prot_id] = (prot_specific_mapping_dump, structure_infos_split_b)
                    cost_map[prot_id] = package_cost - split_cost
                else:
                    done.append(prot_id)
            elif chunk_cost + package_cost <= optimal_chunk_cost:
                chunk.append(protein_packages[prot_id])
                done.append(prot_id)
                chunk_cost += package_cost
                if chunk_cost >= optimal_chunk_cost:
                    alignment_results.append(align.remote(mapping_dump, chunk))
                    chunk = []
                    chunk_cost = 0
                    started_processes += 1
            if started_processes >= config.proc_n:
                break

        for prot_id in done:
            prots_todo.remove(prot_id)

    if config.verbosity >= 2:
        t4 = time.time()
        print("Alignment Part 4: %s" % (str(t4 - t3)))


    alignment_insertion_list = []
    structure_insertion_list = set()
    warn_map = set()

    safe_complexes = set()
    safe_structures = set()
    database_structure_list = None

    while True:

        if started_processes >= config.proc_n or len(prots_todo) == 0:
            ready, not_ready = ray.wait(alignment_results)
        else:
            ready, not_ready = ray.wait(alignment_results, timeout = 0.01)

        if len(ready) > 0:

            t_i_0 = 0.
            t_i_1 = 0.
            t_i_2 = 0.
            t_i_3 = 0.
            t_i_4 = 0.
            t_i_5 = 0.
            t_i_6 = 0.
            t_i_7 = 0.
            t_i_8 = 0.
            t_i_9 = 0.

            if config.verbosity >= 2:
                t_get_0 = time.time()
            align_outs = ray.get(ready)
            if config.verbosity >= 2:
                t_get_1 = time.time()

            if config.verbosity >= 3:
                print('Times for get:', (t_get_1 - t_get_0))

            for out_chunks in align_outs:
                for out in out_chunks:
                    t_i_0 += time.time()
                    if len(out) == 3:
                        config.errorlog.add_error('Illegal alignment output: %s' % (str(out)))
                        t_i_1 += time.time()
                        continue

                    if len(out) == 4:
                        (prot_id, pdb_id, chain, warn_text) = out
                        proteins.remove_annotation(prot_id, pdb_id, chain)
                        if not skip_inserts:
                            sus_complexes.add(pdb_id)
                            sus_structures.add((pdb_id, chain))
                        if prot_id not in warn_map:
                            config.errorlog.add_warning(warn_text)
                        elif config.verbosity >= 3:
                            print('Alignment failed:', prot_id, pdb_id, chain, warn_text)
                        warn_map.add(prot_id)
                        t_i_1 += time.time()
                        continue

                    if len(out) == 5:
                        (prot_id, pdb_id, chain, sub_infos, seq_id) = out
                        proteins.remove_annotation(prot_id, pdb_id, chain)
                        if not skip_inserts:
                            sus_complexes.add(pdb_id)
                            sus_structures.add((pdb_id, chain))
                        if config.verbosity >= 4:
                            print('Alignment got filtered:', prot_id, pdb_id, chain, len(sub_infos), seq_id)
                        t_i_1 += time.time()
                        continue
                    t_i_1 += time.time()
                    t_i_2 += time.time()

                    (prot_id, pdb_id, chain, alignment, seq_id, coverage, interaction_partners, chain_type_map,
                     oligo, sub_infos, atom_count, last_residue, first_residue, chainlist, rare_residues) = out

                    proteins.set_coverage(prot_id, pdb_id, chain, coverage)
                    proteins.set_sequence_id(prot_id, pdb_id, chain, seq_id)
                    proteins.set_sub_infos(prot_id, pdb_id, chain, sub_infos)
                    proteins.set_atom_count(pdb_id, atom_count)
                    if indel_analysis_follow_up:
                        proteins.set_alignment(prot_id, pdb_id, chain, alignment)

                    t_i_3 += time.time()

                    proteins.set_interaction_partners(pdb_id, interaction_partners)
                    proteins.set_chain_type_map(pdb_id, chain_type_map, chainlist)

                    t_i_4 += time.time()

                    proteins.set_oligo(pdb_id, chain, oligo)
                    proteins.set_last_residue(pdb_id, chain, last_residue)
                    proteins.set_first_residue(pdb_id, chain, first_residue)

                    t_i_5 += time.time()

                    if (pdb_id, chain) not in safe_structures:
                        structure_insertion_list.add((pdb_id, chain))

                    t_i_6 += time.time()

                    safe_complexes.add(pdb_id)
                    safe_structures.add((pdb_id, chain))

                    t_i_7 += time.time()

                    prot_db_id = proteins.get_protein_db_id(prot_id)
                    alignment_insertion_list.append((prot_id, prot_db_id, pdb_id, chain, alignment))

                    t_i_8 += time.time()

                    config.rare_residues.update(rare_residues)

                    t_i_9 += time.time()

                started_processes -= 1

            if config.verbosity >= 3:
                print('Times for data integration:', (t_i_1 - t_i_0), (t_i_3 - t_i_2), (t_i_4 - t_i_3), (t_i_5 - t_i_4), (t_i_6 - t_i_5), (t_i_7 - t_i_6), (t_i_8 - t_i_7), (t_i_9 - t_i_8))

        t_new_0 = time.time()

        new_alignment_results = []
        if len(prots_todo) > 0:
            done = []
            new_alignment_results = []
            for prot_id in prots_todo:
                package_cost = cost_map[prot_id]

                if package_cost >= optimal_chunk_cost:  #packages that are greater than the optimal_chunk_cost get split
                    structure_infos_split_a = []
                    structure_infos_split_b = []
                    cost_unit = (len(proteins.get_sequence(prot_id)))**2
                    (prot_specific_mapping_dump, structure_infos) = protein_packages[prot_id]
                    split_cost = 0
                    for structure_info in structure_infos:
                        if split_cost < optimal_chunk_cost:
                            structure_infos_split_a.append(structure_info)
                            split_cost += cost_unit
                        else:
                            structure_infos_split_b.append(structure_info)

                    if config.verbosity >= 3:
                        print('Splitting protein alignment:', prot_id, 'split_cost:', split_cost)


                    new_alignment_results.append(align.remote(mapping_dump, [(prot_specific_mapping_dump, structure_infos_split_a)]))
                    started_processes += 1
                    if len(structure_infos_split_b) > 0:
                        protein_packages[prot_id] = (prot_specific_mapping_dump, structure_infos_split_b)
                        cost_map[prot_id] = package_cost - split_cost
                    else:
                        done.append(prot_id)

                else:
                    overload = min([optimal_chunk_cost - chunk_cost, optimal_chunk_cost//2])
                    if chunk_cost + package_cost <= optimal_chunk_cost + overload: #chunks may have greater costs than optimal_chunk_cost
                        chunk.append(protein_packages[prot_id])
                        done.append(prot_id)
                        chunk_cost += package_cost
                        if chunk_cost >= optimal_chunk_cost:
                            new_alignment_results.append(align.remote(mapping_dump, chunk))
                            chunk = []
                            chunk_cost = 0
                            started_processes += 1
                        
                if started_processes >= config.proc_n:
                    break

            for prot_id in done:
                prots_todo.remove(prot_id)

            if len(chunk) > 0 and started_processes < config.proc_n:
                new_alignment_results.append(align.remote(mapping_dump, chunk))
                chunk = []
                chunk_cost = 0
                started_processes += 1
            if len(prots_todo) > 0 and started_processes < config.proc_n:
                prot_id = prots_todo.pop()
                new_alignment_results.append(align.remote(mapping_dump, [protein_packages[prot_id]]))
                started_processes += 1

        alignment_results = new_alignment_results + not_ready

        t_new_1 = time.time()
        if config.verbosity >= 3:
            print('Times for starting new processes:', (t_new_1 - t_new_0))
            print('Proteins to do:', len(prots_todo), 'Chunk:', len(chunk), 'Pending processes:', len(alignment_results))

        if len(alignment_results) == 0 and len(prots_todo) == 0 and len(chunk) == 0:
            if started_processes != 0 or len(prots_todo) != 0:
                print('\nWarning: left alignment loop while stuff was still to do', started_processes, prots_todo,'\n')
            break


    if config.verbosity >= 2:
        t5 = time.time()
        print("Alignment Part 5: %s" % (str(t5 - t4)))
        if len(config.rare_residues) > 0:
            print('######### Detected some rare residues:', config.rare_residues)

    if not skip_db and not skip_inserts:
        database_structure_list = database.insertStructures(structure_insertion_list,
                                                            proteins, config, results=database_structure_list,
                                                            return_results=config.low_mem_system)

    if config.verbosity >= 2:
        t6 = time.time()
        print("Alignment Part 6: %s" % (str(t6 - t5)))

    if not skip_db and not skip_inserts:
        if config.verbosity >= 2:
            t7 = time.time()
            print("Alignment Part 7: %s" % (str(t7 - t6)))

        database.insertAlignments(alignment_insertion_list, proteins, config)

        if config.verbosity >= 2:
            t8 = time.time()
            print("Alignment Part 8: %s" % (str(t8 - t7)))

    # Due the removal of annotations in the previous loop, we might to remove some structures and complexes
    proteins.remove_structures(sus_structures - safe_structures)
    complexes_to_remove = sus_complexes - safe_complexes
    if config.verbosity >= 3:
        print('Len of sus_complexes:', len(sus_complexes), 'Len of safe complexes:', len(safe_complexes), 'Len of complexes_to_remove:', len(complexes_to_remove))
        if len(complexes_to_remove) < 50:
            print('Remove complexes:', complexes_to_remove)
    proteins.remove_complexes(complexes_to_remove)

    if skip_db:
        # even lite mode checks for stored structures
        database.structureCheck(proteins, config)


@ray.remote(max_calls = 1)
def paraMap(mapping_dump, u_ac, pdb_id, chain, structure_id, target_seq, template_seq, aaclist, prot_id):
    ray_hack()

    config = mapping_dump

    pdb_path = config.pdb_path

    template_page, atom_count = pdbParser.standardParsePDB(pdb_id, pdb_path, obsolete_check=True)

    seq_res_map, last_residue, first_residue = globalAlignment.createTemplateFasta(template_page, pdb_id, chain, config, onlySeqResMap=True)

    sub_out = globalAlignment.getSubPos(config, u_ac, target_seq, template_seq, aaclist, seq_res_map)

    if isinstance(sub_out, str):
        config.errorlog.add_error("%s %s %s\n%s\n%s" % (sub_out, pdb_id, chain, template_seq.replace('-', ''), seq_res_map))

    sub_infos, aaclist = sub_out

    return (u_ac, pdb_id, chain, sub_infos, atom_count, last_residue, first_residue)


@ray.remote(max_calls = 1)
def align(align_dump, package, model_path=None):
    ray_hack()

    config = align_dump

    results = []

    for prot_specific_mapping_dump, structure_infos in package:

        (u_ac, seq, aaclist) = ray.get(prot_specific_mapping_dump, timeout = 600)
        for pdb_id, chain, oligo in structure_infos:

            pdb_path = config.pdb_path
            option_seq_thresh = config.option_seq_thresh

            try:
                parse_out = pdbParser.getStandardizedPdbFile(pdb_id, pdb_path, oligo=oligo, verbosity=config.verbosity, model_path=model_path)

                if parse_out is None:
                    results.append((u_ac, pdb_id, chain, 'pdbParser failed for %s' % str(pdb_id)))
                    continue

                (template_page, interaction_partners, chain_type_map, oligo, atom_count, chainlist, rare_residues) = parse_out

                align_out = globalAlignment.alignBioPython(config, u_ac, seq, pdb_id, template_page, chain, aaclist, rare_residues=rare_residues)

                if isinstance(align_out, str):
                    align_out = '%s %s' % (align_out, str(model_path))
                    results.append((u_ac, pdb_id, chain, align_out))
                    continue

                (coverage, seq_id, sub_infos, alignment_pir, times, aaclist, last_residue, first_residue) = align_out

                if sub_infos is None or seq_id is None:
                    results.append((u_ac, pdb_id, chain, sub_infos, seq_id))
                    continue

                if 100.0 * seq_id >= option_seq_thresh:
                    results.append((u_ac, pdb_id, chain, alignment_pir, seq_id, coverage, interaction_partners, chain_type_map, oligo, sub_infos, atom_count, last_residue, first_residue, chainlist, rare_residues))
                else:
                    results.append((u_ac, pdb_id, chain, sub_infos, seq_id))

            except:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
                results.append((u_ac, pdb_id, chain, '%s,%s\n%s\n%s\n%s' % (pdb_id, chain, e, str(f), g)))

    return results


@ray.remote(max_calls = 1)
def para_classify_remote_wrapper(classification_dump, package, package_size):
    ray_hack()

    return (para_classify(classification_dump, unpack(package), para=True), package_size)


def para_classify(classification_dump, package, para=False):
    t0 = time.time()
    if para:
        config, complexes = classification_dump
    else:
        config = classification_dump

    outs = []
    for u_ac, classification_inp in package:
        for pos, mappings, disorder_score, disorder_region in classification_inp:

            mappings_obj = sdsc.Mappings()

            for mapp in mappings:

                if para:
                    id_triple, quality_measures, mapping = mapp

                    (pdb_id, chain, res_nr) = id_triple
                    (seq_id, cov, identical_aa) = quality_measures

                    (chains, resolution) = complexes[pdb_id]
                else:
                    id_triple, quality_measures, structure_infos = mapp

                    (pdb_id, chain, res_nr) = id_triple
                    (seq_id, cov, identical_aa) = quality_measures

                    structure_info, mapping = structure_infos
                    (chains, resolution) = structure_info

                (rsa, mc_rsa, sc_rsa, ssa, profile_or_str, centralities_or_str,
                 phi, psi, intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
                 inter_chain_median_kd, inter_chain_dist_weighted_kd,
                 inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                 intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                 inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                 intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                 b_factor, modres,
                 lig_dists, chain_distances, homomer_distances, res_aa) = mapping

                gsd_return = sdsc.get_shortest_distances(chains, lig_dists, chain_distances, homomer_distances)

                if gsd_return is None:
                    continue
                homo_dist, lig_dist, metal_dist, ion_dist, chain_dist, rna_dist, dna_dist, min_lig, min_metal, min_ion, iacs = gsd_return

                if isinstance(profile_or_str, str):
                    profile = rin.Interaction_profile(profile_str=profile_or_str)
                else:
                    profile = profile_or_str

                if isinstance(centralities_or_str, str):
                    centralities = rin.Centrality_scores(code_str=centralities_or_str)
                else:
                    centralities = centralities_or_str

                loc, mc_loc, sc_loc = sdsc.triple_locate(rsa, mc_rsa, sc_rsa, config)

                rin_class, rin_simple_class = sdsc.rin_classify(profile, sc_loc)

                qual = templateFiltering.qualityScore(resolution, cov, seq_id)

                mapping = (qual, seq_id, cov, rsa, mc_rsa, sc_rsa, ssa, lig_dist, metal_dist, ion_dist, chain_dist,
                           rna_dist, dna_dist, homo_dist, profile, centralities,
                           phi, psi, intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
                           inter_chain_median_kd, inter_chain_dist_weighted_kd,
                           inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                           intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                           inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                           intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                           b_factor, modres, rin_class, rin_simple_class,
                           identical_aa, resolution, res_aa)

                mappings_obj.add_mapping((pdb_id, chain, res_nr), mapping)

            mappings_obj.weight_all(config, disorder_score, disorder_region)

            mapping_results = mappings_obj.get_raw_result()

            outs.append((u_ac, pos, mapping_results))

    if para:
        outs = pack(outs)

    t1 = time.time()

    return outs, t1 - t0


def get_res_info_from_store(structure, res_nr):
    res_aa = structure.residues[res_nr].aa

    centralities_or_str = structure.residues[res_nr].get_centralities(get_whats_there=True)
    modres = structure.residues[res_nr].modres
    b_factor = structure.residues[res_nr].b_factor
    rsa, mc_rsa, sc_rsa = structure.get_residue_rsa_triple(res_nr)
    ssa = structure.residues[res_nr].SSA
    profile_or_str = structure.residues[res_nr].get_interaction_profile(get_whats_there=True)

    phi = structure.residues[res_nr].phi
    psi = structure.residues[res_nr].psi

    intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower = structure.get_residue_link_information(res_nr)

    (inter_chain_median_kd, inter_chain_dist_weighted_kd,
     inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
     intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa) = structure.get_residue_milieu(res_nr)

    (inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
     intra_chain_interactions_median, intra_chain_interactions_dist_weighted) = structure.residues[res_nr].get_interface_milieu()

    lig_dists = structure.get_residue_sld(res_nr)
    chain_distances = structure.get_residue_scd(res_nr)
    homomer_distances = structure.get_residue_homomer_dists(res_nr)

    res_info = (rsa, mc_rsa, sc_rsa, ssa, profile_or_str, centralities_or_str,
                phi, psi, intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
                inter_chain_median_kd, inter_chain_dist_weighted_kd,
                inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                b_factor, modres, lig_dists, chain_distances, homomer_distances, res_aa)
    return res_info


def get_res_info(proteins, pdb_id, chain, res_nr):

    residue_obj = proteins.structures[(pdb_id, chain)].residues[res_nr]

    centralities_or_str = residue_obj.get_centrality_str()

    rsa, mc_rsa, sc_rsa = residue_obj.get_rsa(splitted = True)
    profile_or_str = residue_obj.get_interaction_profile_str()

    intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower = residue_obj.get_residue_link_information()

    (inter_chain_median_kd, inter_chain_dist_weighted_kd,
     inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
     intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa) = residue_obj.get_milieu()

    (inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
     intra_chain_interactions_median, intra_chain_interactions_dist_weighted) = residue_obj.get_interface_milieu()

    lig_dists = residue_obj.get_ligand_distances()
    chain_distances = residue_obj.get_chain_distances()
    homomer_distances = residue_obj.get_homomer_dists()

    res_info = (rsa, mc_rsa, sc_rsa, residue_obj.SSA, profile_or_str, centralities_or_str,
                residue_obj.phi, residue_obj.psi, intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
                inter_chain_median_kd, inter_chain_dist_weighted_kd,
                inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                residue_obj.b_factor, residue_obj.modres, lig_dists, chain_distances, homomer_distances, residue_obj.aa)
    return res_info


def add_structure_to_store(structure_store, proteins, pdb_id, chain):
    chains = proteins.get_complex_chains(pdb_id)
    resolution = proteins.get_resolution(pdb_id)
    if resolution is None:
        if config.verbosity >= 3:
            print('Skipped classification of', u_ac, 'due to', pdb_id, 'got no resolution')
        return
    res_map = {}
    for res_nr in proteins.structures[(pdb_id, chain)].residues:
        res_map[res_nr] = get_res_info(proteins, pdb_id, chain, res_nr)

    structure_store[(pdb_id, chain)] = ray.put((chains, resolution, res_map))


def classification(proteins, config, indel_analysis_follow_up=False, custom_structure_annotations=None):

    def pack_packages(package_size, package, send_packages, total_pending_package_size, classification_results, proteins, config, size_sorted, processed_positions):
        if config.verbosity >= 3:
            print('Call of pack_packages with:', package_size, len(package), send_packages)
            print('Current CPU load:', psutil.cpu_percent())
        t_pack_0 = 0
        t_pack_01 = 0
        t_pack_02 = 0
        t_pack_03 = 0
        t_pack_1 = 0
        t_pack_11 = 0
        t_pack_12 = 0
        t_pack_13 = 0
        t_pack_14 = 0
        t_pack_15 = 0
        t_pack_16 = 0
        t_pack_21 = 0
        t_pack_22 = 0

        snap_size_sorted = [x for x in size_sorted]
        for u_ac, size in snap_size_sorted:
            classification_inp = []
            positions = proteins.get_position_ids(u_ac)
            annotation_list = proteins.get_protein_annotation_list(u_ac)

            if not u_ac in processed_positions:
                processed_positions[u_ac] = set()

            protein_obj = proteins.protein_map[u_ac]

            for pos in positions:
                if pos in processed_positions[u_ac]:
                    continue
                if protein_obj.is_position_stored(pos):
                    continue
                t_pack_0 += time.time()

                mappings = []
                aacbase = protein_obj.get_aac_base(pos)
                disorder_score = protein_obj.get_disorder_score(pos)
                disorder_region = protein_obj.get_disorder_region(pos)

                t_pack_01 += time.time()

                for (pdb_id, chain) in annotation_list:

                    t_pack_11 += time.time()

                    try:
                        sub_info = protein_obj.structure_annotations[(pdb_id, chain)].get_sub_info(pos)

                        t_pack_12 += time.time()
                    except:
                        if config.verbosity >= 3:
                            print('Skipped classification of', u_ac, 'due to pos', pos, 'was not in sub_infos (len:', len(sub_infos), ')')

                        t_pack_12 += time.time()
                        continue


                    res_nr = sub_info[0]

                    if res_nr is None:
                        if config.verbosity >= 4:
                            print('Skipped classification of', u_ac, pos, 'due to res_nr is None in ', pdb_id, chain)
                        continue

                    t_pack_13 += time.time()

                    seq_id = protein_obj.structure_annotations[(pdb_id, chain)].get_sequence_id()

                    t_pack_14 += time.time()

                    if seq_id is None:
                        if config.verbosity >= 3:
                            print('Skipped classification of', u_ac, 'due to', pdb_id, 'got no sequence identity')
                        continue

                    t_pack_15 += time.time()

                    cov = protein_obj.structure_annotations[(pdb_id, chain)].get_coverage()

                    t_pack_16 += time.time()

                    try:
                        res_aa = proteins.get_residue_aa(pdb_id, chain, res_nr)
                    except:
                        if config.verbosity >= 3:
                            print('Skipped classification of', u_ac, pos, 'due to res_nr', res_nr, 'was not contained in ', pdb_id, chain)
                        continue
                    identical_aa = res_aa == aacbase[0]

                    t_pack_21 += time.time()

                    if not para:
                        resolution = proteins.get_resolution(pdb_id)
                        if resolution is None:
                            if config.verbosity >= 3:
                                print('Skipped classification of', u_ac, 'due to', pdb_id, 'got no resolution')
                            continue
                        package_size += 1
                        chains = proteins.get_complex_chains(pdb_id)
                        mappings.append(((pdb_id, chain, res_nr), (seq_id, cov, identical_aa), ((chains, resolution), get_res_info(proteins, pdb_id, chain, res_nr))))
                    else:
                        package_size += 1
                        mappings.append(((pdb_id, chain, res_nr), (seq_id, cov, identical_aa), get_res_info(proteins, pdb_id, chain, res_nr)))

                    t_pack_22 += time.time()

                t_pack_02 += time.time()

                classification_inp.append((pos, mappings, disorder_score, disorder_region))

                t_pack_03 += time.time()

                processed_positions[u_ac].add(pos)

                t_pack_1 += time.time()

                if package_size >= max_package_size and para:

                    package.append((u_ac, classification_inp))

                    t_send_0 = time.time()
                    classification_results.append(para_classify_remote_wrapper.remote(classification_dump, pack(package), package_size))
                    t_send_1 = time.time()

                    total_pending_package_size += package_size

                    send_packages += 1
                    if config.verbosity >= 3:
                        print('Start remote classifcation', u_ac, 'Package size:', package_size,'Amount of different proteins in the package:', len(package))
                        print('Amount of positions in last protein:', len(classification_inp))
                        print('Pending processes:', send_packages, 'Size of all packages in all pending processes combined:', total_pending_package_size)
                        print('RAM memory % used:', psutil.virtual_memory()[2])
                        print('Time for sending the process:', (t_send_1 - t_send_0), 'Time for packaging:', (t_pack_1 - t_pack_0), (t_pack_01 - t_pack_0), (t_pack_02 - t_pack_01), (t_pack_03 - t_pack_02), (t_pack_1 - t_pack_03), (t_pack_22 - t_pack_21))
                        print('More times:', (t_pack_12 - t_pack_11), (t_pack_14 - t_pack_13), (t_pack_16 - t_pack_15), (t_pack_21 - t_pack_16))
                        t_pack_0 = 0
                        t_pack_01 = 0
                        t_pack_02 = 0
                        t_pack_03 = 0
                        t_pack_1 = 0
                        t_pack_11 = 0
                        t_pack_12 = 0
                        t_pack_13 = 0
                        t_pack_14 = 0
                        t_pack_15 = 0
                        t_pack_16 = 0
                        t_pack_21 = 0
                        t_pack_22 = 0
                    package_size = 0
                    package = []
                    classification_inp = []

                    if send_packages >= config.proc_n:
                        return package_size, package, send_packages, total_pending_package_size, classification_results, size_sorted, processed_positions
            package.append((u_ac, classification_inp))
            del size_sorted[0]
            del processed_positions[u_ac]
        if len(package) > 0 and para:
            classification_results.append(para_classify_remote_wrapper.remote(classification_dump, pack(package), package_size))
            send_packages += 1
            package = []
            total_pending_package_size += package_size
            package_size = 0

        return package_size, package, send_packages, total_pending_package_size, classification_results, size_sorted, processed_positions

    if config.verbosity >= 2:
        t0 = time.time()

    classification_results = []

    size_sorted = []

    if custom_structure_annotations is None:
        u_acs = proteins.get_protein_ids()
    else:
        u_acs = custom_structure_annotations.keys()

    total_mappings = 0
    for u_ac in u_acs:
        if proteins.is_completely_stored(u_ac):
            continue
        if custom_structure_annotations is None:
            annotation_list = proteins.get_protein_annotation_list(u_ac)
        else:
            annotation_list = custom_structure_annotations[u_ac]
        size = len(annotation_list) * len(proteins[u_ac].positions)
        size_sorted.append((u_ac, size))

        total_mappings += size

    if config.verbosity >= 3:
        print('Total mappings:', total_mappings)

    size_sorted = sorted(size_sorted, reverse=True, key=lambda x: x[1])

    packages = []
    package_counter = 0

    para = (total_mappings > 20000)

    max_package_size = min([50000, total_mappings // (config.proc_n * 4)])
    package_size = 0
    package = []

    send_packages = 0
    total_pending_package_size = 0

    if para:
        complex_store = {}

        t_init_0 = time.time()

        for u_ac, size in size_sorted:
            annotation_list = proteins.get_protein_annotation_list(u_ac)

            for (pdb_id, chain) in annotation_list:
                if not pdb_id in complex_store:
                    complex_store[pdb_id] = (proteins.complexes[pdb_id].chains, proteins.complexes[pdb_id].resolution)

        t_init_1 = time.time()
        if config.verbosity >= 2:
            print('Time for para classification preprocessing 1:', (t_init_1 - t_init_0))

        classification_dump = ray.put((config, complex_store))

        t_init_2 = time.time()
        if config.verbosity >= 2:
            print('Time for para classification preprocessing 2:', (t_init_2 - t_init_1))

    processed_positions = {}
    package_size, package, send_packages, total_pending_package_size, classification_results, size_sorted, processed_positions = pack_packages(package_size, package, send_packages, total_pending_package_size, classification_results, proteins, config, size_sorted, processed_positions)

    # para = send_packages > 0  # only use ray, when we have more than one package

    if not para and len(package) > 0:
        classification_results.append(para_classify(config, package))
    elif not para and len(package) == 0:
        return

    if config.verbosity >= 2:
        t11 = time.time()
        print('Time for classification part 1.1:', t11 - t0)

    if config.verbosity >= 2:
        t1 = time.time()
        print('Time for classification part 1:', t1 - t0, para, max_package_size)

    if para:
        max_comp_time = 0
        total_comp_time = 0.
        max_comp_time_pos = None
        amount_of_positions = 0

        total_processed_mappings = 0

        while True:
            #if config.low_mem_system:
            #t_gc_0 = time.time()
            #gc.collect()
            #t_gc_1 = time.time()

            #if config.verbosity >= 3:
            #    print('Time for garbage collection:', (t_gc_1 - t_gc_0))

            ready, not_ready = ray.wait(classification_results)

            if len(ready) > 0:

                t_integrate_0 = 0.
                t_integrate_1 = 0.
                t_integrate_2 = 0.
                send_packages -= len(ready)

                t_get_0 = time.time()

                remote_returns = ray.get(ready)

                t_get_1 = time.time()

                for returned_package, size_of_returned_package in remote_returns:
                    (outs, comp_time) = returned_package
                    total_pending_package_size -= size_of_returned_package
                    total_processed_mappings += size_of_returned_package
                    t_integrate_0 += time.time()
                    outs = unpack(outs)
                    t_integrate_1 += time.time()
                    for u_ac, pos, mapping_results in outs:
                        proteins.protein_map[u_ac].positions[pos].mappings = sdsc.Mappings(raw_results=mapping_results)

                    total_comp_time += comp_time
                    amount_of_positions += 1

                    if comp_time > max_comp_time:
                        max_comp_time = comp_time
                        max_comp_time_pos = u_ac
                    t_integrate_2 += time.time()

                classification_results = not_ready

                t_pack_0 = time.time()
                package_size, package, send_packages, total_pending_package_size, classification_results, size_sorted, processed_positions = pack_packages(package_size, package, send_packages, total_pending_package_size, classification_results, proteins, config, size_sorted, processed_positions)
                t_pack_1 = time.time()

                if config.verbosity >= 3:
                    print('Time for getting:', (t_get_1 - t_get_0))
                    print('Time for unpacking:', (t_integrate_1 - t_integrate_0))
                    print('Time for result integration:', (t_integrate_2 - t_integrate_1))
                    print('Time for packaging new processes:', (t_pack_1 - t_pack_0))
                    print('Progress:', total_processed_mappings, 'of', total_mappings)
            else:
                classification_results = not_ready

            if len(classification_results) == 0:
                break
    else:
        gc.collect()
        for outs, comp_time in classification_results:
            for u_ac, pos, mapping_results in outs:
                proteins.protein_map[u_ac].positions[pos].mappings = sdsc.Mappings(raw_results=mapping_results)

    if not indel_analysis_follow_up:
        for u_ac in u_acs:
            proteins.remove_protein_annotations(u_ac)
        proteins.semi_deconstruct()


    if config.verbosity >= 2:
        t2 = time.time()
        print('Time for classification part 2:', t2 - t1)
        if para:
            print('Longest computation for:', max_comp_time_pos, 'with:', max_comp_time, 'In total', amount_of_positions, 'proteins', 'Accumulated time:', total_comp_time)


def core(protein_list, indels, multi_mutation_objects, config, session, outfolder, session_name, out_objects):
    if len(protein_list) == 0:
        return None, 0

    ray_init(config)

    #monitor_ray_store()

    background_insert_residues_process = None
    background_process_MS = None

    # transform the protein map into a Proteins object
    proteins = sdsc.Proteins(protein_list, indels, multi_mutation_objects, lite=config.lite)

    if config.verbosity >= 4:
        print('Proteins state after object initialization:')
        proteins.print_protein_state()

    try:
        t1 = time.time()
        if config.verbosity >= 1:
            print("Before protCheck")
        # check for already stored genes and make all the necessary database interactions
        # sets the fields database_id and stored in the protein objects as well as the stored and non stored ids in the proteins object

        if not config.lite:
            database.protCheck(proteins, session, config)

        if not config.only_wt:
            for wt_prot_id in proteins.multi_mutations:
                if config.verbosity >= 3:
                    print('Mutating', len(proteins.multi_mutations[wt_prot_id]), 'multi mutations for:', wt_prot_id)
                for multi_mutation in proteins.multi_mutations[wt_prot_id]:
                    multi_mutation.mutate(proteins, config)

        if config.verbosity >= 4:
            print('Proteins state after protCheck:')
            proteins.print_protein_state()

        t2 = time.time()
        if config.verbosity >= 2:
            print("Time for protCheck: %s" % (str(t2 - t1)))
        # check for already stored mutations or position twins and make all the necessary database interactions
        # sets the fields database_id and stored in the position objects

        if not config.lite:
            background_process_MS = database.positionCheck(proteins, session, config)

            database.indelCheck(proteins, session, config)

            database.insertMultiMutations(proteins, session, config)

        if config.verbosity >= 4:
            print('Proteins state after positionCheck:')
            proteins.print_protein_state()

        t3 = time.time()
        if config.verbosity >= 2:
            print("Time for positionCheck: %s" % (str(t3 - t2)))

        if config.verbosity >= 1:
            print("Before getSequences")

        background_iu_process = getSequences(proteins, config, skip_db=config.lite)

        t4 = time.time()
        if config.verbosity >= 2:
            print("Time for getSequences: %s" % (str(t4 - t3)))

        if config.verbosity >= 1:
            print("Before autoTemplateSelection")
        autoTemplateSelection(config, proteins)

        t5 = time.time()
        if config.verbosity >= 2:
            print("Time for Template Selection: %s" % (str(t5 - t4)))

        indel_analysis_follow_up = len(proteins.indels) > 0 and not config.skip_indel_analysis

        if config.verbosity >= 1:
            print("Before paraAlignment")
        paraAlignment(config, proteins, skip_db=config.lite, indel_analysis_follow_up=indel_analysis_follow_up)

        t6 = time.time()
        if config.verbosity >= 2:
            print("Time for Alignment: %s" % (str(t6 - t5)))

        if background_iu_process is not None:
            background_iu_process.join()  # Disorder values have to finished before the classification happens

        if config.verbosity >= 1:
            print("Before paraAnnotate")

        background_insert_residues_process, amount_of_structures = templateFiltering.paraAnnotate(config, proteins, indel_analysis_follow_up=indel_analysis_follow_up, lite=config.lite)

        t7 = time.time()
        if config.verbosity >= 2:
            print("Time for Annotation: %s" % (str(t7 - t6)))

        if not config.lite and not config.skip_indel_analysis:
            indel_analysis.para_indel_analysis(proteins, config)
        else:
            n_indel_types = None
            results_tuples_not_final = None
            out_objects = output.appendOutput(proteins, outfolder, session_name, out_objects=out_objects)

        t1 = time.time()
        if config.verbosity >= 2:
            print("Time for Indelanalysis: %s" % (str(t1 - t7)))

        # join the background inserts
        if background_process_MS is not None:
            background_process_MS.join()
        if background_insert_residues_process is not None:
            background_insert_residues_process.join()

        t2 = time.time()
        if config.verbosity >= 2:
            print('Resttime for background inserts: ', t2 - t1)

        if config.verbosity >= 3:
            for name, size in sorted(((name, sys.getsizeof(value)) for name, value in locals().items()), key=lambda x: -x[1])[:10]:
                print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))
            print(list_fds())

    # Error-Handling for a whole input line
    except:

        [e, f, g] = sys.exc_info()
        g = traceback.format_exc()
        # print "Pipeline Core Error: ",e,f,g
        errortext = '\n'.join([str(e), str(f), str(g)]) + '\n\n'
        config.errorlog.add_error(errortext)
        if background_process_MS is not None:
            background_process_MS.join()
        if background_insert_residues_process is not None:
            background_insert_residues_process.join()
        amount_of_structures = 0

    ray.shutdown()

    return out_objects, amount_of_structures


# @profile
def main(filename, config):
    n_of_cores = config.proc_n
    mrna_fasta = config.mrna_fasta
    num_of_cores = config.proc_n
    verbose = config.verbose
    session = 0  # This can later be used, structman_main.py could give a specific session id and the pipeline can then expand that session

    if config.verbosity >= 2:
        print('==== Start main pipeline ====')

    config.errorlog.start(filename, config.outfolder)

    t0 = time.time()

    # need structman package path for ray
    ray_init(config)

    if config.verbosity >= 2:
        print('ray init successful')

    '''
    mem_tracked_processes = []
    for proc in psutil.process_iter(['pid', 'name', 'username']):
        if proc.info['username'] == 'agr18':
            p_id = proc.info['pid']
            mem_tracked_processes.append(psutil.Process())
    '''

    if mrna_fasta is not None:
        if not os.path.exists(mrna_fasta):
            raise NameError("mRNA path not found: %s" % mrna_fasta)

    if isinstance(filename, str):
        # annovar-pipeline in case of vcf-file
        if filename.rsplit(".", 1)[1] == "vcf":
            anno_db = "%s_annovar" % db_name.rsplit("_", 1)[0]  # BUG: undefined variable
            if config.verbosity >= 1:
                print('Convert vcf file format using Annovar')
            if mrna_fasta is not None:
                '... and using mrna file: ', mrna_fasta
            nfname = annovar.annovar_pipeline(filename, config.tax_id, config.annovar_path, config.db_address, config.db_user_name, config.db_password, anno_db, mrna_fasta, ref_id=config.ref_genome_id)
        else:
            nfname = filename

    # Single line input check
    else:
        nfname = 'Single line input'
        single_line_inputs = filename
        if config.verbosity >= 1:
            print('=== Single line input mode ===')

    t01 = time.time()

    if config.verbosity >= 2:
        print("Time for preparation before buildQueue: %s" % (str(t01 - t0)))

    try:
        os.stat("%s/tmp_structman_pipeline" % (config.outfolder))
    except:
        os.mkdir("%s/tmp_structman_pipeline" % (config.outfolder))
    os.chdir("%s/tmp_structman_pipeline" % (config.outfolder))
    cwd = "%s/tmp_structman_pipeline" % (config.outfolder)

    config.temp_folder = cwd

    chunksize = config.chunksize

    if config.verbosity >= 1:
        print("Call buildQueue with chunksize: %s and file: %s" % (str(chunksize), nfname))

    if nfname != 'Single line input':
        if config.fasta_input:
            proteins_chunks = parseFasta(config, nfname)
            temp_infiles = [None]
        else:
            proteins_chunks, temp_infiles = buildQueue(config, nfname)
    else:
        proteins_chunks, temp_infiles = buildQueue(config, single_line_inputs)
        nfname = '/%s.' % (' '.join(single_line_inputs))

    t02 = time.time()
    if config.verbosity >= 2:
        print("Time for buildQueue: %s" % (str(t02 - t01)))
    if config.verbosity >= 1:
        print("Number of chunks: ", len(proteins_chunks))

    newsession = False
    if session == 0 and not config.lite:
        starttime = SQLDateTime()
        session = database.insertSession(starttime, nfname, config)
        newsession = True
    session_name = (nfname.rsplit("/", 1)[1]).rsplit(".", 1)[0]

    out_objects = None

    total_amount_of_analyzed_structures = 0

    for nr_temp_file, temp_infile in enumerate(temp_infiles):
        if temp_infile is not None:
            if config.verbosity >= 1:
                print('Infile splitting due to low memory system, processing infile split nr.:', nr_temp_file + 1, 'out of', len(temp_infiles))
            proteins_chunks, nothing = buildQueue(config, temp_infile, already_split=True)
            os.remove(temp_infile)

        chunk_nr = 1
        for protein_list, indels in proteins_chunks:

            config.indels_given_by_input = (len(indels) > 0)

            if config.verbosity >= 1:
                print("Chunk %s/%s" % (str(chunk_nr), str(len(proteins_chunks))))
            chunk_nr += 1

            protein_list, indels, multi_mutation_objects = sequenceScan(config, protein_list, indels)

            out_objects, amount_of_structures = core(protein_list, indels, multi_mutation_objects, config, session, config.outfolder, session_name, out_objects)

            total_amount_of_analyzed_structures += amount_of_structures

    os.chdir(config.outfolder)

    if total_amount_of_analyzed_structures == 0:
        if config.verbosity >= 1:
            print("=== Nothing got processed, probably the input didn't contain any usable protein ID. ===")

    if config.verbosity >= 2:
        t03 = time.time()

    try:
        shutil.rmtree(cwd)
    except:
        pass

    if config.verbosity >= 2:
        t04 = time.time()
        print('Time for folder cleanup:', t04 - t03)

    ray.shutdown()

    config.errorlog.stop()

    if newsession:
        endtime = SQLDateTime()
        database.updateSession(session, endtime, config)

    tend = time.time()
    if config.verbosity >= 1:
        print(total_amount_of_analyzed_structures, 'structures in total got analyzed')
        print('Total runtime of StructMAn:', (tend - t0))

    '''
    total_memory_peak = 0
    for p in mem_tracked_processes:
        total_memory_peak += p.memory_info().rss

    total_memory_peak = total_memory_peak/1024./1024./1024.

    print('Accumulated memory peak',total_memory_peak,'Gb')
    '''

    return session
