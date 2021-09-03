import time

from Bio import pairwise2

from structman.lib import sdsc


def toOne(res_name, only_natural=False):
    if only_natural:
        if res_name in sdsc.AA_TO_ONE:
            return sdsc.AA_TO_ONE[res_name]
        # if res_name in sdsc.NUCLEOTIDE_TO_ONE:
        #    return sdsc.NUCLEOTIDE_TO_ONE[res_name]
        return '.'
    if res_name in sdsc.THREE_TO_ONE:
        return sdsc.THREE_TO_ONE[res_name][0]
    elif res_name in sdsc.NUCLEOTIDE_TO_ONE:
        return sdsc.NUCLEOTIDE_TO_ONE[res_name]
    else:
        return 'X'


def consecutive(last_residue, res_nr):
    if last_residue is None:
        return True
    try:
        a = int(last_residue)
        b = int(res_nr)
    except:  # insertion codes
        return True
    return (b - a) == 1

# called by serializedPipeline, sdsc


def createTemplateFasta(template_page, template_name, chain, config, onlySeqResMap=False, seqAndMap=False, for_modeller=False, could_be_empty=False, rare_residues=None):
    lines = template_page.split('\n')
    if rare_residues is None:
        rare_residues = set()
    seq = ""
    seq_res_map = []
    used_res = set()
    last_residue = None
    first_residue = None
    res_atom_count = 0
    for line in lines:
        record_name = line[0:6].replace(" ", "")
        atom_nr = line[6:11].replace(" ", "")
        atom_name = line[12:16].replace(" ", "")
        # ignore water
        if len(atom_name) > 0:
            if atom_name[0] == 'H':
                continue
        res_name = line[17:20].replace(" ", "")
        if len(line) > 21:
            chain_id = line[21]
            if chain_id == ' ':
                chain_id = chain
            res_nr = line[22:27].replace(" ", "")  # this includes the insertion code

            if record_name == 'SEQRES':
                for tlc in line[19:].split():
                    tlc = tlc.strip()
                    if len(tlc) != 3:
                        continue
                    if tlc not in sdsc.THREE_TO_ONE:
                        rare_residues.add(tlc)

            if chain != chain_id:
                continue

            if record_name == "ATOM" or record_name == 'HETATM':
                if record_name == 'HETATM':
                    last_residue = res_nr
                    if (res_name not in sdsc.THREE_TO_ONE) and (res_name not in rare_residues):
                        continue
                if for_modeller:
                    if res_nr not in used_res:
                        if len(seq) > 0:
                            last_res = seq[-1]
                            if last_res != '-' and last_res != '.':
                                if sdsc.CORRECT_COUNT[last_res] > res_atom_count:
                                    #seq = '%s-' % seq[:-1]
                                    seq += '-'
                                elif sdsc.CORRECT_COUNT[last_res] < res_atom_count:
                                    seq = '%s.' % seq[:-1]
                            res_atom_count = 0
                        if record_name == "ATOM":
                            aa = toOne(res_name, only_natural=True)
                        else:
                            aa = '.'

                        res_atom_count += 1
                        # if consecutive(last_residue,res_nr):
                        seq = seq + aa
                        # else:
                        #    seq += '-%s' % aa
                        seq_res_map.append(res_nr)
                        used_res.add(res_nr)
                        last_residue = res_nr
                        if first_residue is None:
                            first_residue = res_nr
                    elif atom_name != 'OXT':
                        res_atom_count += 1
                else:
                    if res_nr not in used_res:
                        aa = toOne(res_name)
                        if aa not in sdsc.ONE_TO_THREE:
                            aa = 'X'
                        seq = seq + aa
                        seq_res_map.append(res_nr)
                        used_res.add(res_nr)
                        last_residue = res_nr
                        if first_residue is None:
                            first_residue = res_nr

    if for_modeller:
        if len(seq) > 0:
            last_res = seq[-1]
            if last_res != '-' and last_res != '.':
                if sdsc.CORRECT_COUNT[last_res] > res_atom_count:
                    #seq = '%s-' % seq[:-1]
                    seq += '-'
                elif sdsc.CORRECT_COUNT[last_res] < res_atom_count:
                    seq = '%s.' % seq[:-1]

    if seq_res_map == [] and not could_be_empty:
        config.errorlog.add_warning('Warning: seq_res_map empty: %s:%s' % (template_name, chain))

    if onlySeqResMap:
        return seq_res_map, last_residue, first_residue
    if seqAndMap:
        return seq_res_map, seq, last_residue, first_residue

    page = ">" + template_name + "\n" + seq

    template_fasta = "tmp.%s.fasta" % template_name
    f = open(template_fasta, "wb")
    f.write(page)
    f.close()

    return (template_fasta, seq_res_map)


# gets also called by serializedPipeline
def getSubPos(config, u_ac, target_aligned_sequence, template_aligned_sequence, aaclist, seq_res_map, ignore_gaps=False, lock=None):
    target_aligned_sequence = target_aligned_sequence.replace("\n", "")
    template_aligned_sequence = template_aligned_sequence.replace("\n", "")

    sub_infos = {}

    errors = []
    align_map = {}
    n = 0
    tar_n = 0
    tem_n = 0

    for char in target_aligned_sequence:
        tem_char = template_aligned_sequence[n]
        n += 1
        if tem_char != '-':
            tem_n += 1
        if char != '-':
            tar_n += 1
            if tem_char == '-':
                align_map[tar_n] = ((None, '-', None), char)
            else:
                if tem_n - 1 >= len(seq_res_map):
                    return ("Seq_res_map too short: %s, %s, %s" % (u_ac, tem_n, len(seq_res_map)))
                align_map[tar_n] = ((seq_res_map[tem_n - 1], tem_char, tem_n), char)

    for aac_base in aaclist:
        target_pos = int(aac_base[1:])
        target_aa = aac_base[0]
        if target_pos not in align_map:
            error = 'Mutation not inside target sequence: %s %s' % (u_ac, aac_base)
            config.errorlog.add_warning(error, lock=lock)
            continue
        sub_info = align_map[target_pos][0]

        if align_map[target_pos][1] != target_aa:

            if target_aa != '?':
                if not (target_aa == 'U' and align_map[target_pos][1] == 'C'):
                    if not (target_aa == 'O' and align_map[target_pos][1] == 'K'):
                        if not (target_aa == 'J' and align_map[target_pos][1] == 'I'):
                            config.errorlog.add_warning('Amino acid of Mutation does not match with amino acid in the query protein: %s; %s%s (found),%s (given)' % (u_ac, align_map[target_pos][1], str(target_pos), aac_base), lock=lock)

        if ignore_gaps and sub_info[0] is None:
            continue
        sub_infos[target_pos] = sub_info

    return (sub_infos, aaclist)


def getCovSI(full_length, target_seq, template_seq):
    target_length = len(target_seq.replace("-", ""))
    template_length = float((target_length - template_seq.count("-")))
    if template_length == 0.:
        return None, None
    aln_length = template_length / float(full_length)
    i = 0
    identical = 0
    for res_a in target_seq:
        if i == len(template_seq):
            print(target_seq + '\n')
            print(template_seq)
        if template_seq[i] == res_a:
            identical += 1
        i += 1
    seq_id = float(identical) / template_length

    return (aln_length, seq_id)


def BPalign(config, u_ac, target_seq, template_seq, aaclist, seq_res_map, ignore_gaps=False, lock=None):
    target_seq = target_seq.replace('U', 'C').replace('O', 'K').replace('J', 'I')
    if config.verbosity >= 4:
        print('Aligning:')
        print(target_seq)
        print(template_seq)
    align_out = pairwise2.align.globalds(target_seq, template_seq, sdsc.BLOSUM62, -10.0, -0.5, one_alignment_only=True, penalize_end_gaps=False)
    if len(align_out) > 0:
        (target_aligned_sequence, template_aligned_sequence, a, b, c) = align_out[0]
    else:
        if config.verbosity >= 4:
            print('Alignment failed')
        return 'Alignment produced no Output, %s:\n%s\n%s\n' % (u_ac, target_seq, template_seq)

    sub_infos, aaclist = getSubPos(config, u_ac, target_aligned_sequence, template_aligned_sequence, aaclist, seq_res_map, ignore_gaps=ignore_gaps, lock=lock)

    return (target_aligned_sequence, template_aligned_sequence, sub_infos, aaclist)


# truncates aseq by looking at terminal gaps of bseq
# called by modelling
def truncateSequences(aseq, bseq):
    aseq = aseq.replace("\n", "")
    bseq = bseq.replace("\n", "")
    full = len(bseq)
    bseq = bseq.lstrip("-")
    left = full - len(bseq)
    bseq = bseq.rstrip("-")
    right = full - left - len(bseq)
    aseq = aseq[left:(len(aseq) - right)]
    bseq = bseq
    return (aseq, bseq)


def createAlignmentPir(target_name, target_aligned_sequence, template_name, template_aligned_sequence, chain, startres, endres, outfile=''):
    # print outfile
    page = ">P1;%s\nsequence:%s:::::::0.00: 0.00\n%s*\n>P1;%s\nstructureX:%s:%s:%s:%s:%s:::-1.00:-1.00\n%s*" % (target_name, target_name, target_aligned_sequence, template_name, template_name, startres, chain, endres, chain, template_aligned_sequence)
    # print(page)
    #f = open(outfile, "wb")
    # f.write(page)
    # f.close()
    # return(outfile)
    return page


# called by serializedPipeline
def alignBioPython(config, target_name, wildtype_sequence, template_name, template_page, chain, aaclist, ignore_gaps=False, lock=None, rare_residues=None):
    # preparing the alignment of the target and the template, by:
    t0 = time.time()
    (seq_res_map, template_seq, last_residue, first_residue) = createTemplateFasta(template_page, template_name, chain, config, seqAndMap=True, rare_residues=rare_residues)
    if len(seq_res_map) == 0:
        return 'Unable to create template fasta %s %s %s' % (target_name, template_name, chain)
    startres = seq_res_map[0]
    endres = seq_res_map[-1]

    t1 = time.time()
    align_out = BPalign(config, target_name, wildtype_sequence, template_seq, aaclist, seq_res_map, ignore_gaps=ignore_gaps, lock=lock)

    if isinstance(align_out, str):
        return align_out

    (target_aligned_sequence, template_aligned_sequence, sub_infos, aaclist) = align_out
    t2 = time.time()
    # write the alignment into the pir format, which can be used by the modeller
    (truncated_target_sequence, truncated_template_sequence) = truncateSequences(target_aligned_sequence, template_aligned_sequence)
    t3 = time.time()
    (aln_length, seq_id) = getCovSI(len(target_aligned_sequence), truncated_target_sequence, truncated_template_sequence)
    t4 = time.time()
    alignment_pir = createAlignmentPir(target_name, target_aligned_sequence, template_name, template_aligned_sequence, chain, startres, endres)
    t5 = time.time()

    times = [t1 - t0, t2 - t1, t3 - t2, t4 - t3, t5 - t4]
    return (aln_length, seq_id, sub_infos, alignment_pir, times, aaclist, last_residue, first_residue)
