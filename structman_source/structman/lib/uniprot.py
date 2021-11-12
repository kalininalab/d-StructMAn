import os
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

from Bio import Entrez
import pymysql as MySQLdb

from structman.lib import database, sdsc
from structman.utils import unpack

def is_connected():
    try:
        # connect to the host -- tells us if the host is actually
        # reachable
        socket.create_connection(("1.1.1.1", 53))
        return True
    except OSError:
        pass
    return False


def connection_sleep_cycle(verbosity):
    while not is_connected():
        if verbosity >= 1:
            print('No connection, sleeping a bit and then try again')
        time.sleep(30)


def getUniprotId(query, querytype, verbosity=0):
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
        'from': '%s' % (querytype),
        'to': 'ID',
        'format': 'tab',
        'query': '%s' % (query)
    }
    # print params
    connection_sleep_cycle(verbosity)
    data = urllib.parse.urlencode(params).encode('utf-8')
    request = urllib.request.Request(url, data)
    contact = ""  # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    try:
        response = urllib.request.urlopen(request, timeout=60)
    except:
        return None
    page = response.read(200000).decode('utf-8')
    # print page
    try:
        lines = page.split("\n")
        line = lines[1].split()
        uniprot_id = line[1]
    # If unusable result, try without version number
    except:
        query = query.split(".")
        if len(query) > 1:
            query = query[0]
            return getUniprotId(query, querytype)
        else:
            return None

    return uniprot_id


def tag_update(tag_map, u_ac, new_entry):
    if u_ac not in tag_map:
        tag_map[u_ac] = new_entry
        return tag_map
    else:
        old_entry = tag_map[u_ac]

    for aac in new_entry:
        if aac not in old_entry:
            tag_map[u_ac][aac] = new_entry[aac]
        else:
            tags = set([])
            for tag in old_entry[aac].split(','):
                tags.add(tag)
            for tag in new_entry[aac].split(','):
                tags.add(tag)
            tag_map[u_ac][aac] = ','.join(tags)
    return tag_map

# deactive, needs update
def updateMappingDatabase(u_acs, db, config):
    cursor = db.cursor()
    ac_id_values = []
    ac_ref_values = []
    seq_values = []
    for u_ac in u_acs:
        seq_out = getSequence(u_ac, config, return_id=True)
        if seq_out is None:
            continue

        seq, refseqs, go_terms, pathways, u_id = seq_out
        if u_id is None:  # This can happen for uniprot entries, which got deleted from uniprot
            print("Warning: Uniprot entry:", u_ac, " not found, most probably the entry got deleted from uniprot")
            continue
        ac_id_values.append("('%s','%s','%s','%s')" % (u_ac, u_ac[-2:], u_id, u_id[:2]))
        for refseq in refseqs.split(','):
            ac_ref_values.append("('%s','%s','%s','%s')" % (u_ac, u_ac.split('-')[0][-2:], refseq, refseq[:2]))
        seq_values.append("('%s','%s')" % (u_ac, seq))

    # Don't insert into database if in lite mode
    if config.lite:
        return

    try:

        if len(ac_id_values) > 0:

            sql = "INSERT IGNORE INTO AC_ID(Uniprot_Ac,Ac_Hash,Uniprot_Id,Id_Hash) VALUES %s " % (','.join(ac_id_values))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Error in updateMappingDatabase: %s\n%s" % (sql, f))

        if len(ac_ref_values) > 0:

            sql = "INSERT IGNORE INTO AC_Refseq(Uniprot_Ac,Ac_Hash,Refseq,Refseq_Hash) VALUES %s" % (','.join(ac_ref_values))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Error in updateMappingDatabase: %s\n%s" % (sql, f))

        if len(seq_values) > 0:

            sql = "INSERT IGNORE INTO Sequences(Uniprot_Ac,Sequence) VALUES %s " % (','.join(seq_values))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e, f, g] = sys.exc_info()
                raise NameError("Error in updateMappingDatabase: %s\n%s" % (sql, f))
    except:
        # this happens for Devs without insert rights to the mapping DB, just ignore for the moment
        pass


# called by serializedPipeline
def u_ac_isoform_search(gene_sequence_map, stems, ref_stem_map, config):
    isoform_map = get_all_isoforms(stems, config)
    isoform_specific_id_map = {}
    for nm_ref in gene_sequence_map:
        if nm_ref not in ref_stem_map:
            continue
        seq = gene_sequence_map[nm_ref]
        u_ac_stem = ref_stem_map[nm_ref]
        if u_ac_stem not in isoform_map:
            continue
        for u_ac in isoform_map[u_ac_stem]:
            if seq == isoform_map[u_ac_stem][u_ac][0]:
                isoform_specific_id_map[nm_ref] = u_ac
                break
    return isoform_specific_id_map

# called by serializedPipeline


def integrate_protein(config, proteins, indels, primary_protein_id, input_id, prot_map, u_ac=None, u_id=None, ref=None, pdb_id = None, other_ids={}, is_pdb_input = False):

    if primary_protein_id not in proteins:
        protein = sdsc.Protein(config.errorlog, primary_protein_id=primary_protein_id, u_ac=u_ac, u_id=u_id, ref_id=ref, positions=prot_map[input_id][0], input_id=input_id, pdb_id = pdb_id)
        proteins[primary_protein_id] = protein
    else:
        proteins[primary_protein_id].u_id = u_id
        if is_pdb_input:
            proteins[primary_protein_id].add_residues(prot_map[input_id][0])
        else:
            proteins[primary_protein_id].add_positions(prot_map[input_id][0])
    for other_id_type in other_ids:
        proteins[primary_protein_id].add_other_ids(other_id_type, other_ids[other_id_type])

    if len(prot_map[input_id][1]) > 0:
        indel_insert(config, proteins, indels, prot_map[input_id][1], primary_protein_id)

    for multi_mutation in prot_map[input_id][2]:
        if len(multi_mutation) > 1:
            proteins[primary_protein_id].add_multi_mutation(multi_mutation)
    return


def indel_insert(config, proteins, indel_map, indels, ac):
    if ac not in indel_map:
        indel_map[ac] = []
    for indel in indels:
        indel.wt_prot = ac
        indel_protein_name = indel.create_protein_name(ac)
        if indel_protein_name in proteins:
            broken = False
            while indel_protein_name in proteins:
                indel_protein_name = indel.create_other_protein_name(ac, indel_protein_name)
                if indel_protein_name is None:  # Duplicate entry
                    #config.errorlog.add_error('All indel protein names are reserved, duplicate indels in input?')
                    broken = True
                    break
            if broken:
                continue
        indel_mut_protein = sdsc.Protein(config.errorlog, primary_protein_id=indel_protein_name)
        proteins[indel_protein_name] = indel_mut_protein
        indel.set_proteins(ac, indel_protein_name)
        indel_map[ac].append(indel)
    return

# called by serializedPipeline


def IdMapping(config, ac_map, id_map, np_map, pdb_map, hgnc_map, nm_map):

    proteins = {}
    indel_map = {}

    for ac in ac_map:
        integrate_protein(config, proteins, indel_map, ac, ac, ac_map, u_ac=ac)

    # Step one: map everything to uniprot-ac
    if len(id_map) > 0:
        if config.mapping_db_is_set:

            results = database.select(config, ['Uniprot_Ac', 'Uniprot_Id'], 'UNIPROT', in_rows={'Uniprot_Id':id_map.keys()}, from_mapping_db = True)

            stored_ids = set()
            for row in results:
                u_ac = row[0]
                u_id = row[1]
                stored_ids.add(u_id)
                integrate_protein(config, proteins, indel_map, u_ac, u_id, id_map, u_ac=u_ac, u_id=u_id)

            unstored_ids = []
            for u_id in id_map:
                if u_id not in stored_ids:
                    unstored_ids.append(u_id)
            update_acs = []
            if len(unstored_ids) > 0:  # whenever ids are left in the dict, go to uniprot and download all unmapped entries (this happens for newer uniprot entries, which are not yet in the local mapping database)

                # This part is identical to the part, when no local database is used
                id_ac_map = getUniprotIds(config, unstored_ids, 'ID', target_type="ACC")
                if not id_ac_map is None:
                    for u_id in id_ac_map:
                        u_ac = id_ac_map[u_id]
                        integrate_protein(config, proteins, indel_map, u_ac, u_id, id_map, u_ac=u_ac, u_id=u_id)

                        # This part is different
                        if not sdsc.is_mutant_ac(u_ac):
                            update_acs.append(u_ac)
                    # updateMappingDatabase(update_acs,db,config)

        else:
            id_ac_map = getUniprotIds(config, list(id_map.keys()), 'ID', target_type="ACC")
            for u_id in id_ac_map:
                u_ac = id_ac_map[u_id]
                integrate_protein(config, proteins, indel_map, u_ac, u_id, id_map, u_ac=u_ac, u_id=u_id)

    if len(np_map) > 0:
        if config.mapping_db_is_set:
            results = database.select(config, ['Uniprot_Ac', 'Refseq'], 'UNIPROT', in_rows={'Refseq': np_map.keys()}, from_mapping_db=True)

            ref_u_ac_map = {}  # different u_acs may have the same refseq, try to choose the right one, prefering u_acs containing '-'
            gene_id_snap = set(proteins.keys())  # snapshot of u_acs not originating from refseq-mapping

            stored_refs = set()

            for row in results:
                u_ac = row[0]
                ref = row[1]
                stored_refs.add(ref)
                if ref not in ref_u_ac_map:
                    ref_u_ac_map[ref] = u_ac
                    integrate_protein(config, proteins, indel_map, ref, ref, np_map, u_ac=u_ac)

                elif u_ac in gene_id_snap:
                    if ref_u_ac_map[ref].count('-') == 0 and u_ac.count('-') > 0:
                        ref_u_ac_map[ref] = u_ac
                    integrate_protein(config, proteins, indel_map, ref, ref, np_map, u_ac=u_ac)

                elif ref_u_ac_map[ref].count('-') == 0:
                    if u_ac.count('-') > 0:
                        # if the current u_ac does not contain a '-' and the new found u_ac contains a '-': replace the corresponding ids
                        old_ac = ref_u_ac_map[ref]
                        ref_u_ac_map[ref] = u_ac
                        del proteins[old_ac]
                        integrate_protein(config, proteins, indel_map, ref, ref, np_map, u_ac=u_ac)
            # similar to uniprot-id mapping, we have to go to uniprot to get search for unstored refseq entries and if we find them, we have to update the local mapping database
            unstored_refs = []
            for ref in np_map:
                if ref not in stored_refs:
                    unstored_refs.append(ref)
            update_acs = []
            if len(unstored_refs) > 0:
                np_ac_map = getUniprotIds(config, unstored_refs, 'P_REFSEQ_AC', target_type="ACC")
                for ref in np_ac_map:
                    u_ac = np_ac_map[ref]
                    integrate_protein(config, proteins, indel_map, ref, ref, np_map, u_ac=u_ac)
                    if not sdsc.is_mutant_ac(u_ac):
                        update_acs.append(u_ac)
                # updateMappingDatabase(update_acs,db,config)

        else:
            np_ac_map = getUniprotIds(config, list(np_map.keys()), 'P_REFSEQ_AC', target_type="ACC")
            for ref in np_ac_map:
                u_ac = np_ac_map[ref]
                integrate_protein(config, proteins, indel_map, ref, ref, np_map, u_ac=u_ac)

    if len(nm_map) > 0:
        iso_unspec_nm_keys = [x.split('.')[0] for x in nm_map.keys()]
        nm_ac_map = getUniprotIds(config, iso_unspec_nm_keys, 'REFSEQ_NT_ID', target_type="ACC")
        for ref in nm_map:
            iso_unspec_nm = ref.split('.')[0]
            if iso_unspec_nm in nm_ac_map:
                u_ac = nm_ac_map[iso_unspec_nm]
                integrate_protein(config, proteins, indel_map, ref, ref, nm_map, u_ac=u_ac)

            else:
                integrate_protein(config, proteins, indel_map, ref, ref, nm_map)

    if len(hgnc_map) > 0:  # No support for mapping DB yet
        hgnc_ac_map = getUniprotIds(config, list(hgnc_map.keys()), 'HGNC_ID', target_type="ACC")
        for hgnc in hgnc_ac_map:
            u_ac = hgnc_ac_map[hgnc]
            integrate_protein(config, proteins, indel_map, u_ac, hgnc, hgnc_map, u_ac=u_ac, other_ids={'HGNC_ID': hgnc})

    # Step two: get uniprot-id and refseqs from uniprot-ac

    ac_iso_map = {}
    id_search = []
    for primary_protein_id in proteins:
        if proteins[primary_protein_id].u_id is not None:
            continue
        if primary_protein_id.count('_') > 0:
            continue
        split = primary_protein_id.split('-')
        if len(split) == 2:
            base, iso = split
            if base not in ac_iso_map:
                ac_iso_map[base] = [iso]
            else:
                ac_iso_map[base].append(iso)

            id_search.append(base)
        else:
            id_search.append(primary_protein_id)

    if len(id_search) > 0:
        if config.mapping_db_is_set:
            stored_u_acs = set()

            results = database.select(config, ['Uniprot_Ac', 'Uniprot_Id'], 'UNIPROT', in_rows={'Uniprot_Ac': id_search}, from_mapping_db=True)

            for row in results:
                u_ac = row[0]
                u_id = row[1]

                stored_u_acs.add(u_ac)
                if u_ac in proteins:
                    proteins[u_ac].u_id = u_id

                if u_ac in ac_iso_map:
                    for iso in ac_iso_map[u_ac]:
                        proteins['%s-%s' % (u_ac, iso)].u_id = u_id
                        stored_u_acs.add('%s-%s' % (u_ac, iso))

            unstored_u_acs = []
            for u_ac in id_search:
                if u_ac not in stored_u_acs:
                    if not sdsc.is_mutant_ac(u_ac):
                        unstored_u_acs.append(u_ac)

            if len(unstored_u_acs) > 0:
                # updateMappingDatabase(unstored_u_acs, db, config)

                results = database.select(config, ['Uniprot_Ac', 'Uniprot_Id'], 'UNIPROT', in_rows={'Uniprot_Ac': unstored_u_acs}, from_mapping_db=True)

                for row in results:
                    u_ac = row[0]
                    u_id = row[1]

                    if u_ac in proteins:
                        proteins[u_ac].u_id = u_id
                    if u_ac in ac_iso_map:
                        for iso in ac_iso_map[u_ac]:
                            proteins['%s-%s' % (u_ac, iso)].u_id = u_id

        else:
            ac_id_map = getUniprotIds(config, id_search, 'ACC', target_type="ID")
            if ac_id_map is not None:
                for u_ac in ac_id_map:
                    u_id = ac_id_map[u_ac]
                    if u_ac in proteins:
                        proteins[u_ac].u_id = u_id
                    if u_ac in ac_iso_map:
                        for iso in ac_iso_map[u_ac]:
                            proteins['%s-%s' % (u_ac, iso)].u_id = u_id

    if config.mapping_db_is_set:

        results = database.select(config, ['Uniprot_Ac', 'Refseq'], 'UNIPROT', in_rows={'Uniprot_Ac': proteins.keys()}, from_mapping_db=True)

        for row in results:
            u_ac = row[0]
            ref = row[1]
            proteins[u_ac].ref_id = ref
    else:
        ac_np_map = getUniprotIds(config, list(proteins.keys()), 'ACC', target_type="P_REFSEQ_AC")
        for u_ac in ac_np_map:
            ref = ac_np_map[u_ac]
            proteins[u_ac].ref_id = ref

    for pdb_tuple in pdb_map:
        integrate_protein(config, proteins, indel_map, pdb_tuple, pdb_tuple, pdb_map, is_pdb_input = True, pdb_id = pdb_tuple)

    return proteins, indel_map


def getUniprotIds(config, query_ids, querytype, target_type = "ID", timeout = 60):
    # print query_ids
    if len(query_ids) == 0:
        return {}

    if config.verbosity >= 4:
        print('Call of getUniprotIds:', query_ids, querytype, target_type)

    query = ' '.join(query_ids)
    url = 'https://www.uniprot.org/uploadlists/'
    connection_sleep_cycle(config.verbosity)
    params = {
        'from': '%s' % (querytype),
        'to': '%s' % (target_type),
        'format': 'tab',
        'query': '%s' % (query)
    }
    # print params
    data = urllib.parse.urlencode(params).encode('utf-8')
    request = urllib.request.Request(url, data)
    contact = config.user_mail  # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    try:
        response = urllib.request.urlopen(request, timeout = timeout)
    except:
        e, f, g = sys.exc_info()
        config.errorlog.add_warning("Uniprot did not answer: %s\n%s" % (str(e), str(f)))
        return {}

    page = response.read(2000000).decode('utf-8')

    if config.verbosity >= 4:
        print('Uniprot answer:\n',page)

    uniprot_ids = {}
    try:
        error = True
        lines = page.split("\n")
        for line in lines[1:]:
            words = line.split()
            if len(words) > 1:
                error = False
                quer = words[0]
                target = words[1]
                if quer not in query_ids:
                    error = True
                    break
                if quer == target:
                    target = words[3].split(';')[0]
                uniprot_ids[quer] = target
        if error:

            if len(query_ids) == 1:
                quer = query_ids.pop()
                if querytype == 'ACC' and quer.count('-'):
                    int_quer = set([quer.split('-')[0]])
                    intermediate_map = getUniprotIds(config, int_quer, querytype, target_type=target_type)
                    return {quer: intermediate_map[quer.split('-')[0]]}
                else:
                    return {quer: None}
            # try a divide and conqer solution
            set_A = set()
            set_B = set()
            while len(query_ids) > 0:
                set_A.add(query_ids.pop())
                if len(query_ids) > 0:
                    set_B.add(query_ids.pop())
            map_A = getUniprotIds(config, set_A, querytype, target_type=target_type)
            map_B = getUniprotIds(config, set_B, querytype, target_type=target_type)
            return map_A.update(map_B)

    # If unusable result, try without version number
    except:
        # This seems not up-to-date
        '''
        if query.count(".") > 0:
            ids = query.split("\t")
            nids = []
            for i in ids:
                ni = i.split(".")[0]
                nids.append(ni)
            nquery = "\t".join(nids)
            if nquery != query:
                return getUniprotId(nquery,querytype,target_type=target_type)
            else:
                return {}
        '''
        config.errorlog.add_warning(f'Parse error of uniprot ID-mapping webservice\n{page}')
        return {}
    return uniprot_ids

# First version of refseq sequence retrieval. Uses biopython.
# called by serializePipeline


def get_refseq_sequences(refseqs, config, seq_type='nucleotide'):
    Entrez.email = config.user_mail

    ret_type = 'fasta_cds_aa'
    if seq_type == 'protein':
        ret_type = 'fasta'

    net_handle = Entrez.efetch(
        db=seq_type, id=refseqs, rettype=ret_type, retmode="text"
    )
    page = net_handle.read()
    net_handle.close()
    right_split = '_prot'
    left_split = '|'
    if seq_type == 'protein':
        right_split = ' '
        left_split = None
    seq_map = sdsc.parseFasta(page=page, left_split=left_split, right_split=right_split)

    return seq_map


def translateGNSMap(gene_nuc_sequence_map):
    gene_sequence_map = {}
    for ref in gene_nuc_sequence_map:
        nuc_seq = gene_nuc_sequence_map[ref]
        aa_seq = sdsc.translate(nuc_seq)
        gene_sequence_map[ref] = aa_seq
    return gene_sequence_map


def get_all_isoforms(u_ac_stems, config):

    could_find_more = True
    isoform_counter = 1
    isoform_counter_steps = 20
    isoform_map = {}
    while could_find_more:
        u_acs = []
        last_isoforms = {}
        for u_ac_stem in u_ac_stems:
            isoform_map[u_ac_stem] = {}
            stem_isoform_counter = isoform_counter
            for i in range(isoform_counter_steps):
                if stem_isoform_counter == 1:
                    u_acs.append(u_ac_stem)
                else:
                    isoform = '%s-%s' % (u_ac_stem, str(stem_isoform_counter))
                    u_acs.append(isoform)
                    if stem_isoform_counter == ((isoform_counter + isoform_counter_steps) - 1):
                        last_isoforms[u_ac_stem] = isoform
                stem_isoform_counter += 1
        isoform_counter += isoform_counter_steps

        gene_sequence_map = getSequencesPlain(u_acs, config, save_errors=False, skip_missing_routine=True)
        for u_ac in gene_sequence_map:
            u_ac_stem = u_ac.split('-')[0]
            isoform_map[u_ac_stem][u_ac] = gene_sequence_map[u_ac]
        could_find_more = False
        for u_ac_stem in last_isoforms:
            last_isoform = last_isoforms[u_ac_stem]
            if last_isoform in isoform_map[u_ac_stem]:
                could_find_more = True
            else:
                u_ac_stems.remove(u_ac_stem)

    return isoform_map

# called by serializedPipeline


def getSequencesPlain(u_acs, config, max_seq_len=None, filtering_db=None, save_errors=True, skip_missing_routine=False):
    gene_sequence_map = {}
    filtered_set = set()

    missing_set = set()

    if filtering_db is not None:
        filtering_db_path, dbs = filtering_db

    if config.mapping_db_is_set:
        if len(u_acs) == 0:
            return {}
        if config.verbosity >= 2:
            t0 = time.time()

        results = database.select(config, ['Uniprot_Ac', 'Sequence'], 'UNIPROT', in_rows={'Uniprot_Ac': u_acs}, from_mapping_db = True)

        if config.verbosity >= 2:
            t1 = time.time()
            print("getSequencesPlain Part 1: ", str(t1 - t0))

        n = 0
        for row in results:
            u_ac = row[0]
            if u_ac not in u_acs:
                continue
            if row[1] is None:
                config.errorlog.add_warning('Sequence field is None in Mapping DB, for: %s' % u_ac) 
                #gene_sequence_map[u_ac] = None, None, None
                continue
            try:
                seq = unpack(row[1])
            except:
                config.errorlog.add_warning('Sequence field is defect in Mapping DB, for: %s' % u_ac)
                #gene_sequence_map[u_ac] = None, None, None
                continue

            if max_seq_len is not None:
                if len(seq) > max_seq_len:
                    filtered_set.add(u_ac)
                    n += 1
                    continue
            
            disorder_scores = None
            disorder_regions_datastruct = None

            gene_sequence_map[u_ac] = seq, disorder_scores, disorder_regions_datastruct

        if n > 0 and config.verbosity >= 2:
            print('Filtered ', n, ' Sequences due to max length: ', max_seq_len)

        if config.verbosity >= 2:
            t2 = time.time()
            print("getSequencesPlain Part 2: ", str(t2 - t1))

    else:
        missing_set = set(u_acs)

    if config.verbosity >= 2:
        t2 = time.time()

    in_db = set()
    for u_ac in u_acs:
        if filtering_db is not None:
            inside_all = True
            for db_name in dbs:
                folder_key = u_ac.split('-')[0][-2:]
                filename = '%s/%s/%s_%s_gpw.fasta.gz' % (filtering_db_path, folder_key, u_ac, db_name)
                if not os.path.isfile(filename):
                    inside_all = False
            if inside_all:
                in_db.add(u_ac)
                continue

        if u_ac in filtered_set:
            continue

        if u_ac not in gene_sequence_map:
            missing_set.add(u_ac)
            # print u_ac

    if config.verbosity >= 2:
        t3 = time.time()
        print("getSequencesPlain Part 3: ", str(t3 - t2))

    if not skip_missing_routine:

        if config.verbosity >= 2:
            print('Size of missing set: ', len(missing_set))
        if config.verbosity >= 3:
            print(missing_set)

        for u_ac in missing_set:
            seq_out = getSequence(u_ac, config)
            if seq_out is None and save_errors:
                config.errorlog.add_warning('getSequence output is None for %s' % u_ac)
                gene_sequence_map[u_ac] = 0, None, None
                continue
            elif seq_out is None:
                continue
            seq, refseqs, go_terms, pathways = seq_out
            gene_sequence_map[u_ac] = seq, None, None

        if config.verbosity >= 2:
            t4 = time.time()
            print("getSequencesPlain Part 4: ", str(t4 - t3))

    if filtering_db is not None:
        return gene_sequence_map, in_db
    return gene_sequence_map


# called by serializePipeline
def getSequences(proteins, config):

    t0 = time.time()
    u_acs = set()
    iso_map = {}
    protein_map = proteins.get_protein_map()
    for prot_id in protein_map:
        if protein_map[prot_id].sequence is None:
            u_acs.add(prot_id)

        if prot_id.count('-') == 1:
            [base, iso] = prot_id.split('-')
        else:
            base = prot_id
            iso = 'c'
        if base not in iso_map:
            iso_map[base] = [iso]
        else:
            iso_map[base].append(iso)

    gene_sequence_map = getSequencesPlain(u_acs, config)

    for u_ac in u_acs:
        protein_map[u_ac].sequence = gene_sequence_map[u_ac][0]
        proteins.set_disorder_scores(u_ac, gene_sequence_map[u_ac][1])
        proteins.set_disorder_regions(u_ac, gene_sequence_map[u_ac][2])

    info_map_path = '%s/human_info_map.tab' % config.human_id_mapping_path

    f = open(info_map_path, 'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        words = line[:-1].split('\t')
        u_ac = words[0]
        if u_ac in iso_map:
            go_terms = {}
            pathways = {}
            go_str = words[1]
            if go_str != '':
                go_doubles = go_str[:-1].split(';')
                for go_double in go_doubles:
                    go_id = (':'.join(go_double.split(':')[:2])).strip()
                    go_term = ':'.join(go_double.split(':')[2:])
                    go_terms[go_id] = go_term
            path_str = words[2]
            if path_str != '':
                path_doubles = path_str[:-1].split(';')
                for path_double in path_doubles:
                    [reac_id, pathway] = path_double.split(':', 1)
                    pathways[reac_id] = pathway
            for iso in iso_map[u_ac]:
                if iso == 'c':
                    iso_u_ac = u_ac
                else:
                    iso_u_ac = '%s-%s' % (u_ac, iso)
                protein_map[iso_u_ac].go_terms = go_terms
                protein_map[iso_u_ac].pathways = pathways

def get_last_version(u_ac):
    url = 'https://www.uniprot.org/uniprot/%s?version=*' % u_ac
    #connection_sleep_cycle(config.verbosity)
    try:
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request)
        page = response.read(9000000).decode('utf-8')
    except:
        e, f, g = sys.exc_info()
        config.errorlog.add_error('Error trying to reach: %s\n%s\n%s' % (url, str(e), str(f)))
        return None
    print(page)

def get_obsolete_sequence(u_ac, config, return_id=False, tries = 0):
    uniparc_id_map = getUniprotIds(config, [u_ac], 'ACC', target_type="UPARC", timeout = 60 * (tries + 1))
    if not u_ac in uniparc_id_map:
        if tries < 3:
             return get_obsolete_sequence(u_ac, config, return_id=return_id, tries = (tries + 1))
        config.errorlog.add_warning(f'Couldnt find uniparc id for: {u_ac}, {uniparc_id_map}')
        return None
    uniparc_id = uniparc_id_map[u_ac]
    if config.verbosity >= 3:
        print('Uniparc ID of potential obsolete uniprot entry (', u_ac, ') is:', uniparc_id)
    return getSequence(uniparc_id, config, return_id=return_id, obsolete_try = True)

def getSequence(uniprot_ac, config, tries=0, return_id=False, obsolete_try = False):
    if config.verbosity >= 3:
        print('uniprot.getSequence for ', uniprot_ac)

    if sdsc.is_mutant_ac(uniprot_ac):
        config.errorlog.add_error('Cannot call getSequence with a mutant protein: %s' % uniprot_ac)
        return None

    # new part just for the sequence
    if len(uniprot_ac) < 2:
        return None

    if not uniprot_ac[0:3] == 'UPI':
        url = 'https://www.uniprot.org/uniprot/%s.fasta' % uniprot_ac
    else:
        url = 'https://www.uniprot.org/uniparc/%s.fasta' % uniprot_ac
    connection_sleep_cycle(config.verbosity)
    try:
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request, timeout=(tries + 1) * 10)
        page = response.read(9000000).decode('utf-8')
    except:
        if tries < 3:
            return getSequence(uniprot_ac, config, tries=tries + 1, return_id=return_id)
        else:
            return get_obsolete_sequence(uniprot_ac, config, return_id=return_id)

    if config.verbosity >= 3 and obsolete_try:
        print('Obsolete entry try:', uniprot_ac, url)
        if config.verbosity >= 4:
            print('Online return:\n', page)

    lines = page.split("\n")

    wildtype_sequences = []

    for line in lines:
        if line == '':
            continue
        if line[0] == '>':
            continue
        wildtype_sequences.append(line)

    wildtype_sequence = ("".join(wildtype_sequences)).replace(" ", "").replace("\n", "")

    if wildtype_sequence == '':
        if obsolete_try: #Unless it could come to an endless loop
            return None
        if config.verbosity >= 3:
            print('Try to search for obsolete uniprot entry:', uniprot_ac)
        return get_obsolete_sequence(uniprot_ac, config, return_id=return_id)

    if uniprot_ac[0:3] == 'UPI':
        if return_id:
            return (wildtype_sequence, {}, {}, {}, None)
        return (wildtype_sequence, {}, {}, {})

    # old part, now just for refseqs,go and reactome

    url = 'https://www.uniprot.org/uniprot/%s.txt' % uniprot_ac
    try:
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request, timeout=(tries + 1) * 10)
        page = response.read(9000000).decode('utf-8')
    except:
        if tries < 3:
            return getSequence(uniprot_ac, config, tries=tries + 1, return_id=return_id)

        else:
            # print uniprot_ac
            return None

    lines = page.split("\n")
    # print(lines)
    refseqs = {}
    go_terms = {}
    pathways = {}
    u_id = None
    for line in lines:
        if line == '':
            continue
        words = line.split()
        if len(words) == 0:
            print(uniprot_ac)
            return None
        if words[0] == 'ID':
            u_id = words[1]
        if words[0] == "DR":
            if len(words) > 1:
                if words[1] == "GO;":
                    words_ = line.split(";")
                    go_id = words_[1].strip(";").strip()
                    go_name = words_[2].strip(";")
                    go_terms[go_id] = go_name
                if words[1] == "Reactome;":
                    split = line.split(";")
                    reac_id = split[1].replace(" ", "")
                    pathway = split[2][1:-1]
                    pathways[reac_id] = pathway
                # TODO filter out isoforms
                if words[1] == "RefSeq;":
                    if words[-1].count('[') > 0:
                        u_ac_iso = words[-1][1:-1]
                        refs = words[2:-1]
                    else:
                        u_ac_iso = uniprot_ac.split('-')[0]
                        refs = words[2:]
                    refs = [x[:-1] for x in refs]
                    if u_ac_iso not in refseqs:
                        refseqs[u_ac_iso] = []
                    refseqs[u_ac_iso] += refs
    if uniprot_ac in refseqs:
        refseqs = ",".join(refseqs[uniprot_ac])
    else:
        refseqs = ''
    if return_id:
        return (wildtype_sequence, refseqs, go_terms, pathways, u_id)
    return (wildtype_sequence, refseqs, go_terms, pathways)
