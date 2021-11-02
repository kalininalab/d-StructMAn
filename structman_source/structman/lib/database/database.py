import cProfile
import multiprocessing
import random
import re
import sys
import time
import traceback

import ray

from structman.utils import ray_init, calc_checksum, pack, unpack
from structman.lib import pdbParser, rin, sdsc, serializedPipeline, templateFiltering


def binningSelect(keys, rows, table, config, binning_function='median_focus', density=0.5):
    # Use this for huge selects, the first entry of rows has to be the key by which the entries are selected
    key_name = rows[0]
    t0 = time.time()
    if binning_function == 'split_fusion':
        singletons, bins = split_fusion_binning(keys, density=density, fusion=True)
    elif binning_function == 'split':
        singletons, bins = split_fusion_binning(keys, density=density, fusion=False)
    elif binning_function == 'median_focus':
        singletons, bins = median_focus_binning(keys)
    else:
        config.errorlog.add_error('Unknown binning function:', binning_function)
        return []

    if config.verbosity >= 6:
        print('\nbinningSelect keys:\n', keys, 'and binning results:\n', singletons, '\n', bins)

    t1 = time.time()
    if config.verbosity >= 3:
        print('Time for binning in binningSelect:', t1 - t0)

    if len(singletons) > 0:
        if len(singletons) == 1:
            equals_rows = {key_name: singletons[0]}
            total_results = list(select(config, rows, table, equals_rows=equals_rows))
        else:
            in_rows = {key_name: singletons}
            total_results = list(select(config, rows, table, in_rows=in_rows))
    else:
        total_results = []

    t2 = time.time()
    if config.verbosity >= 3:
        print('Time for singleton select in binningSelect:', t2 - t1, 'Amount of singletons:', len(singletons))

    t3 = 0.
    t4 = 0.
    t5 = 0.
    for ids, min_id, max_id in bins:
        t3 += time.time()
        between_rows = {key_name: (min_id, max_id)}

        results = select(config, rows, table, between_rows=between_rows)

        t4 += time.time()
        for row in results:
            if not row[0] in ids:
                continue
            total_results.append(row)
        t5 += time.time()

    if config.verbosity >= 3:
        print('Time for between select in binningSelect:', t4 - t3, 'Amount of bins:', len(bins))
        print('Time for id check in binningSelect:', t5 - t4)

    return total_results


def select(config, rows, table, between_rows={}, in_rows={}, equals_rows={}, null_columns=set(), n_trials=3, from_mapping_db=False):
    if len(rows) == 0:
        row_str = '*'
    else:
        row_str = ','.join(rows)

    params = []

    if len(between_rows) == 0 and len(in_rows) == 0 and len(equals_rows) == 0:
        where_str = ''
    else:
        where_parts = []

        for equals_row in equals_rows:
            params.append(equals_rows[equals_row])
            where_parts.append(equals_row + ' = %s')

        for null_column in null_columns:
            where_parts.append(null_column + ' IS NULL')

        for in_row in in_rows:
            for param in in_rows[in_row]:
                params.append(param)
            where_parts.append(in_row + ' IN (%s)' % ','.join(['%s'] * len(in_rows[in_row])))  # There have to be as many %s placeholders in the statement as there are parameters for the IN clasue

        for bet_row in between_rows:
            (low, high) = between_rows[bet_row]
            params.append(low)
            params.append(high)
            where_parts.append(bet_row + ' BETWEEN %s AND %s')

        where_str = ' WHERE %s' % ' AND '.join(where_parts)

        if len(params) == 0:
            return []

    statement = 'SELECT %s FROM %s%s' % (row_str, table, where_str)

    n = 0
    while n < n_trials:  # Repeat the querry if fails for n_trials times
        db, cursor = config.getDB(mapping_db=from_mapping_db)
        try:
            cursor.execute(statement, params)
            results = cursor.fetchall()
            db.commit()
            break
        except:
            if n == 0:
                [e, f, g] = sys.exc_info()
                g = traceback.format_exc()
            n += 1
        db.close()
    if n == n_trials:
        raise NameError('Invalid Select: %s\nParam size:%s\n%s\n%s, %s\n%s\n%s\n%s' % (statement[:500], str(len(params)), str(params[:50]), from_mapping_db, config.mapping_db, e, str(f), g))
    return results


def size_estimation(config, values):
    # Select 100 random rows and calculate their str size
    n_of_r_values = 100
    if len(values) <= n_of_r_values:
        #size_of_values = sys.getsizeof(str(values))
        size_of_values = sdsc.total_size(values)
    else:
        #size_of_values = (sys.getsizeof(str(random.sample(values,n_of_r_values)))*len(values))/n_of_r_values
        size_of_values = (sdsc.total_size(random.sample(values, n_of_r_values)) * len(values)) / n_of_r_values

    # Add X% too the size estimation, just to be sure
    size_of_values = size_of_values * 2.5
    return size_of_values


def insert(table, columns, values, config, n_trials=3, mapping_db = False):
    params = []

    columns_str = ','.join(columns)

    parts = []
    value_strs = []

    if len(values) == 0:
        if config.verbosity >= 3:
            print('Empty insert try', table)
        return

    max_package_size = config.max_package_size
    size_of_values = size_estimation(config, values)

    number_of_packages = (size_of_values // max_package_size)
    if not size_of_values % max_package_size == 0:
        number_of_packages += 1

    package_length = (len(values) // number_of_packages)
    if not len(values) % number_of_packages == 0:
        package_length += 1

    if config.verbosity >= 3:
        print('Insert to', table, 'with total estimated size', size_of_values / 1024. / 1024., 'Mb,since max size is',
              max_package_size / 1024. / 1024., 'Mb this makes', number_of_packages,
              'packages in total (size per package:', size_of_values / 1024. / 1024. / number_of_packages, ')')

    n = 0
    for value_list in values:
        for value in value_list:
            params.append(value)
        value_str_part = '(%s)' % ','.join(['%s'] * len(value_list))
        value_strs.append(value_str_part)
        n += 1
        if n == package_length:
            n = 0
            value_str = ','.join(value_strs)
            parts.append((value_str, params))
            value_strs = []
            params = []
    if value_strs != []:
        value_str = ','.join(value_strs)
        parts.append((value_str, params))

    for value_str, params in parts:
        statement = 'INSERT IGNORE INTO %s (%s) VALUES %s' % (table, columns_str, value_str)

        if config.verbosity >= 2:
            print('Insert with', len(params), 'parameters')

        if config.verbosity >= 3:
            print('Better size estimation of the package:', sys.getsizeof(str(params)) + sys.getsizeof(statement))

        n = 0
        while n < n_trials:  # Repeat the querry if fails for n_trials times
            db, cursor = config.getDB(mapping_db = mapping_db)
            try:
                cursor.execute(statement, params)
                db.commit()
                db.close()
                break
            except:
                db.close()
                n += 1
                if n == 1:
                    [e, f, g] = sys.exc_info()
                    g = traceback.format_exc()
        if n == n_trials:
            raise NameError('Invalid Insert: %s\nParam size:%s\n%s\n%s\n%s\n%s' % (statement[:500], str(len(params)), str(params[:500]), e, str(f), g))


def update(config, table, columns, values, mapping_db = False):
    # Updating with an insert statement (see: https://stackoverflow.com/questions/20255138/sql-update-multiple-records-in-one-query)
    params = []

    column_str = ','.join(columns)

    update_str_parts = []
    for column in columns:
        update_str_parts.append('%s=VALUES(%s)' % (column, column))

    update_str = ','.join(update_str_parts)

    parts = []
    n = 0
    value_strs = []

    if len(values) == 0:
        return

    max_package_size = config.max_package_size
    size_of_values = size_estimation(config, values)

    number_of_packages = (size_of_values // max_package_size)
    if not size_of_values % max_package_size == 0:
        number_of_packages += 1

    package_length = (len(values) // number_of_packages)
    if not len(values) % number_of_packages == 0:
        package_length += 1

    n = 0
    for value_list in values:
        for value in value_list:
            params.append(value)
        value_str_part = '(%s)' % ','.join(['%s'] * len(value_list))
        value_strs.append(value_str_part)
        n += 1
        if n == package_length:
            n = 0
            value_str = ','.join(value_strs)
            parts.append((value_str, params))
            value_strs = []
            params = []
    if value_strs != []:
        value_str = ','.join(value_strs)
        parts.append((value_str, params))

    for value_str, params in parts:
        if len(params) == 0:
            continue
        statement = 'INSERT IGNORE INTO %s (%s) VALUES %s ON DUPLICATE KEY UPDATE %s' % (table, column_str, value_str, update_str)
        db, cursor = config.getDB(mapping_db = mapping_db)
        try:
            cursor.execute(statement, params)
            db.commit()
        except:
            [e, f, g] = sys.exc_info()
            g = traceback.format_exc()
            raise NameError('Invalid Update: %s\nParam size:%s\n%s\n%s\n%s\n%s' % (statement[:500], str(len(params)), str(params[:500]), e, f, g))
        db.close()


# called by output.classification
def getProteinDict(prot_id_list, session_id, config, includeSequence=False):
    if len(prot_id_list) == 0:
        return {}

    table = 'RS_Protein_Session'
    columns = ['Protein', 'Session', 'Input_Id']

    results = select(config, columns, table, equals_rows={'Session': session_id})
    db, cursor = config.getDB()

    prot_input_dict = {}

    for row in results:
        if row[1] != session_id:
            continue
        prot_id = row[0]
        if prot_id in prot_id_list:
            prot_input_dict[prot_id] = row[2]

    if not includeSequence:
        rows = ['Protein_Id', 'Primary_Protein_Id', 'Uniprot_Ac', 'RefSeq_Ids', 'Uniprot_Id', 'Error_Code', 'Error']
    else:
        rows = ['Protein_Id', 'Primary_Protein_Id', 'Uniprot_Ac', 'RefSeq_Ids', 'Uniprot_Id', 'Error_Code', 'Error', 'Sequence']
    table = 'Protein'
    results = binningSelect(prot_id_list, rows, table, config)

    protein_dict = {}
    for row in results:
        prot_id = row[0]
        if prot_id in prot_input_dict:
            if not includeSequence:
                protein_dict[prot_id] = (row[1], row[2], row[3], row[4], row[5], row[6], prot_input_dict[prot_id])
            else:
                protein_dict[prot_id] = (row[1], row[2], row[3], row[4], row[5], row[6], prot_input_dict[prot_id], row[7])

    return protein_dict


# called from serializedPipeline
def protCheck(proteins, session_id, config):
    # Scan the database for stored Proteins
    if proteins.isEmpty():
        return {}, {}

    results = select(config, ['Protein_Id', 'Primary_Protein_Id', 'Sequence'], 'Protein')

    if config.verbosity >= 3:
        print('Just after protCheck selection')

    prot_id_list = set([])
    prot_ids_mutants_excluded = set()

    max_database_id = 0

    for row in results:
        database_id = row[0]
        prot_id = row[1]
        if not proteins.contains(prot_id):
            continue
        proteins.set_protein_stored(prot_id, True)
        proteins.set_protein_db_id(prot_id, database_id)
        proteins.set_protein_sequence(prot_id, row[2])

        prot_id_list.add(database_id)
        if not sdsc.is_mutant_ac:
            prot_ids_mutants_excluded.add(database_id)

        if database_id > max_database_id:
            max_database_id = database_id

    proteins.set_stored_ids(prot_id_list, prot_ids_mutants_excluded)

    # Insert the new proteins into the database
    new_prots = set()
    prot_ids = proteins.get_protein_ids()
    for prot_id in prot_ids:
        if not proteins.is_protein_stored(prot_id):
            new_prots.add(prot_id)

    if len(new_prots) > 0:
        values = []
        for prot_id in new_prots:
            u_ac = proteins.get_u_ac(prot_id)
            ref_id = proteins.get_ref_id(prot_id)
            u_id = proteins.get_u_id(prot_id)
            values.append((prot_id, u_ac, ref_id, u_id, session_id))

        insert('Protein', ['Primary_Protein_Id', 'Uniprot_Ac', 'RefSeq_Ids', 'Uniprot_Id', 'Original_Session'], values, config)

        # Retrieve the Protein-Ids from the new added proteins

        db, cursor = config.getDB()

        sql = "SELECT Protein_Id,Primary_Protein_Id FROM Protein WHERE Protein_Id > %s" % str(max_database_id)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e, f, g] = sys.exc_info()
            raise NameError("Error in protCheck: %s,\n%s" % (sql, f))
        db.close()
    else:
        results = ()

    new_ids = set()
    for row in results:
        database_id = row[0]
        prot_id = row[1]
        if not proteins.contains(prot_id):
            continue
        proteins.set_protein_db_id(prot_id, database_id)

        if database_id not in prot_id_list:
            new_ids.add(database_id)

    proteins.set_not_stored_ids(new_ids)

    proteins.generate_id_map()

    # Insert the Protein-Session-Connections into the database
    values = []
    for prot_id in prot_ids:
        database_id = proteins.get_protein_database_id(prot_id)
        input_id = proteins[prot_id].input_id
        values.append((database_id, session_id, input_id))

    insert('RS_Protein_Session', ['Protein', 'Session', 'Input_Id'], values, config)


def getLastAutoInc(db, cursor):
    sql = "SELECT LAST_INSERT_ID()"

    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        for row in results:
            auto_id = row[0]
        db.commit()
    except:
        print("NameError in getLastAutoInc")
        db.rollback()
    return auto_id


# called by babel
def getUniprotAcFromId(Protein_Id, db, cursor):
    sql = "SELECT Uniprot_Ac FROM Protein WHERE Protein_Id = '%s'" % (str(Protein_Id))

    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("NameError in getUniprotFromId")
        db.rollback()

    if results == ():
        return 0
    row = results[0]
    return row[0]


# called by babel
def getAAChange(mutation_id, db, cursor):
    sql = "SELECT Amino_Acid_Change FROM Position WHERE Position_Id = '%s'" % str(mutation_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getAAChange: %s" % sql)
        db.rollback()

    if results == ():
        return ""
    for row in results:
        aac = row[0]
        return aac
    return ""


def insertMultiMutations(proteins, session, config):
    db_ids = []
    for wt_prot in proteins.multi_mutations:
        db_ids.append(proteins.get_protein_database_id(wt_prot))

    if len(db_ids) > 0:
        columns = ['Wildtype_Protein', 'Multi_Mutation_Id', 'SNVs', 'Indels']
        table = 'Multi_Mutation'

        results = binningSelect(db_ids, columns, table, config)

        for row in results:
            wt_database_id = row[0]
            wt_prot_id = proteins.getU_acByDbId(wt_database_id)
            if wt_prot_id not in proteins.multi_mutations:
                continue
            snv_db_ids = row[2].split(',')
            indel_db_ids = row[3].split(',')

            for multi_mutation in proteins.multi_mutations[wt_prot]:
                if snv_db_ids != multi_mutation.get_snv_db_ids():
                    continue

                if indel_db_ids != multi_mutation.get_indel_db_ids():
                    continue
                multi_mutation.database_id = row[1]
                multi_mutation.stored = True

    wt_prot_db_ids = []
    values = []
    for wt_prot in proteins.multi_mutations:
        wt_db_id = proteins.get_protein_database_id(wt_prot)
        wt_prot_db_ids.append(wt_db_id)
        for multi_mutation in proteins.multi_mutations[wt_prot]:
            if multi_mutation.stored:
                continue
            mut_db_id = proteins[multi_mutation.mut_prot].database_id
            snv_db_ids = ','.join([str(x) for x in multi_mutation.get_snv_db_ids()])
            indel_db_ids = ','.join([str(x) for x in multi_mutation.get_indel_db_ids()])
            values.append((wt_db_id, mut_db_id, snv_db_ids, indel_db_ids))

    if len(values) > 0:
        columns = ['Wildtype_Protein', 'Mutant_Protein', 'SNVs', 'Indels']
        table = 'Multi_Mutation'
        insert(table, columns, values, config)

    columns = ['Wildtype_Protein', 'Mutant_Protein', 'Multi_Mutation_Id']
    table = 'Multi_Mutation'

    results = binningSelect(wt_prot_db_ids, columns, table, config)
    for row in results:
        wt_db_id = row[0]
        wt_prot_id = proteins.getU_acByDbId(wt_db_id)
        mut_db_id = row[1]
        mut_prot_id = proteins.getU_acByDbId(mut_db_id)
        for multi_mutation in proteins.multi_mutations[wt_prot_id]:
            if multi_mutation.mut_prot != mut_prot_id:
                continue
            multi_mutation.database_id = row[2]

    columns = ['Multi_Mutation', 'Session', 'Tags']
    table = 'RS_Multi_Mutation_Session'

    values = []
    for wt_prot in proteins.multi_mutations:
        for multi_mutation in proteins.multi_mutations[wt_prot]:
            values.append((multi_mutation.database_id, session, ','.join(multi_mutation.tags)))

    if len(values) > 0:
        insert(table, columns, values, config)

    return


# called by serializedPipeline
def indelCheck(proteins, session, config):
    if len(proteins.indels) == 0:
        return
    stored_wt_ids = proteins.get_stored_ids(exclude_indel_mutants=True)

    if len(stored_wt_ids) > 0:
        rows = ['Wildtype_Protein', 'Indel_Notation', 'Indel_Id']
        table = 'Indel'

        results = binningSelect(stored_wt_ids, rows, table, config)

        for row in results:
            wt_prot_id = row[0]
            indel_notation = row[1]
            database_id = row[2]
            u_ac = proteins.getU_acByDbId(wt_prot_id)
            proteins.indels[u_ac][indel_notation].set_database_id(database_id)
            proteins.indels[u_ac][indel_notation].set_stored(True)

    values = []
    all_wt_ids = set()
    for u_ac in proteins.indels:
        for indel_notation in proteins.indels[u_ac]:
            if proteins.indels[u_ac][indel_notation].stored:
                continue
            wt_prot_id = proteins.get_protein_database_id(proteins.indels[u_ac][indel_notation].wt_prot)
            mut_prot_id = proteins.get_protein_database_id(proteins.indels[u_ac][indel_notation].mut_prot)
            if wt_prot_id is None:
                config.errorlog.add_error('Wildtype protein id is not allowed to be None')
                continue
            values.append((wt_prot_id, mut_prot_id, indel_notation))
            all_wt_ids.add(wt_prot_id)

    if len(values) > 0:
        columns = ['Wildtype_Protein', 'Mutant_Protein', 'Indel_Notation']
        insert('Indel', columns, values, config)

    session_values = []

    if len(all_wt_ids) > 0:
        rows = ['Wildtype_Protein', 'Indel_Notation', 'Indel_Id']
        table = 'Indel'
        results = binningSelect(all_wt_ids, rows, table, config)
        for row in results:
            indel_notation = row[1]
            wt_prot_id = row[0]
            u_ac = proteins.getU_acByDbId(wt_prot_id)
            if u_ac not in proteins.indels:
                continue
            if indel_notation not in proteins.indels[u_ac]:
                continue
            database_id = row[2]
            proteins.indels[u_ac][indel_notation].set_database_id(database_id)
            session_values.append((session, database_id, ','.join(proteins.indels[u_ac][indel_notation].tags)))

    if len(session_values) > 0:
        columns = ['Session', 'Indel', 'Tags']
        insert('RS_Indel_Session', columns, session_values, config)


# called by serializedPipeline
def positionCheck(proteins, database_session, config):
    verbose = config.verbose

    # search for stored positions
    stored_positions = []
    pos_map = {}
    if len(proteins.get_stored_ids()) > 0:
        columns = ['Protein', 'Position_Number', 'Position_Id', 'Recommended_Structure_Data']
        table = 'Position'

        results = binningSelect(proteins.get_stored_ids(), columns, table, config)

        for row in results:
            prot_id = row[0]

            pos = row[1]
            position_database_id = row[2]
            if not proteins.position_in_protein_by_db_id(prot_id, pos):
                continue

            proteins.set_position_stored(prot_id, pos, True)
            proteins.set_position_database_id(prot_id, pos, position_database_id)
            stored_positions.append(position_database_id)
            pos_map[position_database_id] = (prot_id, pos)
            if row[3] is not None:
                recommended_structure_tuple = unpack(row[3])
                recommended_structure, seq_id, cov, resolution = sdsc.process_recommend_structure_str(recommended_structure_tuple[0])
                proteins.getByDbId(prot_id).positions[pos].recommended_structure = recommended_structure

    # search for stored SNVs
    if len(stored_positions) > 0:
        columns = ['Position', 'New_AA', 'SNV_Id']
        table = 'SNV'
        results = binningSelect(stored_positions, columns, table, config)

        for row in results:
            prot_database_id, pos = pos_map[row[0]]
            prot_id = proteins.getU_acByDbId(prot_database_id)
            if not row[1] in proteins[prot_id].positions[pos].mut_aas:
                continue
            proteins[prot_id].positions[pos].mut_aas[row[1]].database_id = row[2]
            proteins[prot_id].positions[pos].mut_aas[row[1]].stored = True

    # insert new positions
    prot_ids = proteins.get_protein_ids()
    values = []
    for prot_id in prot_ids:
        prot_database_id = proteins.get_protein_database_id(prot_id)
        positions = proteins.get_position_ids(prot_id)
        all_stored = True
        for pos in positions:
            if proteins.is_position_stored(prot_id, pos):
                continue
            all_stored = False
            aac_base = proteins.get_aac_base(prot_id, pos)

            res_id = proteins.get_res_id(prot_id, pos)

            values.append((prot_database_id, int(aac_base[1:]), res_id, aac_base[0]))

        if all_stored:
            proteins.set_completely_stored(prot_id)

    if len(values) > 0:
        columns = ['Protein', 'Position_Number', 'Residue_Id', 'Wildtype_Residue']
        insert('Position', columns, values, config)

    # retrieve the database ids of the new positions
    columns = ['Protein', 'Position_Number', 'Position_Id']
    table = 'Position'

    fused_prot_ids = proteins.get_not_stored_ids() | proteins.get_stored_ids()

    results = binningSelect(fused_prot_ids, columns, table, config)

    position_database_ids = []
    prot_back_map = {}
    for row in results:
        prot_database_id = row[0]

        pos = row[1]

        position_database_id = row[2]
        if not proteins.position_in_protein_by_db_id(prot_database_id, pos):
            continue
        proteins.set_position_database_id(prot_database_id, pos, position_database_id)
        position_database_ids.append(position_database_id)
        prot_back_map[position_database_id] = (prot_database_id, pos)

    # insert the new SNVs
    values = []
    for prot_id in prot_ids:
        for pos in proteins[prot_id].positions:
            pos_database_id = proteins[prot_id].positions[pos].database_id
            for new_aa in proteins[prot_id].positions[pos].mut_aas:
                values.append((pos_database_id, new_aa))
    if len(values) > 0:
        table = 'SNV'
        columns = ['Position', 'New_AA']
        insert(table, columns, values, config)

    # retrieve the database ids of the new SNVs
    columns = ['Position', 'New_AA', 'SNV_Id']
    table = 'SNV'
    results = binningSelect(position_database_ids, columns, table, config)

    for row in results:
        position_database_id = row[0]
        if position_database_id not in prot_back_map:
            continue
        (prot_database_id, pos) = prot_back_map[position_database_id]
        prot_id = proteins.getU_acByDbId(prot_database_id)
        new_aa = row[1]
        if new_aa not in proteins[prot_id].positions[pos].mut_aas:
            continue
        proteins[prot_id].positions[pos].mut_aas[new_aa].database_id = row[2]

    # insert the the SNV session connections
    values = []
    for prot_id in prot_ids:
        for pos in proteins[prot_id].positions:
            for new_aa in proteins[prot_id].positions[pos].mut_aas:
                tags = ','.join(proteins[prot_id].positions[pos].mut_aas[new_aa].tags)
                snv_database_id = proteins[prot_id].positions[pos].mut_aas[new_aa].database_id
                values.append((database_session, snv_database_id, tags))

    if len(values) > 0:
        table = 'RS_SNV_Session'
        columns = ['Session', 'SNV', 'Tag']
        insert(table, columns, values, config)

    # insert the position session connections
    values = []
    for prot_id in prot_ids:
        positions = proteins.get_position_ids(prot_id)
        for pos in positions:
            pos_id = proteins.get_position_database_id(prot_id, pos)
            pos_tags = proteins.get_pos_tags(prot_id, pos)
            values.append((database_session, pos_id, ','.join(pos_tags)))

    if len(values) > 10000:
        process = multiprocessing.Process(target=backgroundInsertMS, args=(values, config))
        try:
            process.start()
        except:
            [e, f, g] = sys.exc_info()
            raise NameError("Error in PositionCheck: %s" % (f))
        return process
    else:
        insert('RS_Position_Session', ['Session', 'Position', 'Tag'], values, config)
        return None


def backgroundInsertMS(values, config):
    insert('RS_Position_Session', ['Session', 'Position', 'Tag'], values, config)


# called by serializedPipeline
def addIupred(proteins, config):
    values = []
    u_acs = proteins.get_protein_ids()
    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        scores = proteins.get_disorder_scores(u_ac)
        regions = proteins.get_disorder_regions(u_ac)
        method = proteins.get_disorder_tool(u_ac)

        positions = proteins.get_position_ids(u_ac)

        for pos in positions:
            pos_id = proteins.get_position_database_id(u_ac, pos)

            if method == 'MobiDB3.0' or method == 'mobidb-lite':
                pos_region_type = 'globular'
            else:
                pos_region_type = 'disorder'
            if regions is None:
                continue
            for [a, b, region_type] in regions:
                if int(pos) > int(a) and int(pos) < int(b):
                    pos_region_type = region_type

            if scores is None:
                config.errorlog.add_error('IUpred scores are None for: %s' % u_ac)
                break

            if pos not in scores:
                continue
            iupred_score = scores[pos]
            values.append((pos_id, iupred_score, pos_region_type))

    process = multiprocessing.Process(target=backgroundIU, args=(values, config))
    process.start()
    return process


def backgroundIU(values, config):
    if not len(values) == 0:
        update(config, 'Position', ['Position_Id', 'IUPRED', 'IUPRED_Glob'], values)
    return


# called by babel
def getPDB(template_id, db, cursor):
    sql = "SELECT Name FROM Template WHERE Template_Id = '%s'" % str(template_id)

    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        print("NameError in getPDB")
        db.rollback()

    return results[0][0]


# called by babel
def checkMutationSession(mutation_id, session_id, db, cursor):
    sql = "SELECT Position FROM RS_Position_Session WHERE Position = '%s' AND Session = '%s'" % (str(mutation_id), str(session_id))
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in checkMutationSession")
        db.rollback()
    if results == ():
        return False
    else:
        return True


# called by output.classification
def getMutationDict(mutation_id_list, config):
    if len(mutation_id_list) == 0:
        return {}

    rows = ['Position_Id', 'Position_Number', 'Protein', 'IUPRED', 'IUPRED_Glob', 'Residue_Id', 'Wildtype_Residue']
    table = 'Position'
    results = binningSelect(mutation_id_list, rows, table, config)

    mutation_dict = {}

    for row in results:
        mut_id = row[0]
        if mut_id in mutation_id_list:
            mutation_dict[mut_id] = (row[1], row[2], row[3], row[4], row[5], row[6])
    return mutation_dict


# called by serializedPipeline
def updateSession(session_id, time, config):
    db, cursor = config.getDB()
    sql = "UPDATE IGNORE Session SET End = '%s' WHERE Session_Id = '%s'" % (str(time), str(session_id))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Commit your changes in the database
        db.commit()
    except:
        print(("Couldn't update Session '%s'") % (str(session_id)))
        # Rollback in case there is any NameError
        db.rollback()
    db.close()


# called by serializedPipeline
def insertSession(time, ori_filename, config):
    db, cursor = config.getDB()
    checksum = calc_checksum(ori_filename)
    sql = "INSERT IGNORE INTO Session(Input_File,Checksum,Start) VALUES ('%s', %d,'%s');" % (ori_filename, checksum, str(time))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Couldn't insert Session '%s'") % (sql)
        # Rollback in case there is any NameError
        db.rollback()
    session_id = getLastAutoInc(db, cursor)
    db.close()
    return session_id


# called by output
# called by structman
def getSessionId(infile, config):
    db, cursor = config.getDB()
    checksum = calc_checksum(infile)
    sql = "SELECT Session_Id FROM Session WHERE Input_File = '%s' AND Checksum = %d;" % (infile, checksum)

    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getSessionId: %s" % sql)
        # Rollback in case there is any NameError
    db.close()

    if results == ():
        return None
    try:
        sid = results[0][0]
    except:
        raise NameError("Input_File not in Session: %s" % sql)
    return results[0][0]


# called by babel
# called by output
def getSessionFile(session_id, db, cursor):
    sql = "SELECT Input_File FROM Session WHERE Session_Id = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getSessionFile: '%s'" % sql)
        # Rollback in case there is any NameError

    if results == ():
        return None
    return results[0][0]


# called by serializedPipeline
def addProtInfos(proteins, config):
    ref_value_strs = []
    go_term_ids = set()
    reac_ids = set()
    seq_values = []
    u_acs = proteins.get_protein_ids()

    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        go_terms = proteins.get_go_terms(u_ac)
        pathways = proteins.get_pathways(u_ac)
        for go_id in go_terms:
            go_term_ids.add(go_id)
        for reac_id in pathways:
            reac_ids.add(reac_id)
        prot_id = proteins.get_protein_database_id(u_ac)
        seq = proteins.get_sequence(u_ac)
        seq_values.append((prot_id, seq))

    stored_go_terms = {}
    if len(go_term_ids) > 0:
        results = select(config, ['Id', 'GO_Term_Id'], 'GO_Term', in_rows={'ID': go_term_ids})

        for row in results:
            stored_go_terms[row[0]] = row[1]

    stored_pathways = {}
    if len(reac_ids) > 0:
        results = select(config, ['Reactome_Id', 'Pathway_Id'], 'Pathway', in_rows={'Reactome_Id': reac_ids})

        for row in results:
            stored_pathways[row[0]] = row[1]

    new_go_terms = {}
    go_values = []
    new_pathways = {}
    pathway_values = []
    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        go_terms = proteins.get_go_terms(u_ac)
        pathways = proteins.get_pathways(u_ac)
        for go_id in go_terms:
            if go_id not in stored_go_terms:
                new_go_terms[go_id] = go_terms[go_id]
                go_values.append((go_terms[go_id].replace("'", "''"), go_id))
                stored_go_terms[go_id] = None
        for reac_id in pathways:
            if reac_id not in stored_pathways:
                new_pathways[reac_id] = pathways[reac_id]
                pathway_values.append((pathways[reac_id].replace("'", "''"), reac_id))
                stored_pathways[reac_id] = None

    if len(go_values) > 0:
        insert('GO_Term', ['Name', 'Id'], go_values, config)

    if len(pathway_values) > 0:
        insert('Pathway', ['Name', 'Reactome_Id'], pathway_values, config)

    stored_go_terms = {}
    if len(go_term_ids) > 0:
        results = select(config, ['Id', 'GO_Term_Id', ], 'GO_Term', in_rows={'Id': go_term_ids})

        for row in results:
            stored_go_terms[row[0]] = row[1]

    stored_pathways = {}
    if len(reac_ids) > 0:
        results = select(config, ['Reactome_Id', 'Pathway_Id'], 'Pathway', in_rows={'Reactome_Id': reac_ids})

        for row in results:
            stored_pathways[row[0]] = row[1]

    go_values = []
    pathway_values = []
    for u_ac in u_acs:
        if proteins.is_protein_stored(u_ac):
            continue
        prot_id = proteins.get_protein_database_id(u_ac)
        go_terms = proteins.get_go_terms(u_ac)
        pathways = proteins.get_pathways(u_ac)
        for go_id in go_terms:
            go_values.append((prot_id, stored_go_terms[go_id]))
        for reac_id in pathways:
            pathway_values.append((prot_id, stored_pathways[reac_id]))

    if len(go_values) > 0:
        insert('RS_Protein_GO_Term', ['Protein', 'GO_Term'], go_values, config)

    if len(pathway_values) > 0:
        insert('RS_Protein_Pathway', ['Protein', 'Pathway'], pathway_values, config)

    if len(seq_values) > 0:
        update(config, 'Protein', ['Protein_Id', 'Sequence'], seq_values)


# called by serializedPipeline
def insertComplexes(proteins, config):

    smiles_path = config.smiles_path
    inchi_path = config.inchi_path
    pdb_path = config.pdb_path

    stored_ligands = {}
    results = select(config, ['Ligand_Id', 'Name'], 'Ligand')
    for row in results:
        ligand_name = row[1]
        stored_ligands[ligand_name] = row[0]

    values = []

    ligand_db = pdbParser.parseLigandDB(smiles_path, inchi_path)

    lig_values = []
    update_ligands = []
    new_ligands = set()
    complex_list = proteins.get_complex_list()

    for pdb_id in complex_list:
        if proteins.is_complex_stored(pdb_id):
            continue

        resolution = proteins.get_resolution(pdb_id)
        chains_str = proteins.get_chains_str(pdb_id)
        homomers_str = proteins.get_homomers_str(pdb_id)
        lig_str = proteins.get_lig_str(pdb_id)
        metal_str = proteins.get_metal_str(pdb_id)
        ion_str = proteins.get_ion_str(pdb_id)
        cc_str = proteins.get_chain_chain_str(pdb_id)

        values.append((pdb_id, resolution, chains_str, homomers_str, lig_str, metal_str, ion_str, cc_str))

        interaction_partners = proteins.get_interaction_partners(pdb_id)

        for iap in interaction_partners:
            ia_type = iap[0]
            if ia_type != "Ligand":
                continue
            name = iap[1]
            if name in stored_ligands:
                continue
            if name == "UNK" or name == "UNX":
                continue
            if name in new_ligands:
                continue
            res = iap[2]
            chain = iap[3]
            if name in ligand_db:
                (smiles, inchi) = ligand_db[name]
            elif len(name) < 3:
                if len(name) == 1:
                    smiles = "[%s]" % name
                else:
                    smiles = "[%s%s]" % (name[0], name[1].lower())
                inchi = "InChI=1S/%s" % smiles
            else:
                (smiles, inchi) = pdbParser.getSI(pdb_id, name, res, chain, pdb_path, config)
                update_ligands.append((name, smiles, inchi))
            new_ligands.add(name)
            lig_values.append((name, smiles, inchi))

    if len(update_ligands) > 0:
        pdbParser.updateLigandDB(update_ligands, smiles_path, inchi_path)

    if len(lig_values) > 0:
        try:
            insert('Ligand', ['Name', 'Smiles', 'Inchi'], lig_values, config)
        except:
            # There is a (or more than one) ligand that results in an sql error 1241. We need to find it, to figure out why.
            config.errorlog.add_error('Error trying inserting ligands\n%s' % '\n'.join([str(x) for x in lig_values]))

    if len(values) > 0:
        insert('Complex', ['PDB', 'Resolution', 'Chains', 'Homooligomers', 'Ligand_Profile', 'Metal_Profile', 'Ion_Profile', 'Chain_Chain_Profile'], values, config)

    stored_ligands = {}
    results = select(config, ['Ligand_Id', 'Name'], 'Ligand')
    for row in results:
        ligand_name = row[1]
        stored_ligands[ligand_name] = row[0]

    values = []
    results = select(config, ['Complex_Id', 'PDB'], 'Complex')
    for row in results:
        if not proteins.contains_complex(row[1]):
            continue
        if proteins.is_complex_stored(row[1]):
            continue
        proteins.set_complex_db_id(row[1], row[0])

        interaction_partners = proteins.get_interaction_partners(row[1])

        for iap in interaction_partners:
            ia_type = iap[0]
            if ia_type != "Ligand":
                continue
            name = iap[1]
            if name == "UNK" or name == "UNX" or name not in stored_ligands:
                continue
            lig_id = stored_ligands[name]
            res = iap[2]
            chain = iap[3]
            values.append((lig_id, row[0], chain, res))

    if len(values) > 0:
        insert('RS_Ligand_Structure', ['Ligand', 'Complex', 'Chain', 'Residue'], values, config)


def getComplexMap(config, pdb_ids=None):
    table = 'Complex'
    rows = ['Complex_Id', 'PDB', 'Resolution', 'Chains', 'Homooligomers', 'Ligand_Profile', 'Metal_Profile', 'Ion_Profile', 'Chain_Chain_Profile']

    results = select(config, rows, table)

    complex_map = {}
    for row in results:
        pdb_id = row[1]
        if pdb_ids is not None:
            if pdb_id not in pdb_ids:
                continue
        comp_id = row[0]
        resolution = row[2]
        chains_str = row[3]
        homooligomers = row[4]
        lig_profile_str = row[5]
        metal_profile_str = row[6]
        ion_profile_str = row[7]
        cc_profile_str = row[8]

        lig_profile = {}
        if lig_profile_str != '':
            for lig_profile_part in lig_profile_str.split(','):
                part1, part2 = lig_profile_part.split(':')
                chain, res = part1.split('_')
                deg, score = part2.split('_')

                lig_profile[(chain, res)] = int(deg), float(score)

        metal_profile = {}

        if metal_profile_str != '':
            for metal_profile_part in metal_profile_str.split(','):
                part1, part2 = metal_profile_part.split(':')
                chain, res = part1.split('_')
                deg, score = part2.split('_')

                metal_profile[(chain, res)] = int(deg), float(score)

        ion_profile = {}

        if ion_profile_str != '':
            for ion_profile_part in ion_profile_str.split(','):
                part1, part2 = ion_profile_part.split(':')
                chain, res = part1.split('_')
                deg, score = part2.split('_')

                ion_profile[(chain, res)] = int(deg), float(score)

        cc_profile = {}

        if cc_profile_str != '':
            for cc_profile_part in cc_profile_str.split(','):
                part1, part2 = cc_profile_part.split(':')
                chain, chain_b = part1.split('_')
                deg, score = part2.split('_')

                cc_profile[(chain, chain_b)] = int(deg), float(score)

        complex_map[pdb_id] = (comp_id, resolution, chains_str, homooligomers, lig_profile, metal_profile, ion_profile, cc_profile)

    return complex_map


# called by serializedPipeline
def structureCheck(proteins, config):
    table = 'Structure'
    rows = ['Structure_Id', 'PDB', 'Chain', 'Homooligomer']
    results = select(config, rows, table)

    stored_complexes = set()

    for row in results:
        s_id = row[0]
        pdb_id = row[1]
        chain = row[2]
        oligos = row[3]
        if not proteins.contains_structure(pdb_id, chain):
            continue
        proteins.set_structure_db_id(pdb_id, chain, s_id)
        proteins.set_structure_stored(pdb_id, chain, True)  # all structures, mapped or not go into this dictionary, this is important for not reinserting residues from interacting structures


def draw_complexes(config, proteins, stored_complexes=[], draw_all=False):
    results = select(config, ['Complex_Id', 'PDB', 'Resolution', 'Chains', 'Ligand_Profile', 'Metal_Profile', 'Ion_Profile', 'Chain_Chain_Profile', 'Homooligomers'], 'Complex')

    for row in results:
        pdb_id = row[1]
        if pdb_id not in stored_complexes and not draw_all:
            continue

        compl = sdsc.Complex(pdb_id, resolution=float(row[2]), chains_str=row[3], lig_profile_str=row[4], metal_profile_str=row[5], ion_profile_str=row[6], chain_chain_profile_str=row[7], stored=True, database_id=row[0], homomers_str=row[8])

        proteins.add_complex(pdb_id, compl)


# called by serializedPipeline
def insertStructures(structurelist, proteins, config, results=None, return_results=False):

    table = 'Structure'
    rows = ['Structure_Id', 'PDB', 'Chain']
    if results is None:
        results = select(config, rows, table)
        already_called = False
    else:
        already_called = True

    stored_complexes = set()
    del_list = []

    for pos, row in enumerate(results):
        pdb_id = row[1]
        chain = row[2]
        if not proteins.contains_structure(pdb_id, chain):
            if return_results:
                del_list.append(pos)
            continue
        if not already_called:
            s_id = row[0]
            proteins.set_structure_db_id(pdb_id, chain, s_id)
            proteins.set_structure_stored(pdb_id, chain, True)  # all structures, mapped or not go into this dictionary, this is important for not reinserting residues from interacting structures
            stored_complexes.add(pdb_id)
        if (pdb_id, chain) in structurelist:
            structurelist.remove((pdb_id, chain))

    if return_results:
        results = list(results)
        for pos in reversed(del_list):
            del results[pos]
        ret_results = results
    else:
        ret_results = None

    if not already_called:
        draw_complexes(config, proteins, stored_complexes=stored_complexes)

    values = []

    for (pdb_id, chain) in structurelist:
        oligos = proteins.get_oligo(pdb_id, chain)
        oligos = ''.join(oligos)
        values.append((pdb_id, chain, oligos))

    if len(values) > 0:
        insert('Structure', ['PDB', 'Chain', 'Homooligomer'], values, config)

    if len(structurelist) > 0:
        table = 'Structure'
        rows = ['Structure_Id', 'PDB', 'Chain']
        results = select(config, rows, table)

        for row in results:
            s_id = row[0]
            pdb_id = row[1]
            chain = row[2]

            if not (pdb_id, chain) in structurelist:
                continue

            if return_results:
                ret_results.append(row)

            proteins.set_structure_db_id(pdb_id, chain, s_id)

    return ret_results


# called by serializedPipeline
def insertInteractingChains(interaction_structures, proteins, config):
    if len(interaction_structures) == 0:

        return {}

    interacting_structure_ids = {}

    results = select(config, ['Structure_Id', 'PDB', 'Chain'], 'Structure')

    for row in results:
        pdb_id = row[1]
        chain = row[2]
        if not (pdb_id, chain) in interaction_structures:
            continue
        #s_id = row[0]
        #interacting_structure_ids[(pdb_id,chain)] = s_id
        interaction_structures.remove((pdb_id, chain))

    values = []

    for (pdb_id, chain) in interaction_structures:

        homomers = proteins.get_complex_homomers(pdb_id, chain)

        oligos = ''.join(homomers)

        values.append((pdb_id, chain, oligos))

    if len(values) > 0:
        insert('Structure', ['PDB', 'Chain', 'Homooligomer'], values, config)

    if len(interaction_structures) > 0:
        results = select(config, ['Structure_Id', 'PDB', 'Chain'], 'Structure')

        for row in results:
            pdb_id = row[1]
            chain = row[2]
            if not (pdb_id, chain) in interaction_structures:
                continue
            s_id = row[0]

            interacting_structure_ids[(pdb_id, chain)] = s_id

    return interacting_structure_ids


# called by serializedPipeline
def insertAlignments(alignment_list, proteins, config):
    values = []
    if config.verbosity >= 2:
        t0 = time.time()
    for (u_ac, prot_id, pdb_id, chain, alignment_pir) in alignment_list:
        s_id = proteins.get_structure_db_id(pdb_id, chain)
        seq_id = proteins.get_sequence_id(u_ac, pdb_id, chain)
        coverage = proteins.get_coverage(u_ac, pdb_id, chain)
        values.append((prot_id, s_id, seq_id, coverage, pack(alignment_pir)))
    if config.verbosity >= 2:
        t1 = time.time()
        print('Time for insertAlignments, part 1: ', t1 - t0)
    if len(values) > 0:
        insert('Alignment', ['Protein', 'Structure', 'Sequence_Identity', 'Coverage', 'Alignment'], values, config)
    if config.verbosity >= 2:
        t2 = time.time()
        print('Time for insertAlignments, part 2: ', t2 - t1)


def background_insert_residues(values, config):
    """
    columns = ['Structure', 'Number', 'Amino_Acid', 'Sub_Lig_Dist', 'Sub_Chain_Distances', 'Relative_Surface_Access',
               'Relative_Surface_Access_Main_Chain', 'Relative_Surface_Access_Side_Chain', 'Secondary_Structure_Assignment',
               'Homomer_Distances', 'Interaction_Profile', 'Centralities', 'B_Factor', 'Modres', 'PHI', 'PSI', 'Intra_SSBOND', 'SSBOND_Length',
               'Intra_Link', 'Link_Length', 'CIS_Conformation', 'CIS_Follower', 'Inter_Chain_Median_KD', 'Inter_Chain_Dist_Weighted_KD', 'Inter_Chain_Median_RSA',
               'Inter_Chain_Dist_Weighted_RSA', 'Intra_Chain_Median_KD', 'Intra_Chain_Dist_Weighted_KD', 'Intra_Chain_Median_RSA', 'Intra_Chain_Dist_Weighted_RSA',
               'Inter_Chain_Interactions_Median', 'Inter_Chain_Interactions_Dist_Weighted',
               'Intra_Chain_Interactions_Median', 'Intra_Chain_Interactions_Dist_Weighted',
               'Interacting_Chains', 'Interacting_Ligands'
               ]
    """
    columns = ['Structure', 'Residue_Data']
    insert('Residue', columns, values, config)


# called by serializedPipeline
def insertResidues(structural_analysis, interacting_structure_ids, proteins, config):

    if config.verbosity >= 2:
        t0 = time.time()

    values = []

    if len(structural_analysis) == 0:
        return

    structure_ids = {}

    t_0 = 0.
    t_1 = 0.
    t_2 = 0.
    t_3 = 0.
    t_4 = 0.
    t_5 = 0.
    t_6 = 0.
    t_7 = 0.
    t_8 = 0.
    t_9 = 0.
    t_10 = 0.
    t_11 = 0.
    t_12 = 0.
    t_13 = 0.

    for (pdb_id, chain) in structural_analysis:
        if not (proteins.contains_structure(pdb_id, chain) or (pdb_id, chain) in interacting_structure_ids):
            continue
        if proteins.contains_structure(pdb_id, chain) and proteins.is_structure_stored(pdb_id, chain):  # all residues belonging to stored structures must not inserted twice
            continue

        analysis_map = structural_analysis[(pdb_id, chain)]
        if (pdb_id, chain) in interacting_structure_ids:
            s_id = interacting_structure_ids[(pdb_id, chain)]
        else:
            s_id = proteins.get_structure_db_id(pdb_id, chain)
            if s_id is None:
                continue
            structure_ids[s_id] = (pdb_id, chain)
        for res_id in analysis_map:
            t_0 += time.time()
            residue = analysis_map[res_id]

            one_letter = residue.get_aa()
            lig_dist_str = residue.get_lig_dist_str()
            chain_dist_str = residue.get_chain_dist_str()

            t_1 += time.time()

            (rsa, relative_main_chain_acc, relative_side_chain_acc) = residue.get_rsa(splitted=True)
            ssa = residue.get_ssa()

            t_2 += time.time()

            homo_str = residue.get_homo_dist_str()

            t_3 += time.time()

            profile_str = residue.get_interaction_profile_str()

            t_4 += time.time()

            interacting_chains_str = residue.get_interacting_chains_str()
            interacting_ligands_str = residue.get_interacting_ligands_str()

            t_5 += time.time()

            #if interacting_chains is None:
            #    interacting_chains_str = None
            #else:
            #    interacting_chains_str = ','.join(interacting_chains)

            #if interacting_ligands is None:
            #    interacting_ligands_str = None
            #else:
            #    interacting_ligands_str = ','.join(interacting_ligands)

            t_6 += time.time()

            centrality_score_str = residue.get_centrality_str()

            t_7 += time.time()

            b_factor = residue.get_b_factor()
            modres = residue.get_modres()
            phi, psi = residue.get_angles()

            t_8 += time.time()

            intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower = residue.get_residue_link_information()

            t_9 += time.time()

            (inter_chain_median_kd, inter_chain_dist_weighted_kd, inter_chain_median_rsa, inter_chain_dist_weighted_rsa,
                intra_chain_median_kd, intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa) = residue.get_milieu()

            t_10 += time.time()

            (inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                intra_chain_interactions_median, intra_chain_interactions_dist_weighted) = residue.get_interface_milieu()

            t_11 += time.time()

            packed_res_info = pack((res_id, one_letter, lig_dist_str, chain_dist_str, rsa, relative_main_chain_acc, relative_side_chain_acc,
                           ssa, homo_str, profile_str,
                           centrality_score_str, b_factor, modres, phi, psi, intra_ssbond, ssbond_length, intra_link, link_length,
                           cis_conformation, cis_follower, inter_chain_median_kd, inter_chain_dist_weighted_kd,
                           inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                           intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                           inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                           intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                           interacting_chains_str, interacting_ligands_str))
            
            values.append([s_id, packed_res_info])
            t_12 += time.time()

            if not (pdb_id, chain) in interacting_structure_ids:
                proteins.add_residue(pdb_id, chain, res_id, residue)

            t_13 += time.time()

    if config.verbosity >= 2:
        t1 = time.time()
        print('Time for insertResidues part 1:', t1 - t0)
        print('Individual parts: 1 -', (t_1 - t_0), '2 -', (t_2 - t_1), '3 -', (t_3 - t_2), '4 -', (t_4 - t_3), '5 -', (t_5 - t_4), '6 -', (t_6 - t_5),
                '7 -', (t_7 - t_6), '8 -', (t_8 - t_7), '9 -', (t_9 - t_8), '10 -', (t_10 - t_9), '11 -', (t_11 - t_10), '12 -', (t_12 - t_11), '13 -', (t_13 - t_12))

    process = multiprocessing.Process(target=background_insert_residues, args=(values, config))
    process.start()

    return process


# called by serializedPipeline
def getAlignments(proteins, config, get_all_alignments=False):
    if config.verbosity >= 2:
        t0 = time.time()

    prot_db_ids = proteins.get_stored_ids(exclude_completely_stored=(not get_all_alignments))

    for prot_id in proteins.indels:
        if proteins[prot_id].stored:
            prot_db_ids.add(proteins[prot_id].database_id)

    rows = ['Protein', 'Structure', 'Sequence_Identity', 'Coverage', 'Alignment']
    table = 'Alignment'
    results = binningSelect(prot_db_ids, rows, table, config)

    if config.verbosity >= 2:
        t1 = time.time()
        print('Time for part 1 in getAlignments:', (t1 - t0), ', number of stored proteins', len(prot_db_ids), ', number of stored alignments:', len(results))

    prot_structure_alignment_map = {}
    structure_ids = set()

    for row in results:

        prot_db_id = row[0]

        structure_id = row[1]
        seq_id = row[2]
        coverage = row[3]
        alignment = unpack(row[4])

        structure_ids.add(structure_id)

        target_seq, template_seq = sdsc.process_alignment_data(alignment)

        if prot_db_id not in prot_structure_alignment_map:
            prot_structure_alignment_map[prot_db_id] = {}

        prot_structure_alignment_map[prot_db_id][structure_id] = (target_seq, template_seq, coverage, seq_id)

    if config.verbosity >= 2:
        t2 = time.time()
        print("Time for part 2 in getAlignments: %s" % (str(t2 - t1)))

    structure_map, id_structure_map = getStructure_map(structure_ids, config)
    pdb_ids = set()
    for (pdb_id, chain) in structure_map:
        pdb_ids.add(pdb_id)

    complex_map = getComplexMap(config, pdb_ids=pdb_ids)

    for prot_db_id in prot_structure_alignment_map:
        prot_id = proteins.getByDbId(prot_db_id).primary_protein_id
        for structure_id in prot_structure_alignment_map[prot_db_id]:
            (pdb_id, chain) = id_structure_map[structure_id]

            (target_seq, template_seq, coverage, seq_id) = prot_structure_alignment_map[prot_db_id][structure_id]

            struct_anno = sdsc.StructureAnnotation(prot_id, pdb_id, chain, alignment=(target_seq, template_seq), stored=True)
            proteins.add_annotation(prot_id, pdb_id, chain, struct_anno)

            if not proteins.contains_structure(pdb_id, chain):
                oligo = structure_map[(pdb_id, chain)][1]
                struct = sdsc.Structure(pdb_id, chain, oligo=oligo, mapped_proteins=[prot_id], database_id=structure_id)
                proteins.add_structure(pdb_id, chain, struct)
            else:
                proteins.add_mapping_to_structure(pdb_id, chain, prot_id)

            if not proteins.contains_complex(pdb_id) and (pdb_id in complex_map):
                (comp_id, resolution, chains_str, homooligomers, lig_profile, metal_profile, ion_profile, cc_profile) = complex_map[pdb_id]
                compl = sdsc.Complex(pdb_id, resolution=resolution, chains_str=chains_str, lig_profile=lig_profile,
                                     metal_profile=metal_profile, ion_profile=ion_profile,
                                     chain_chain_profile=cc_profile, stored=True, database_id=comp_id, homomers_str=homooligomers)
                proteins.add_complex(pdb_id, compl)

            proteins.set_structure_db_id(pdb_id, chain, structure_id)

            proteins.set_coverage_by_db_id(prot_db_id, pdb_id, chain, coverage)
            proteins.set_sequence_id_by_db_id(prot_db_id, pdb_id, chain, seq_id)
            proteins.set_annotation_db_id_by_db_id(prot_db_id, pdb_id, chain, True)

    if config.verbosity >= 2:
        t3 = time.time()
        print("Time for part 3 in getAlignments: %s" % (str(t3 - t2)))


def getStructure_map(structure_ids, config):
    structure_map = {}
    id_structure_map = {}
    if len(structure_ids) > 0:
        results = binningSelect(structure_ids, ['Structure_Id', 'PDB', 'Chain', 'Homooligomer'], 'Structure', config)

        for row in results:
            s_id = row[0]

            pdb_id = row[1]
            chain = row[2]
            oligo = row[3]
            structure_map[(pdb_id, chain)] = (s_id, oligo)
            id_structure_map[s_id] = (pdb_id, chain)
    return structure_map, id_structure_map


# called by indel_analysis
def insert_indel_results(proteins, config):
    table = 'Indel'
    columns = ['Indel_Id', 'Wildtype_Protein', 'Mutant_Protein', 'Indel_Notation', 'Analysis_Results']

    values = []
    for prot_id in proteins.indels:
        for indel_notation in proteins.indels[prot_id]:
            indel_obj = proteins.indels[prot_id][indel_notation]
            wt_prot_id = proteins.get_protein_database_id(proteins.indels[prot_id][indel_notation].wt_prot)
            mut_prot_id = proteins.get_protein_database_id(proteins.indels[prot_id][indel_notation].mut_prot)
            analysis_results = pack((indel_obj.size, indel_obj.delta_delta_classification, indel_obj.wt_aggregates, indel_obj.mut_aggregates,
                                        indel_obj.left_flank_wt_aggregates, indel_obj.left_flank_mut_aggregates,
                                        indel_obj.right_flank_wt_aggregates, indel_obj.right_flank_mut_aggregates))
            values.append((indel_obj.database_id, wt_prot_id, mut_prot_id, indel_notation, analysis_results))
    update(config, table, columns, values)
    return


# called by serializedPipeline
def insertClassifications(proteins, config):

    if config.verbosity >= 2:
        t1 = time.time()

    table = 'Position'
    columns = ['Position_Id', 'Recommended_Structure_Data', 'Position_Data']

    values = createClassValues(proteins, config)

    if config.verbosity >= 2:
        t2 = time.time()
        print('insertClassifications part 1: %s' % (t2 - t1), 'Update Classifications of', len(values), 'positions')

    update(config, table, columns, values)

    if config.verbosity >= 2:
        t3 = time.time()
        print('insertClassifications part 2: %s' % (t3 - t2))


def createClassValues(proteins, config):
    values = []

    for u_ac in proteins.get_protein_ids():
        if proteins.is_completely_stored(u_ac):
            continue
        positions = proteins.get_position_ids(u_ac)
        for pos in positions:
            aachange = proteins.get_aac_base(u_ac, pos)
            pos = int(aachange[1:])

            if proteins.is_position_stored(u_ac, pos):
                continue

            m = proteins.get_position_database_id(u_ac, pos)

            position = proteins.get_position(u_ac, pos)
            mappings = position.mappings

            values.append((m, pack((mappings.recommended_res, mappings.max_seq_res)),
                           pack((mappings.weighted_location, mappings.weighted_mainchain_location,
                           mappings.weighted_sidechain_location, mappings.weighted_surface_value,
                           mappings.weighted_mainchain_surface_value, mappings.weighted_sidechain_surface_value,
                           mappings.Class, mappings.rin_class, mappings.simple_class, mappings.rin_simple_class,
                           str(mappings.interaction_recommendations), mappings.classification_conf, mappings.weighted_ssa, len(mappings.qualities), mappings.get_weighted_profile_str(),
                           mappings.weighted_modres, mappings.weighted_b_factor, mappings.get_weighted_centralities_str(),
                           mappings.weighted_phi, mappings.weighted_psi, mappings.weighted_intra_ssbond, mappings.weighted_inter_ssbond,
                           mappings.weighted_intra_link, mappings.weighted_inter_link, mappings.weighted_cis_conformation,
                           mappings.weighted_cis_follower,
                           mappings.weighted_inter_chain_median_kd, mappings.weighted_inter_chain_dist_weighted_kd,
                           mappings.weighted_inter_chain_median_rsa, mappings.weighted_inter_chain_dist_weighted_rsa,

                           mappings.weighted_intra_chain_median_kd, mappings.weighted_intra_chain_dist_weighted_kd,
                           mappings.weighted_intra_chain_median_rsa, mappings.weighted_intra_chain_dist_weighted_rsa,

                           mappings.weighted_inter_chain_interactions_median, mappings.weighted_inter_chain_interactions_dist_weighted,
                           mappings.weighted_intra_chain_interactions_median, mappings.weighted_intra_chain_interactions_dist_weighted))))

    return values


def para_residue_init(rows):
    t0 = time.time()
    outs = []
    for row in rows:
        # Those residue inits include decoding of interaction profile and centrality score strings and thus takes some resources. For that a para function
        residue = sdsc.Residue(row[2], aa=row[3], lig_dist_str=row[4], chain_dist_str=row[5], RSA=row[6],
                               SSA=row[7], homo_dist_str=row[8], interaction_profile_str=row[9], centrality_score_str=row[10],
                               modres=row[11], b_factor=row[12], database_id=row[1], stored=True, phi=row[13], psi=row[14],
                               intra_ssbond=row[15], ssbond_length=row[16], intra_link=row[17], link_length=row[18],
                               cis_conformation=row[19], cis_follower=row[20], inter_chain_median_kd=row[21],
                               inter_chain_dist_weighted_kd=row[22], inter_chain_median_rsa=row[23],
                               inter_chain_dist_weighted_rsa=row[24], intra_chain_median_kd=row[25],
                               intra_chain_dist_weighted_kd=row[26], intra_chain_median_rsa=row[27],
                               intra_chain_dist_weighted_rsa=row[28])
        outs.append((row[0], row[2], residue))
    t1 = time.time()
    return outs, t1 - t0

# called by serializedPipeline


def getStoredResidues(proteins, config):
    t0 = time.time()

    stored_ids = proteins.getStoredStructureIds()

    if config.verbosity >= 2:
        t1 = time.time()
        print("Time for getstoredresidues 1: %s" % str(t1 - t0))

    if len(stored_ids) > 0:

        """
        rows = ['Structure', 'Residue_Id', 'Number', 'Amino_Acid', 'Sub_Lig_Dist', 'Sub_Chain_Distances',
                'Relative_Surface_Access', 'Relative_Surface_Access_Main_Chain', 'Relative_Surface_Access_Side_Chain',
                'Secondary_Structure_Assignment', 'Homomer_Distances',
                'Interaction_Profile', 'Centralities', 'Modres', 'B_Factor', 'PHI', 'PSI', 'Intra_SSBOND', 'SSBOND_Length',
                'Intra_Link', 'Link_Length', 'CIS_Conformation', 'CIS_Follower', 'Inter_Chain_Median_KD',
                'Inter_Chain_Dist_Weighted_KD', 'Inter_Chain_Median_RSA',
                'Inter_Chain_Dist_Weighted_RSA', 'Intra_Chain_Median_KD', 'Intra_Chain_Dist_Weighted_KD',
                'Intra_Chain_Median_RSA', 'Intra_Chain_Dist_Weighted_RSA',
                'Inter_Chain_Interactions_Median', 'Inter_Chain_Interactions_Dist_Weighted',
                'Intra_Chain_Interactions_Median', 'Intra_Chain_Interactions_Dist_Weighted',
                'Interacting_Chains', 'Interacting_Ligands'
                ]
        """
        rows = ['Structure', 'Residue_Id', 'Residue_Data']
        table = 'Residue'
        results = binningSelect(stored_ids.keys(), rows, table, config)

        if config.verbosity >= 2:
            t10 = time.time()
            print("Time for getstoredresidues 2.1: %s" % str(t10 - t1))

        for row in results:
            try:
                (res_id, one_letter, lig_dist_str, chain_dist_str, rsa, relative_main_chain_acc, relative_side_chain_acc,
                           ssa, homo_str, profile_str,
                           centrality_score_str, b_factor, modres, phi, psi, intra_ssbond, ssbond_length, intra_link, link_length,
                           cis_conformation, cis_follower, inter_chain_median_kd, inter_chain_dist_weighted_kd,
                           inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                           intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa,
                           inter_chain_interactions_median, inter_chain_interactions_dist_weighted,
                           intra_chain_interactions_median, intra_chain_interactions_dist_weighted,
                           interacting_chains_str, interacting_ligands_str) = unpack(row[2])
            except:
                config.errorlog.add_warning('Defective database entry in Residue table: %s %s' % (str(row[0]), str(row[1])))
                continue

            # Those residue inits include decoding of interaction profile and centrality score strings and thus takes some resources. For that a para function
            residue = sdsc.Residue(res_id, aa=one_letter, lig_dist_str=lig_dist_str, chain_dist_str=chain_dist_str, RSA=rsa,
                                   relative_main_chain_acc=relative_main_chain_acc, relative_side_chain_acc=relative_side_chain_acc,
                                   SSA=ssa, homo_dist_str=homo_str, interaction_profile_str=profile_str, centrality_score_str=centrality_score_str,
                                   modres=modres, b_factor=b_factor, database_id=row[1], stored=True, phi=phi, psi=psi,
                                   intra_ssbond=intra_ssbond, ssbond_length=ssbond_length, intra_link=intra_link, link_length=link_length,
                                   cis_conformation=cis_conformation, cis_follower=cis_follower, inter_chain_median_kd=inter_chain_median_kd,
                                   inter_chain_dist_weighted_kd=inter_chain_dist_weighted_kd, inter_chain_median_rsa=inter_chain_median_rsa,
                                   inter_chain_dist_weighted_rsa=inter_chain_dist_weighted_rsa, intra_chain_median_kd=intra_chain_median_kd,
                                   intra_chain_dist_weighted_kd=intra_chain_dist_weighted_kd, intra_chain_median_rsa=intra_chain_median_rsa,
                                   intra_chain_dist_weighted_rsa=intra_chain_dist_weighted_rsa,
                                   inter_chain_interactions_median=inter_chain_interactions_median, inter_chain_interactions_dist_weighted=inter_chain_interactions_dist_weighted,
                                   intra_chain_interactions_median=intra_chain_interactions_median, intra_chain_interactions_dist_weighted=intra_chain_interactions_dist_weighted,
                                   interacting_chains_str=interacting_chains_str, interacting_ligands_str=interacting_ligands_str)
            s_id = row[0]
            pdb_id, chain = stored_ids[s_id]
            proteins.add_residue(pdb_id, chain, res_id, residue)

    if config.verbosity >= 2:
        t2 = time.time()
        print("Time for getstoredresidues 2: %s" % str(t2 - t1))


def checkLigand(name, db, cursor):
    sql = "SELECT Ligand_Id FROM Ligand WHERE Name = '%s'" % name
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in checkLigand: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()
    if results == ():
        return 0
    elif len(results) == 0:
        return 0
    else:
        row = results[0]
        return row[0]


# called by Babel
def getLigandTemplates(name, db, cursor):
    ligand_id = checkLigand(name, db, cursor)
    sql = "SELECT Template FROM RS_Ligand_Template WHERE Ligand = '%s'" % str(ligand_id)
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getLigandTemplate: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()

    template_ids = []
    for row in results:
        template_ids.append(row[0])
    return template_ids


'''
#called by babel
def createLigandDB(outfile,session_id,db,cursor):
    sql = "SELECT Template,Mutation FROM RS_Annotation_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in createLigandDB: %s,\n%s" % (sql, f))
        db.rollback()

    template_ids = set()
    mutation_ids = set()
    for row in results:
        template_ids.add(row[0])
        mutation_ids.add(row[1])

    filtered_template_ids = set()
    if len(template_ids) > 0:
        sql = "SELECT Mutation,Template FROM RS_Mutation_Template WHERE Error IS NULL"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e, f, g] = sys.exc_info()
            raise NameError("Error in createLigandDB: %s,\n%s" % (sql, f))

        for row in results:
            if not row[1] in template_ids:
                continue
            if row[0] in mutation_ids:
                filtered_template_ids.add(row[1])

    ligand_ids = set()
    if len(filtered_template_ids) > 0:
        sql = "SELECT Ligand,Template FROM RS_Ligand_Template"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e, f, g] = sys.exc_info()
            raise NameError("Error in createLigandDB: %s,\n%s" % (sql, f))

        for row in results:
            if not row[1] in filtered_template_ids:
                continue
            ligand_ids.add(row[0])

    lines = []
    if len(ligand_ids) > 0:
        sql = "SELECT Name,Smiles,Ligand_Id FROM Ligand"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e, f, g] = sys.exc_info()
            raise NameError("Error in createLigandDB: %s,\n%s" % (sql, f))

        for row in results:
            if not row[2] in ligand_ids:
                continue
            lines.append("%s\t%s" % (row[1], row[0]))

    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()
'''


def getClassWT(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc, surface_t_id, chain_t_id, lig_t_id, config):
    if rel_sur_acc is not None:
        if rel_sur_acc > 1.0:
            rel_sur_acc = 0.01 * rel_sur_acc
    if lig_sub_dist is None:
        lig_sub_dist = "NONE"
    if chain_sub_dist is None:
        chain_sub_dist = "NONE"
    t_id = None
    try:
        mut_class = ""
        if (lig_sub_dist == "NONE" and chain_sub_dist == "NONE"):
            if rel_sur_acc is None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > config.surface_threshold:
                mut_class = "Surface isolated chain"
                t_id = surface_t_id
            else:
                mut_class = "Core isolated chain"
                t_id = surface_t_id
        elif lig_sub_dist == "NONE":
            if float(chain_sub_dist) > 5.0:
                if rel_sur_acc is None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > config.surface_threshold:
                    mut_class = "Surface"
                    t_id = surface_t_id
                else:
                    mut_class = "Core"
                    t_id = surface_t_id
            else:
                mut_class = "Contact %s" % chain_type
                t_id = chain_t_id
        elif chain_sub_dist == "NONE":
            if float(lig_sub_dist) > 5.0:
                if rel_sur_acc is None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > config.surface_threshold:
                    mut_class = "Surface"
                    t_id = surface_t_id
                else:
                    mut_class = "Core"
                    t_id = surface_t_id
            else:
                mut_class = "Contact Ligand"
                t_id = lig_t_id

        elif (float(lig_sub_dist) > 5.0 and float(chain_sub_dist) > 5.0):
            if rel_sur_acc is None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > config.surface_threshold:
                mut_class = "Surface"
                t_id = surface_t_id
            else:
                mut_class = "Core"
                t_id = surface_t_id
        elif (float(lig_sub_dist) < 5.0 and float(chain_sub_dist) > 5.0):
            mut_class = "Contact Ligand"
            t_id = lig_t_id
        elif (float(lig_sub_dist) > 5.0 and float(chain_sub_dist) < 5.0):
            mut_class = "Contact %s" % chain_type
            t_id = chain_t_id
        elif (float(lig_sub_dist) < 5.0 and float(chain_sub_dist) < 5.0):
            if (float(lig_sub_dist) <= float(chain_sub_dist)):
                mut_class = "Contact Ligand"
                t_id = lig_t_id
            elif (float(chain_sub_dist) < float(lig_sub_dist)):
                mut_class = "Contact %s" % chain_type
                t_id = chain_t_id

    except:
        mut_class = "Error"
        print("Error in class definition")
        print(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc)
        print(sys.exc_info())
    if mut_class == "":
        print(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc)

    return mut_class, t_id


def getClass(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc, config):
    if rel_sur_acc is not None:
        if rel_sur_acc > 1.0:
            rel_sur_acc = 0.01 * rel_sur_acc
    if lig_sub_dist is None:
        lig_sub_dist = "NONE"
    if chain_sub_dist is None:
        chain_sub_dist = "NONE"
    try:
        mut_class = ""
        if (lig_sub_dist == "NONE" and chain_sub_dist == "NONE"):
            if rel_sur_acc is None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > config.surface_threshold:
                mut_class = "Surface"
            else:
                mut_class = "Core"
        elif lig_sub_dist == "NONE":
            if float(chain_sub_dist) > 5.0:
                if rel_sur_acc is None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > config.surface_threshold:
                    mut_class = "Surface"
                else:
                    mut_class = "Core"
            else:
                mut_class = "Contact %s" % chain_type
        elif chain_sub_dist == "NONE":
            if float(lig_sub_dist) > 5.0:
                if rel_sur_acc is None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > config.surface_threshold:
                    mut_class = "Surface"
                else:
                    mut_class = "Core"
            else:
                mut_class = "Contact Ligand"

        elif (float(lig_sub_dist) > 5.0 and float(chain_sub_dist) > 5.0):
            if rel_sur_acc is None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > config.surface_threshold:
                mut_class = "Surface"
            else:
                mut_class = "Core"
        elif (float(lig_sub_dist) < 5.0 and float(chain_sub_dist) > 5.0):
            mut_class = "Contact Ligand"
        elif (float(lig_sub_dist) > 5.0 and float(chain_sub_dist) < 5.0):
            mut_class = "Contact %s" % chain_type
        elif (float(lig_sub_dist) < 5.0 and float(chain_sub_dist) < 5.0):
            if (float(lig_sub_dist) <= float(chain_sub_dist)):
                mut_class = "Contact Ligand"
            elif (float(chain_sub_dist) < float(lig_sub_dist)):
                mut_class = "Contact %s" % chain_type

    except:
        mut_class = "Error"
        print("Error in class definition")
        print(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc)
        print(sys.exc_info())
    if mut_class == "":
        print(lig_sub_dist, chain_sub_dist, chain_type, rel_sur_acc)

    return mut_class


def majority_vote(secs):
    class_dict = {}
    for (sec, qual) in secs:
        if sec not in class_dict:
            class_dict[sec] = qual
        else:
            class_dict[sec] += qual

    max_sec = None
    max_qual = 0.0
    for sec in class_dict:
        qual = class_dict[sec]
        if qual > max_qual:
            max_qual = qual
            max_sec = sec
    return max_sec


def getChemicalDistance(aac):
    try:
        if aac.count(',') < 1:
            aa1 = aac[0]
            aa2 = aac[-1]

            chemical_distance = sdsc.CHEM_DIST_MATRIX[aa1][aa2]
        else:
            aa1 = aac[0]
            aa2s = aac.split(",")
            aa2s[0] = aa2s[0][-1]
            chem_dists = []
            for aa2 in aa2s:
                chem_dists.append(sdsc.CHEM_DIST_MATRIX[aa1][aa2])
            chemical_distance = float(sum(chem_dists)) / float(len(chem_dists))
    except:
        return None
    return chemical_distance


def getBlosumValue(aac):
    if aac.count(',') < 1:
        try:
            try:
                blosum_value = sdsc.BLOSUM62[(aac[0], aac[-1])]
            except:
                blosum_value = sdsc.BLOSUM62[(aac[-1], aac[0])]
        except:
            blosum_value = 0.0
    else:
        aa1 = aac[0]
        aa2s = aac.split(",")
        aa2s[0] = aa2s[0][-1]
        bvs = []
        for aa2 in aa2s:
            try:
                try:
                    blosum_value = sdsc.BLOSUM62[(aa1, aa2)]
                except:
                    blosum_value = sdsc.BLOSUM62[(aa2, aa1)]
            except:
                blosum_value = 0.0
            bvs.append(blosum_value)
        blosum_value = float(sum(bvs)) / float(len(bvs))
    return blosum_value


'''
#called by babel
def getSLD(mutation_id,template_id,db,cursor):
    sql = "SELECT Sub_Lig_Dist FROM RS_Mutation_Template WHERE Template = '%s' AND Mutation = '%s'" % (str(template_id),str(mutation_id))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getSLD")
        # Rollback in case there is any NameError
        db.rollback()
    return results[0][0]
'''

'''
#called by babel
#structure of chosen_ones: {ligand_three_letter_code:tanimoto_score}
def getLigandAnnotation(chosen_ones,session_id,distance_threshold,db,cursor):
    anno_dict = {}

    ligand_map = {}
    ligand_ids = set()
    sql = "SELECT Ligand_Id,Name FROM Ligand"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in getLigandAnnotation: %s\n%s" % (sql, f))
    for row in results:
        ligand_name = row[1]
        if ligand_name not in chosen_ones:
            continue
        ligand_map[ligand_name] = row[0]
        ligand_ids.add(row[0])

    sql = "SELECT Template,Ligand FROM RS_Ligand_Template"
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getLigandAnnotation: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()

    ligand_template_map = {}
    template_ids = set()
    for row in results:
        if not row[1] in ligand_ids:
            continue
        template_ids.add(row[0])
        if not row[1] in ligand_template_map:
            ligand_template_map[row[1]] = set([row[0]])
        else:
            ligand_template_map[row[1]].add(row[0])

    sql = "SELECT Mutation,Template FROM RS_Annotation_Session WHERE Session = '%s'" % str(session_id)
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in getLigandAnnotation: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()

    template_mutation_map = {}
    for row in results:
        if not row[1] in template_ids:
            continue
        if not row[1] in template_mutation_map:
            template_mutation_map[row[1]] = {}

        template_mutation_map[row[1]][row[0]] = None

    min_t_id = min(template_ids)
    max_t_id = max(template_ids)

    t = 50000
    old_max = None
    if max_t_id - min_t_id > t:
        old_max = max_t_id
        max_t_id = min_t_id + t

    while True:
        sql = "SELECT Mutation,Template,Sub_Lig_Dist FROM RS_Mutation_Template WHERE Template BETWEEN %s AND %s" % (str(min_t_id), str(max_t_id))
        try:
            # Execute the SQL command
            cursor.execute(sql)
            results = cursor.fetchall()
            # Commit your changes in the database
            db.commit()
        except:
            raise NameError("Error in getLigandAnnotation: %s" % sql)
            # Rollback in case there is any NameError
            db.rollback()

        for row in results:
            if not row[1] in template_mutation_map:
                continue
            if not row[0] in template_mutation_map[row[1]]:
                continue
            template_mutation_map[row[1]][row[0]] = row[2]

        if old_max is None:
            break
        else:
            min_t_id = max_t_id + 1
            max_t_id = old_max
            if max_t_id - min_t_id > t:
                old_max = max_t_id
                max_t_id = min_t_id + t
            else:
                old_max = None

    for lig in chosen_ones:
        if lig not in ligand_map:
            continue
        lig_id = ligand_map[lig]
        filtered_annotation_ids = {}
        if lig_id not in ligand_template_map:
            continue
        for template_id in ligand_template_map[lig_id]:
            if template_id not in template_mutation_map:
                continue
            for mutation_id in template_mutation_map[template_id]:

                sub_lig_distances = template_mutation_map[template_id][mutation_id]
                if sub_lig_distances is not None:
                    sub_lig_distances = sub_lig_distances.split(",")
                    good_dists = []
                    for sub_lig_distance in sub_lig_distances:
                        if not sub_lig_distance == "-":
                            infos = sub_lig_distance.split(":")
                            dist = float(infos[1])
                            lig_infos = infos[0].split("_")
                            lig_name = lig_infos[0]
                            if lig_name == lig and dist <= float(distance_threshold):
                                good_dists.append(sub_lig_distance)
                    if len(good_dists) > 0:
                        good_dists = sorted(good_dists, key=lambda x: float(x.rsplit(":")[1]))
                        if template_id in filtered_annotation_ids:
                            filtered_annotation_ids[template_id].append((mutation_id, good_dists))
                        else:
                            filtered_annotation_ids[template_id] = [(mutation_id, good_dists)]
        if len(filtered_annotation_ids) > 0:
            for template_id in filtered_annotation_ids:
                filtered_annotation_ids[template_id] = sorted(filtered_annotation_ids[template_id], key=lambda x: float(x[1][0].split(":")[1]))
            anno_dict[lig] = (chosen_ones[lig], filtered_annotation_ids)
    return anno_dict
'''


def split(id_set, min_id, max_id):
    if len(id_set) == 0:
        return []
    if len(id_set) == 1:
        return [(id_set, min_id, max_id)]
    if len(id_set) == 2:
        return [(set([min_id]), min_id, min_id), (set([max_id]), max_id, max_id)]
    split_id = (min_id + max_id) // 2
    l = set()
    min_l = min_id
    max_l = min_id
    r = set()
    min_r = max_id
    max_r = max_id
    for i in id_set:
        if i < split_id:
            l.add(i)
            if i > max_l:
                max_l = i
        else:
            r.add(i)
            if i < min_r:
                min_r = i
    return [(l, min_l, max_l), (r, min_r, max_r)]


def split_fusion_binning(id_set, tablesize=None, binsize=500000, density=0.5, fusion=True):
    if tablesize == 0:
        return set(), []

    if tablesize is not None:
        if tablesize > 0:
            factor = 1
            density = float(len(id_set)) / (float(tablesize) * factor)
        else:
            density = 0.
    else:
        factor = 1

    bins = {}
    for i in id_set:
        bin_number = i // binsize
        if bin_number not in bins:
            bins[bin_number] = set()
        bins[bin_number].add(i)
    sorted_bins = []
    for bin_number in sorted(bins.keys()):
        id_set = bins[bin_number]
        sorted_bins.append((id_set, min(id_set), max(id_set)))

    dense_enough = False
    while not dense_enough:
        dense_enough = True
        split_bins = []
        for (id_set, min_id, max_id) in sorted_bins:
            min_amount = ((1 + max_id - min_id) * density)
            # split set if smaller than min_amount
            if len(id_set) < min_amount:
                for (iset, mi, ma) in split(id_set, min_id, max_id):
                    if len(iset) > 0:
                        split_bins.append((iset, mi, ma))
                    dense_enough = False
            else:
                split_bins.append((id_set, min_id, max_id))

        sorted_bins = split_bins

    if fusion:

        fusion_done = False
        while not fusion_done:
            fusion_done = True
            fused_bins = []
            last_bin_fused = False
            for pos, (id_set, min_id, max_id) in enumerate(sorted_bins):
                if pos == 0 or last_bin_fused:
                    last_bin_fused = False
                    if pos == (len(sorted_bins) - 1):
                        fused_bins.append((id_set, min_id, max_id))
                    continue
                pre_id_set, pre_min_id, pre_max_id = sorted_bins[pos - 1]
                if ((len(id_set) + len(pre_id_set))) > ((max_id - pre_min_id) * density * 2 * factor):

                    fused_bins.append(((pre_id_set | id_set), pre_min_id, max_id))
                    last_bin_fused = True
                    fusion_done = False
                else:
                    fused_bins.append((pre_id_set, pre_min_id, pre_max_id))
                    last_bin_fused = False
                    if pos == (len(sorted_bins) - 1):
                        fused_bins.append((id_set, min_id, max_id))

            sorted_bins = fused_bins

    singletons = []
    non_singleton_bins = []
    for (id_set, min_id, max_id) in sorted_bins:
        if min_id == max_id:
            singletons.append(min_id)
        else:
            non_singleton_bins.append((id_set, min_id, max_id))
    return singletons, non_singleton_bins


def median_focus_binning(id_set, density_thresh=0.5):
    # small sets are returned as a list of singletons
    if len(id_set) < 10:
        return list(id_set), []

    sorted_ids = sorted(id_set)

    # If the given set is dense enough, just return it as one interval
    density = len(sorted_ids) / (1 + sorted_ids[-1] - sorted_ids[0])
    if density > density_thresh:
        return [], [(id_set, sorted_ids[0], sorted_ids[-1])]

    l_quartile_pos = len(id_set) // 4
    r_quartile_pos = 3 * l_quartile_pos

    avg_quartile_dist = 1 + (2 * (1 + sorted_ids[r_quartile_pos] - sorted_ids[l_quartile_pos]) / len(sorted_ids))

    median_pos = len(id_set) // 2
    median = sorted_ids[median_pos]

    id_set = set([median])
    singletons = []

    l_minus_1_value = median
    r_plus_1_value = median  # needed for giving the boundaries of the set if only median is in the set
    for i in range(median_pos):
        l_value = sorted_ids[median_pos - i]
        l_minus_1_value = sorted_ids[median_pos - i - 1]
        if (l_value - l_minus_1_value) < avg_quartile_dist:
            id_set.add(l_minus_1_value)
        else:
            l_minus_1_value = l_value  # needed for giving the boundaries of the set
            break

    for j in range(median_pos - i):
        singletons.append(sorted_ids[j])

    for i in range(median_pos, len(sorted_ids) - 1):
        r_value = sorted_ids[i]
        r_plus_1_value = sorted_ids[i + 1]
        if (r_plus_1_value - r_value) < avg_quartile_dist:
            id_set.add(r_plus_1_value)
        else:
            r_plus_1_value = r_value
            break

    for j in range(i + 1, len(sorted_ids)):
        singletons.append(sorted_ids[j])

    if l_minus_1_value == r_plus_1_value:
        singletons.append(l_minus_1_value)
        return singletons, []

    return singletons, [(id_set, l_minus_1_value, r_plus_1_value)]


def getProtIdsFromSession(session_id, config, filter_mutant_proteins=False):
    table = 'RS_Protein_Session'
    cols = ['Protein', 'Input_Id']
    eq_cols = {'Session': session_id}
    results = select(config, cols, table, equals_rows=eq_cols)
    prot_ids = {}
    for row in results:
        if row[1] is None:
            continue
        prot_ids[row[0]] = row[1]
    return prot_ids

# called by output


def proteinsFromDb(session, config, with_residues=False, filter_mutant_proteins=False,
                   with_snvs=False, mutate_snvs=False, with_alignments=False,
                   with_complexes=False, keep_ray_alive=False):
    proteins = sdsc.Proteins({}, {}, {})  # create empty Proteins object

    prot_db_ids = getProtIdsFromSession(session, config, filter_mutant_proteins=filter_mutant_proteins)

    cols = ['Protein_Id', 'Primary_Protein_Id', 'Sequence']
    results = binningSelect(prot_db_ids.keys(), cols, 'Protein', config)
    id_prot_id_map = {}

    prot_id_list = set()
    prot_ids_mutants_excluded = set()

    for row in results:

        id_prot_id_map[row[0]] = row[1]
        prot_obj = sdsc.Protein(config.errorlog, primary_protein_id=row[1], database_id=row[0], input_id=prot_db_ids[row[0]], sequence=row[2])
        proteins[row[1]] = prot_obj
        prot_id_list.add(row[0])
        if prot_db_ids[row[0]] is not None:
            prot_ids_mutants_excluded.add(row[0])

    proteins.set_stored_ids(prot_id_list, prot_ids_mutants_excluded)

    cols = ['Protein', 'Position_Number', 'Position_Id', 'Recommended_Structure_Data']
    table = 'Position'

    results = binningSelect(prot_db_ids.keys(), cols, table, config)

    pos_db_map = {}

    for row in results:
        p_id = row[0]

        m_id = row[2]
        pos = row[1]

        recommended_structure, seq_id, cov, resolution = sdsc.process_recommend_structure_str(unpack(row[3])[0])
        pos_obj = sdsc.Position(pos=pos, checked=True, recommended_structure=recommended_structure, database_id=m_id)
        prot_id = id_prot_id_map[p_id]
        proteins[prot_id].add_positions([pos_obj])
        pos_db_map[m_id] = (prot_id, pos)

    if with_complexes:
        draw_complexes(config, proteins, draw_all=True)

    if with_snvs:
        cols = ['Position', 'New_AA']
        table = 'SNV'

        results = binningSelect(pos_db_map.keys(), cols, table, config)
        for row in results:
            snv = sdsc.SNV(row[1])
            (prot_id, pos) = pos_db_map[row[0]]
            proteins[prot_id].positions[pos].mut_aas[row[1]] = snv

    ray_init(config)

    if mutate_snvs:
        prot_ids = list(proteins.get_protein_ids()).copy()
        for prot_id in prot_ids:
            proteins[prot_id].mutate_snvs(proteins, config)

        serializedPipeline.autoTemplateSelection(config, proteins)

    if with_alignments:
        No_Errors = serializedPipeline.paraAlignment(config, proteins, skip_inserts=True, get_all_alignments=True)

    if not keep_ray_alive:
        ray.shutdown()

    if with_residues:
        getStoredResidues(proteins, config)

    return proteins


# method for comparing/developing the new classifiction
def diffSurfs(mutation_surface_dict, g_u_dict, mutation_dict, outfile, config):
    class_dict = {}
    for m in mutation_surface_dict:
        (u_ac, u_id) = g_u_dict[mutation_dict[m][1]]
        aac = mutation_dict[m][0]
        DSC = 0.0  # decision sum core
        DSS = 0.0  # decision sum surface
        n = 0.0
        min_surf = 2.0
        lines = ["Coverage\trASA"]
        for (surface, qual, cov) in mutation_surface_dict[m]:
            if surface < config.surface_threshold:
                DSC += qual * (cov**5)
            else:
                DSS += qual * (cov**5)
            lines.append("%s\t%s" % (str(cov), str(surface)))
            n += 1.0

            if surface < min_surf:
                min_surf = surface

        weighted_surface_value = DSS - 2 * DSC
        if weighted_surface_value > 0:
            weighted_sc = "Surface"
            if min_surf < config.surface_threshold:
                f = open("%s.%s_%s.tsv" % (outfile, u_ac, aac), 'w')
                f.write("\n".join(lines))
                f.close()
        else:
            weighted_sc = "Core"

        conf_sc = (1.0 - 1.0 / (n + 1.0)) * abs(DSS - 2 * DSC) / (DSS + 2 * DSC)


def createStructureDicts(proteins, config):
    if config.profiling:
        profile = cProfile.Profile()
        profile.enable()

    ligand_filter = config.ligand_filter

    u_acs = proteins.get_protein_ids()

    for u_ac in u_acs:
        package = []

        positions = proteins.get_position_ids(u_ac)

        annotation_list = proteins.get_protein_annotation_list(u_ac)

        iupred_map = proteins.get_disorder_scores(u_ac)

        for (pdb_id, chain) in annotation_list:

            sub_infos = proteins.get_sub_infos(u_ac, pdb_id, chain)

            resolution = proteins.get_resolution(pdb_id)
            chains = proteins.get_complex_chains(pdb_id)

            seq_id = proteins.get_sequence_id(u_ac, pdb_id, chain)
            if seq_id is None:
                continue
            cov = proteins.get_coverage(u_ac, pdb_id, chain)

            qual = templateFiltering.qualityScore(resolution, cov, seq_id)

            for pos in positions:
                if pos not in sub_infos:
                    continue
                aacbase = proteins.get_aac_base(u_ac, pos)
                sub_info = sub_infos[pos]
                res_nr = sub_info[0]

                if res_nr is None:
                    continue

                if not proteins.contains_residue(pdb_id, chain, res_nr):
                    continue
                res_aa = proteins.get_residue_aa(pdb_id, chain, res_nr)

                identical_aa = res_aa == aacbase[0]

                sld = proteins.get_residue_sld(pdb_id, chain, res_nr)
                scd = proteins.get_residue_scd(pdb_id, chain, res_nr)
                homomer_dists = proteins.get_residue_homomer_dists(pdb_id, chain, res_nr)
                centralities = proteins.get_residue_centralities(pdb_id, chain, res_nr)
                modres = proteins.get_residue_modres(pdb_id, chain, res_nr)
                b_factor = proteins.get_residue_b_factor(pdb_id, chain, res_nr)
                rsa = proteins.get_residue_rsa(pdb_id, chain, res_nr)
                ssa = proteins.get_residue_ssa(pdb_id, chain, res_nr)
                profile = proteins.get_residue_interaction_profile(pdb_id, chain, res_nr)

                phi = proteins.get_residue_phi(pdb_id, chain, res_nr)
                psi = proteins.get_residue_psi(pdb_id, chain, res_nr)

                intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower = proteins.get_residue_link_information(pdb_id, chain, res_nr)

                (inter_chain_median_kd, inter_chain_dist_weighted_kd,
                 inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                 intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa) = proteins.get_residue_milieu(pdb_id, chain, res_nr)

                min_hd, min_ld, min_md, min_id, min_cd, min_rd, min_dd, min_lig, min_metal, min_ion, iacs = proteins.structures[(pdb_id, chain)].residues[res_nr].get_shortest_distances(chains)

                minimal_distances = []
                if min_cd is not None:
                    minimal_distances.append(min_cd)
                if min_dd is not None:
                    minimal_distances.append(min_dd)
                if min_rd is not None:
                    minimal_distances.append(min_rd)
                if min_ld is not None:
                    minimal_distances.append(min_ld)
                if min_md is not None:
                    minimal_distances.append(min_md)
                if min_id is not None:
                    minimal_distances.append(min_id)

                if len(minimal_distances) == 0:
                    min_minimal_distances = 2.0
                else:
                    min_minimal_distances = min(minimal_distances)

                if min_minimal_distances < 1.2:
                    continue

                if rsa is None:
                    sc = None
                else:
                    if rsa > config.surface_threshold:
                        sc = "Surface"
                    else:
                        sc = "Core"

                raw_rin_class, raw_rin_simple_class = profile.getClass()

                if raw_rin_class == 'No interaction':
                    rin_class = sc
                else:
                    rin_class = raw_rin_class

                if raw_rin_simple_class == 'No interaction':
                    rin_simple_class = sc
                else:
                    rin_simple_class = raw_rin_simple_class

                #Class,conf = getWeightedClass(sc,1.0,min_cd,1.0,min_dd,1.0,min_rd,1.0,min_ld,1.0,min_md,1.0,min_id,1.0)
                #simpleClass =  simplifyClass(Class,sc)
                Class = rin_class
                simpleClass = rin_simple_class

                proteins.add_residue_classification(pdb_id, chain, res_nr, Class, simpleClass)

                mapping = (qual, seq_id, cov, rsa, ssa, min_ld, min_md, min_id, min_cd, min_rd, min_dd, min_hd, profile, centralities,
                           phi, psi, intra_ssbond, ssbond_length, intra_link, link_length, cis_conformation, cis_follower,
                           inter_chain_median_kd, inter_chain_dist_weighted_kd,
                           inter_chain_median_rsa, inter_chain_dist_weighted_rsa, intra_chain_median_kd,
                           intra_chain_dist_weighted_kd, intra_chain_median_rsa, intra_chain_dist_weighted_rsa, b_factor, modres,
                           Class, simpleClass, identical_aa)

                proteins.add_pos_res_mapping(u_ac, pos, pdb_id, chain, res_nr, mapping)

    if config.profiling:
        if cum_stats is not None:  # BUG: undefined variable
            cum_stats.add(profile)
            cum_stats.sort_stats('cumulative')
            cum_stats.print_stats()


def createInterDict(mutation_inter_dict, chain_type='sc'):

    inter_dict = {}

    for m_id in mutation_inter_dict:
        profiles = mutation_inter_dict[m_id]
        total_qual = 0.0
        ion_qual = 0.0
        metal_qual = 0.0
        lig_qual = 0.0
        chain_qual = 0.0
        average_profile = [0.0] * 14
        for profile_str, qual in profiles:
            if None in (profile_str, qual):
                continue
            profile = rin.Interaction_profile(profile_str=profile_str)
            total_qual += qual

            Ion_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'ion')
            Ion_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'ion')
            Metal_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'metal')
            Metal_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'metal')
            Ligand_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'ligand')
            Ligand_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'ligand')
            Chain_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'interchain')
            Chain_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'interchain')
            Short_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'neighbor')
            Short_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'neighbor')
            Medium_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'short')
            Medium_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'short')
            Long_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type, 'long')
            Long_Interaction_Score = profile.getChainSpecificCombiScore(chain_type, 'long')

            if Ligand_Interaction_Degree > 0:
                lig_qual += qual
            if Metal_Interaction_Degree > 0:
                metal_qual += qual
            if Ion_Interaction_Degree > 0:
                ion_qual += qual
            if Chain_Interaction_Degree > 0:
                chain_qual += qual

            average_profile[0] += Ion_Interaction_Degree * qual
            average_profile[1] += Ion_Interaction_Score * qual
            average_profile[2] += Metal_Interaction_Degree * qual
            average_profile[3] += Metal_Interaction_Score * qual
            average_profile[4] += Ligand_Interaction_Degree * qual
            average_profile[5] += Ligand_Interaction_Score * qual
            average_profile[6] += Chain_Interaction_Degree * qual
            average_profile[7] += Chain_Interaction_Score * qual
            average_profile[8] += Short_Interaction_Degree * qual
            average_profile[9] += Short_Interaction_Score * qual
            average_profile[10] += Medium_Interaction_Degree * qual
            average_profile[11] += Medium_Interaction_Score * qual
            average_profile[12] += Long_Interaction_Degree * qual
            average_profile[13] += Long_Interaction_Score * qual

        if total_qual > 0.0:
            average_profile[8] = average_profile[8] / total_qual
            average_profile[9] = average_profile[9] / total_qual
            average_profile[10] = average_profile[10] / total_qual
            average_profile[11] = average_profile[11] / total_qual
            average_profile[12] = average_profile[12] / total_qual
            average_profile[13] = average_profile[13] / total_qual

        if ion_qual > 0.0:
            average_profile[0] = average_profile[0] / ion_qual
            average_profile[1] = average_profile[1] / ion_qual
        if metal_qual > 0.0:
            average_profile[2] = average_profile[2] / metal_qual
            average_profile[3] = average_profile[3] / metal_qual
        if lig_qual > 0.0:
            average_profile[4] = average_profile[4] / lig_qual
            average_profile[5] = average_profile[5] / lig_qual
        if chain_qual > 0.0:
            average_profile[6] = average_profile[6] / chain_qual
            average_profile[7] = average_profile[7] / chain_qual
        inter_dict[m_id] = average_profile

    return inter_dict


def excludeFarClasses(c, sc):
    if c == "Surface" or c == "Core" or c == 'Disorder' or c is None:
        return c

    interactions = re.sub(r' Interaction$', '', re.sub(r'^[^:]*: ', '', c)).split(' and ')
    non_far_interactions = [x for x in interactions if ' far' not in x]

    if (len(non_far_interactions)) == 0:
        return sc

    clas = " and ".join(non_far_interactions)
    if len(non_far_interactions) == 4:
        clas = "Quadruple Interaction: " + clas
    elif len(non_far_interactions) == 3:
        clas = "Triple Interaction: " + clas
    elif len(non_far_interactions) == 2:
        clas = "Double Interaction: " + clas
    elif len(non_far_interactions) == 1:
        clas = clas + " Interaction"

    return clas


def writeInterFile(outfile, inter_dict, mutation_dict, protein_dict, new_aa_map, tag_map, class_dict, header=True):
    startline = "Uniprot\tAAC\tSpecie\tTag\tLigand_Interaction_Degree\tLigand_Interaction_Score\tChain_Interaction_Degree\tChain_Interaction_Score\tShort_Interaction_Degree\tShort_Interaction_Score\tMedium_Interaction_Degree\tMedium_Interaction_Score\tLong_Interaction_Degree\tLong_Interaction_Score\tClass\tComplex class"
    if header:
        lines = [startline]
    else:
        lines = []
    for m in inter_dict:
        aac, Protein_Id = mutation_dict[m][0:2]

        new_aa = new_aa_map[m]
        aac = "%s%s" % (aac.split(',')[0], new_aa)
        #(u_ac,u_id,species) = g_u_dict[mutation_dict[m][1]]
        (u_ac, gpan, u_id, error_code, error, species, input_id) = protein_dict[Protein_Id]

        (Class, conf, weighted_sc, conf_sc, best_res, max_seq_res, amount_of_structures,
         weighted_c, conf_c,
         weighted_d, conf_d,
         weighted_r, conf_r,
         weighted_l, conf_l,
         weighted_m, conf_m,
         weighted_i, conf_i,
         weighted_h, conf_h,
         max_seq_id,
         weighted_raw, weighted_cent, weighted_norm,
         weighted_lig_degree, weighted_lig_score,
         weighted_metal_degree, weighted_metal_score,
         weighted_ion_degree, weighted_ion_score,
         weighted_prot_degree, weighted_prot_score,
         weighted_rna_degree, weighted_rna_score,
         weighted_dna_degree, weighted_dna_score,
         weighted_modres, modres_prop, b_factor,
         intra_ssbond_prop, inter_ssbond_prop,
         intra_link_prop, inter_link_prop,
         cis_prop, cis_follower_prop,
         weighted_inter_chain_median_kd, weighted_inter_chain_dist_weighted_kd,
         weighted_inter_chain_median_rsa, weighted_inter_chain_dist_weighted_rsa,
         weighted_intra_chain_median_kd, weighted_intra_chain_dist_weighted_kd,
         weighted_intra_chain_median_rsa, weighted_intra_chain_dist_weighted_rsa) = class_dict[m]
        simple_class = simplifyClass(Class, weighted_sc)  # BUG: undefined variable

        if m in inter_dict:
            interstr = '\t'.join([str(x) for x in inter_dict[m]])
        else:
            interstr = '\t'.join(([''] * 14))

        line = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (u_ac, aac, species, tag_map[m], interstr, simple_class, Class)
        lines.append(line)

    f = open(outfile, 'a')
    f.write("\n".join(lines))
    f.close()


def writeClassFile(outfile, mutation_surface_dict, mutation_sec_dict, mutation_dict, protein_dict, class_dict, tag_map, header=True):
    startline = "Uniprot-Ac\tUniprot Id\tRefseq\tPDB-ID (Input)\tResidue-Id\tAmino Acid\tPosition\tSpecies\tTag\tWeighted Surface/Core\tClass\tSimple Class\tIndividual Interactions\tConfidence Value\tSecondary Structure\tRecommended Structure\tSequence-ID\tCoverage\tResolution\tMax Seq Id Structure\tMax Sequence-ID\tMax Seq Id Coverage\tMax Seq Id Resolution\tAmount of mapped structures"

    if header:
        lines = [startline]
    else:
        lines = []
    for m in class_dict:
        if m in mutation_sec_dict:
            mv_sec_ass = majority_vote(mutation_sec_dict[m])
        else:
            mv_sec_ass = None

        (Class, conf, weighted_sc, conf_sc, best_res, max_seq_res, amount_of_structures,
         weighted_c, conf_c,
         weighted_d, conf_d,
         weighted_r, conf_r,
         weighted_l, conf_l,
         weighted_m, conf_m,
         weighted_i, conf_i,
         weighted_h, conf_h,
         max_seq_id,
         weighted_raw, weighted_cent, weighted_norm,
         weighted_lig_degree, weighted_lig_score,
         weighted_metal_degree, weighted_metal_score,
         weighted_ion_degree, weighted_ion_score,
         weighted_prot_degree, weighted_prot_score,
         weighted_rna_degree, weighted_rna_score,
         weighted_dna_degree, weighted_dna_score,
         weighted_modres, modres_prop, b_factor,
         intra_ssbond_prop, inter_ssbond_prop,
         intra_link_prop, inter_link_prop,
         cis_prop, cis_follower_prop,
         weighted_inter_chain_median_kd, weighted_inter_chain_dist_weighted_kd,
         weighted_inter_chain_median_rsa, weighted_inter_chain_dist_weighted_rsa,
         weighted_intra_chain_median_kd, weighted_intra_chain_dist_weighted_kd,
         weighted_intra_chain_median_rsa, weighted_intra_chain_dist_weighted_rsa) = class_dict[m]
        simple_class = simplifyClass(Class, weighted_sc)  # BUG: undefined variable

        if best_res is not None:
            [r_id, qual, res_aa, res_nr, pdb_id, chain, resolution, cov, seq_id, rsa, min_lig, min_metal, min_ion, iacs] = best_res

            recommended_structure = '%s:%s %s:%s' % (pdb_id, chain, res_nr, res_aa)
        else:
            resolution = '-'
            cov = '-'
            seq_id = '-'
            recommended_structure = '-'

        if max_seq_res is not None:
            [max_seq_r_id, max_seq_qual, max_seq_res_aa, max_seq_res_nr, max_seq_pdb_id, max_seq_chain, max_seq_resolution, max_seq_cov, max_seq_seq_id, max_seq_rsa, max_min_lig, max_min_metal, max_min_ion, max_iacs] = max_seq_res

            max_seq_structure = '%s:%s %s:%s' % (max_seq_pdb_id, max_seq_chain, max_seq_res_nr, max_seq_res_aa)
        else:
            max_seq_resolution = '-'
            max_seq_cov = '-'
            max_seq_seq_id = '-'
            max_seq_structure = '-'

        aac, Protein_Id = mutation_dict[m][0:2]
        input_res_id = mutation_dict[m][4]
        if input_res_id is None:
            input_res_id = ''

        (u_ac, gpan, u_id, error_code, error, species, input_id) = protein_dict[Protein_Id]

        input_pdb_id = ''
        if len(u_ac) == 6 and u_ac[4] == ':':
            input_pdb_id = u_ac
            u_ac = ''
        interaction_str = '-'  # TODO

        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (u_ac, u_id, gpan, input_pdb_id, input_res_id, aac[0], aac[1:], species, tag_map[m], weighted_sc, Class, simple_class, interaction_str, str(conf), mv_sec_ass, recommended_structure, str(seq_id), str(cov), str(resolution), max_seq_structure, str(max_seq_seq_id), str(max_seq_cov), str(max_seq_resolution), str(amount_of_structures)))
    f = open(outfile, 'a')
    f.write("\n".join(lines))
    f.close()

# called by output, needs update, currently no gene_score


def goTermAnalysis(session_id, outfile, db, cursor):
    sql = "SELECT Protein FROM RS_Protein_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)
    datasetsize = len(results)
    go_dict = {}
    gene_dict = {}
    Protein_Id_list = set([])
    for row in results:
        gene_score = row[1]
        if gene_score is None:
            gene_score = 0.0
        gene_dict[row[0]] = gene_score
        Protein_Id_list.add(row[0])

    sql = "SELECT GO_Term,Protein FROM RS_Protein_GO_Term"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)
    go_term_id_list = set([])
    for row in results:
        if not row[1] in Protein_Id_list:
            continue
        gene_score = gene_dict[row[1]]
        if not row[0] in go_dict:
            go_term_id_list.add(row[0])
            go_dict[row[0]] = ["", "", gene_score, 1.0, gene_score]
        else:
            i = go_dict[row[0]][3] + 1
            avg_score = go_dict[row[0]][2]
            go_dict[row[0]][2] = ((i - 1) / i) * avg_score + (1 / i) * gene_score
            go_dict[row[0]][3] = i
            go_dict[row[0]][4] += gene_score

    sql = "SELECT GO_Term_Id,Name,Id FROM GO_Term"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)

    for row in results:
        if not row[0] in go_term_id_list:
            continue
        go_dict[row[0]][0] = row[1]
        go_dict[row[0]][1] = row[2]

    go_list = sorted(list(go_dict.values()), key=lambda x: x[4], reverse=True)
    lines = ["GO-Term\tGO-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
    for go in go_list:
        normalized_score = go[4] / datasetsize
        lines.append("%s\t%s\t%s\t%s\t%s\t%s" % (str(go[0]), str(go[1]), str(go[4]), str(go[2]), str(int(go[3])), str(normalized_score)))
    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()


# called by output, needs update, currently no gene_score
def pathwayAnalysis(session_id, outfile, db, cursor):
    sql = "SELECT Gene,Gene_Score FROM RS_Protein_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)
    datasetsize = len(results)
    path_dict = {}
    gene_dict = {}
    Protein_Id_list = set([])
    for row in results:
        gene_score = row[1]
        if gene_score is None:
            gene_score = 0.0
        gene_dict[row[0]] = gene_score
        Protein_Id_list.add(row[0])

    sql = "SELECT Pathway,Gene FROM RS_Protein_Pathway"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)
    pathway_id_list = set([])
    for row in results:
        if not row[1] in Protein_Id_list:
            continue
        gene_score = gene_dict[row[1]]
        if not row[0] in path_dict:
            pathway_id_list.add(row[0])
            path_dict[row[0]] = ["", "", gene_score, 1.0, gene_score]
        else:
            i = path_dict[row[0]][3] + 1
            avg_score = path_dict[row[0]][2]
            path_dict[row[0]][2] = ((i - 1) / i) * avg_score + (1 / i) * gene_score
            path_dict[row[0]][3] = i
            path_dict[row[0]][4] += gene_score

    if len(pathway_id_list) == 0:
        lines = ["Pathway\tReactome-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
        page = "\n".join(lines)
        f = open(outfile, "wb")
        f.write(page)
        f.close()
        return

    sql = "SELECT Pathway_Id,Name,Reactome_Id FROM Pathway"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)

    for row in results:
        if not row[0] in pathway_id_list:
            continue
        path_dict[row[0]][0] = row[1]
        path_dict[row[0]][1] = row[2]

    path_list = sorted(list(path_dict.values()), key=lambda x: x[4], reverse=True)
    lines = ["Pathway\tReactome-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
    for path in path_list:
        normalized_score = path[4] / datasetsize
        lines.append("%s\t%s\t%s\t%s\t%s\t%s" % (str(path[0]), str(path[1]), str(path[4]), str(path[2]), str(int(path[3])), str(normalized_score)))
    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()


# called by postAnnoAnno
def updateAnnoAnno(anno_anno_map, session_id, db, cursor):
    valuestrs = []
    for (m_id1, m_id2) in anno_anno_map:
        (min_d, atom, atom2, t_id1, t_id2, chain1, chain2) = anno_anno_map[(m_id1, m_id2)]
        valuestrs.append("('%s','%s','%s','%s','%s','%s','%s','%s','%s %s')" % (str(t_id1), str(t_id2), str(m_id1), str(m_id2), chain1, chain2, str(session_id), str(min_d), str(atom), str(atom2)))
    sql = "INSERT IGNORE INTO RS_Annotation_Annotation(Template_1,Template_2,Mutation_1,Mutation_2,Chain_1,Chain_2,Session,Distance,Atompair) VALUES %s" % ','.join(valuestrs)
    try:
        cursor.execute(sql)
        db.commit()
    except:
        [e, f, g] = sys.exc_info()
        raise NameError("Error in updateAnnoAnno: %s" % (f))


# called by calcAnnoRate
def calculateAnnotationRate(db, cursor, session_id):
    sql = "SELECT Position FROM RS_Position_Session WHERE Session = %s" % session_id
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in calculateAnnotationRate: %s" % sql)
    mut_id_list = set()
    for x in results:
        mut_id_list.add(str(x[0]))

    if len(mut_id_list) > 0:
        sql = "SELECT Position_Id,Gene FROM Position WHERE Position_Id in (%s)" % ",".join(mut_id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        gene_map = {}
        Protein_Id_list = set()
        for row in results:
            Protein_Id = str(row[1])
            mut_id = str(row[0])
            gene_map[mut_id] = Protein_Id
            Protein_Id_list.add(Protein_Id)

        sql = "SELECT Protein_Id,Error_Code FROM Protein WHERE Protein_Id in (%s)" % ",".join(Protein_Id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        error_map = {}
        for row in results:
            error_map[str(row[0])] = str(row[1])

        sql = "SELECT Position,Template,Error_Code FROM RS_Position_Template WHERE Position in (%s)" % ",".join(mut_id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        anno_mut_id_list = set()
        all_error_muts = set()
        fine_templates = set()

        for row in results:
            mut_id = str(row[0])
            if row[2] is None:
                anno_mut_id_list.add(mut_id)
                all_error_muts.discard(mut_id)
                fine_templates.add(row[1])
            elif mut_id not in anno_mut_id_list:
                all_error_muts.add(mut_id)

        sql = "SELECT Gene,Sequence_Identity,Template_Id FROM Template WHERE Gene in (%s)" % ",".join(Protein_Id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        number_of_templates = len(results)
        greater90 = set()
        greater95 = set()
        greater98 = set()
        greater99 = set()
        for row in results:
            if row[2] in fine_templates:
                g_id = str(row[0])
                s_id = float(row[1])
                if s_id > 99.0:
                    greater99.add(g_id)
                if s_id > 98.0:
                    greater98.add(g_id)
                if s_id > 95.0:
                    greater95.add(g_id)
                if s_id > 90.0:
                    greater90.add(g_id)

        muts_without_template = set()
        muts_with_unknown_gene_error = set()
        muts_with_anno = set()
        muts_with_temp_98 = 0
        for m in mut_id_list:
            if m in anno_mut_id_list:
                muts_with_anno.add(m)
                if gene_map[m] in greater98:
                    muts_with_temp_98 += 1
            if error_map[gene_map[m]] == '3':
                muts_without_template.add(m)
            elif error_map[gene_map[m]] == '4':
                muts_with_unknown_gene_error.add(m)

        # print "Total Mutations: ",len(mut_id_list)
        # print "Mutations with Annotation: ",len(muts_with_anno)
        # print "Mutations without template: ",len(muts_without_template)
        # print "Mutations with unknown gene error: ",len(muts_with_unknown_gene_error)
        # print "Mutations, which are mapped to gaps in all templates: ",len(all_error_muts)
        # print "Annotation-Rate: ",float(len(muts_with_anno))/float(len(mut_id_list))

        return (len(mut_id_list), len(muts_with_anno), len(muts_without_template), len(muts_with_unknown_gene_error), len(all_error_muts), float(len(muts_with_anno)) / float(len(mut_id_list)), number_of_templates, len(greater90), len(greater95), len(greater98), len(greater99), muts_with_temp_98)
    else:
        return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
