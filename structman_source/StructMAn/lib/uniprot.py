import urllib.request, urllib.parse, urllib.error,urllib.request,urllib.error,urllib.parse
import time
import pymysql as MySQLdb
import sys
import os
import sdsc
import database

def getUniprotId(query,querytype):
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
    'from':'%s' % (querytype),
    'to':'ID',
    'format':'tab',
    'query':'%s' % (query)
    }
    #print params
    data = urllib.parse.urlencode(params).encode('utf-8')
    request = urllib.request.Request(url, data)
    contact = "" # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    try:
        response = urllib.request.urlopen(request)
    except:
        return "-"
    page = response.read(200000).decode('utf-8')
    #print page
    try:
        lines = page.split("\n")
        line = lines[1].split()
        uniprot_id = line[1]
    #If unusable result, try without version number
    except:
        query = query.split(".")
        if len(query) > 1:
            query = query[0]
            return getUniprotId(query,querytype)
        else:
            return "-"

    return uniprot_id

def tag_update(tag_map,u_ac,new_entry):
    if not u_ac in tag_map:
        tag_map[u_ac] = new_entry
        return tag_map
    else:
        old_entry = tag_map[u_ac]

    for aac in new_entry:
        if not aac in old_entry:
            tag_map[u_ac][aac] = new_entry[aac]
        else:
            tags = set([])
            for tag in old_entry[aac].split(','):
                tags.add(tag)
            for tag in new_entry[aac].split(','):
                tags.add(tag)
            tag_map[u_ac][aac] = ','.join(tags)
    return tag_map

def updateMappingDatabase(u_acs,db,config):
    cursor = db.cursor()
    ac_id_values = []
    ac_ref_values = []
    seq_values = []
    for u_ac in u_acs:
        seq_out = getSequence(u_ac,config,return_id=True)
        if seq_out == None:
            continue

        seq,refseqs,go_terms,pathways,u_id = seq_out
        if u_id == None: #This can happen for uniprot entries, which got deleted from uniprot
            print("Warning: Uniprot entry:",u_ac," not found, most probably the entry got deleted from uniprot")
            continue
        ac_id_values.append("('%s','%s','%s','%s')" % (u_ac,u_ac[-2:],u_id,u_id[:2]))
        for refseq in refseqs.split(','):
            ac_ref_values.append("('%s','%s','%s','%s')" % (u_ac,u_ac.split('-')[0][-2:],refseq,refseq[:2]))
        seq_values.append("('%s','%s')" % (u_ac,seq))

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
                [e,f,g] = sys.exc_info()
                raise NameError("Error in updateMappingDatabase: %s\n%s" % (sql,f))


        if len(ac_ref_values) > 0:

            sql = "INSERT IGNORE INTO AC_Refseq(Uniprot_Ac,Ac_Hash,Refseq,Refseq_Hash) VALUES %s" % (','.join(ac_ref_values))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in updateMappingDatabase: %s\n%s" % (sql,f))

        if len(seq_values) > 0:

            sql = "INSERT IGNORE INTO Sequences(Uniprot_Ac,Sequence) VALUES %s " % (','.join(seq_values))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in updateMappingDatabase: %s\n%s" % (sql,f))
    except:
        #this happens for Devs without insert rights to the mapping DB, just ignore for the moment
        return

    return

#called by serializedPipeline
def IdMapping(config,ac_map,id_map,np_map,pdb_map):
    def indel_insert(indels,ac):
        for indel in indels:
            indel_protein_name = indel.create_protein_name(ac)
            if indel_protein_name in proteins:
                while indel_protein_name in proteins:
                    indel_protein_name = indel.create_other_protein_name(ac,indel_protein_name)
                    if indel_protein_name == None:
                        config.errorlog.add_error('All indel protein names are reserved, duplicate indels in input?')
                        break
            indel_mut_protein = sdsc.Protein(u_ac = indel_protein_name)
            proteins[indel_protein_name] = indel_mut_protein
            indel.set_proteins(ac,indel_protein_name)
            indel_map.append(indel)
        return

    try:
        db_adress = config.db_adress
        db_user_name = config.db_user_name
        db_password = config.db_password
        db = MySQLdb.connect(db_adress,db_user_name,db_password,config.mapping_db)
        cursor = db.cursor()
    except:
        db = None
        cursor = None

    proteins = {}
    indel_map = []

    for ac in ac_map:
        positions = ac_map[ac][0]
        indels = ac_map[ac][1]

        protein = sdsc.Protein(u_ac=ac,positions = positions)
        proteins[ac] = protein
        indel_insert(indels,ac)

    #Step one: map everything to uniprot-ac
    if len(id_map) > 0:
        if db != None:
            sql = "SELECT Uniprot_Ac,Uniprot_Id FROM AC_ID WHERE Uniprot_Id IN ('%s')" % "','".join(list(id_map.keys()))
            try:
                cursor.execute(sql)
                results = cursor.fetchall()
                db.commit()
            except:
                config.errorlog.add_error("Database error in IdMapping")
            stored_ids = set()
            for row in results:
                u_ac = row[0]
                u_id = row[1]
                stored_ids.add(u_id)
                if not u_ac in proteins:
                    protein = sdsc.Protein(u_ac=u_ac,u_id=u_id,positions = id_map[u_id][0])
                    proteins[u_ac] = protein
                else:
                    proteins[u_ac].u_id = u_id
                    proteins[u_ac].add_positions(id_map[u_id][0])
                indel_insert(id_map[u_id][1],u_ac)

            unstored_ids = []
            for u_id in id_map:
                if not u_id in stored_ids:
                    unstored_ids.append(u_id)
            update_acs = []
            if len(unstored_ids) > 0: #whenever ids are left in the dict, go to uniprot and download all unmapped entries (this happens for newer uniprot entries, which are not yet in the local mapping database)

                #This part is identical to the part, when no local database is used
                id_ac_map = getUniprotIds(unstored_ids,'ID',target_type="ACC")
                for u_id in id_ac_map:
                    u_ac = id_ac_map[u_id]
                    if not u_ac in proteins:
                        protein = sdsc.Protein(u_ac=u_ac,u_id=u_id,positions = id_map[u_id][0])
                        proteins[u_ac] = protein

                    else:
                        proteins[u_ac].u_id = u_id
                        proteins[u_ac].add_positions(id_map[u_id][0])
                    indel_insert(id_map[u_id][1],u_ac)
                    #This part is different
                    if not sdsc.is_mutant_ac(u_ac):
                        update_acs.append(u_ac)
                #updateMappingDatabase(update_acs,db,config)

        else:
            id_ac_map = getUniprotIds(list(id_map.keys()),'ID',target_type="ACC")
            for u_id in id_ac_map:
                u_ac = id_ac_map[u_id]
                if not u_ac in proteins:
                    protein = sdsc.Protein(u_ac=u_ac,u_id=u_id,positions = id_map[u_id][0])
                    proteins[u_ac] = protein
                else:
                    proteins[u_ac].u_id = u_id
                    proteins[u_ac].add_positions(id_map[u_id][0])
                indel_insert(id_map[u_id][1],u_ac)

    if len(np_map) > 0:
        if db != None:
            results = database.select(config,['Uniprot_Ac','Refseq'],'AC_Refseq',in_rows={'Refseq':np_map.keys()},from_mapping_db=True)

            ref_u_ac_map = {} #different u_acs may have the same refseq, try to choose the right one, prefering u_acs containing '-'
            gene_id_snap = set(proteins.keys()) #snapshot of u_acs not originating from refseq-mapping

            stored_refs = set()

            for row in results:
                u_ac = row[0]
                ref = row[1]
                stored_refs.add(ref)
                if not ref in ref_u_ac_map:
                    ref_u_ac_map[ref] = u_ac
                    if not u_ac in proteins:
                        protein = sdsc.Protein(u_ac=u_ac,ref_ids=set([ref]),positions = np_map[ref][0])
                        proteins[u_ac] = protein
                    else:
                        proteins[u_ac].add_ref_id(ref)
                        proteins[u_ac].add_positions(np_map[ref][0])
                    indel_insert(np_map[ref][1],u_ac)

                elif u_ac in gene_id_snap:
                    if ref_u_ac_map[ref].count('-') == 0 and u_ac.count('-') > 0:
                        ref_u_ac_map[ref] = u_ac
                    if not u_ac in proteins:
                        protein = sdsc.Protein(u_ac=u_ac,ref_ids=set([ref]),positions = np_map[ref][0])
                        proteins[u_ac] = protein
                    else:
                        proteins[u_ac].add_ref_id(ref)
                        proteins[u_ac].add_positions(np_map[ref][0])
                    indel_insert(np_map[ref][1],u_ac)

                elif ref_u_ac_map[ref].count('-') == 0:
                    if u_ac.count('-') > 0:
                        #if the current u_ac does not contain a '-' and the new found u_ac contains a '-': replace the corresponding ids
                        old_ac = ref_u_ac_map[ref]
                        ref_u_ac_map[ref] = u_ac
                        del proteins[old_ac]
                        if not u_ac in proteins:
                            protein = sdsc.Protein(u_ac=u_ac,ref_ids=set([ref]),positions = np_map[ref][0])
                            proteins[u_ac] = protein
                        else:
                            proteins[u_ac].add_ref_id(ref)
                            proteins[u_ac].add_positions(np_map[ref][0])
                        indel_insert(np_map[ref][1],u_ac)
            #similar to uniprot-id mapping, we have to go to uniprot to get search for unstored refseq entries and if we find them, we have to update the local mapping database
            unstored_refs = []
            for ref in np_map:
                if not ref in stored_refs:
                    unstored_refs.append(ref)
            update_acs = []
            if len(unstored_refs) > 0:
                np_ac_map = getUniprotIds(unstored_refs,'P_REFSEQ_AC',target_type="ACC")
                for ref in np_ac_map:
                    u_ac = np_ac_map[ref]
                    if not u_ac in proteins:
                        protein = sdsc.Protein(u_ac=u_ac,ref_ids=set([ref]),positions = np_map[ref][0])
                        proteins[u_ac] = protein
                    else:
                        proteins[u_ac].add_ref_id(ref)
                        proteins[u_ac].add_positions(np_map[ref][0])
                    indel_insert(np_map[ref][1],u_ac)
                    if not sdsc.is_mutant_ac(u_ac):
                        update_acs.append(u_ac)
                #updateMappingDatabase(update_acs,db,config)

        else:
            np_ac_map = getUniprotIds(list(np_map.keys()),'P_REFSEQ_AC',target_type="ACC")
            for ref in np_ac_map:
                u_ac = np_ac_map[ref]
                if not u_ac in proteins:
                    protein = sdsc.Protein(u_ac=u_ac,ref_ids=set([ref]),positions = np_map[ref][0])
                    proteins[u_ac] = protein
                else:
                    proteins[u_ac].add_ref_id(ref)
                    proteins[u_ac].add_positions(np_map[ref][0])
                indel_insert(np_map[ref][1],u_ac)

    #Step two: get uniprot-id and refseqs from uniprot-ac

    ac_iso_map = {}
    id_search = []
    for u_ac in proteins:
        if proteins[u_ac].u_id != None:
            continue
        split = u_ac.split('-')
        if len(split) == 2:
            base,iso = split
            if not base in ac_iso_map:
                ac_iso_map[base] = [iso]
            else:
                ac_iso_map[base].append(iso)

            id_search.append(base)
        else:
            id_search.append(u_ac)

    if len(id_search) > 0:
        if db != None:
            stored_u_acs = set()
            sql = "SELECT Uniprot_Ac,Uniprot_Id FROM AC_ID WHERE Uniprot_Ac IN ('%s')" % "','".join(id_search)
            try:
                cursor.execute(sql)
                results = cursor.fetchall()
                db.commit()
            except:
                config.errorlog.add_error("Database error in IdMapping")

            for row in results:
                u_ac = row[0]
                u_id = row[1]

                stored_u_acs.add(u_ac)
                if u_ac in proteins:
                    proteins[u_ac].u_id = u_id

                if u_ac in ac_iso_map:
                    for iso in ac_iso_map[u_ac]:
                        proteins['%s-%s' % (u_ac,iso)].u_id = u_id
                        stored_u_acs.add('%s-%s' % (u_ac,iso))

            unstored_u_acs = []
            for u_ac in id_search:
                if not u_ac in stored_u_acs:
                    if not sdsc.is_mutant_ac(u_ac):
                        unstored_u_acs.append(u_ac)

            if len(unstored_u_acs) > 0:
                updateMappingDatabase(unstored_u_acs,db,config)

                sql = "SELECT Uniprot_Ac,Uniprot_Id FROM AC_ID WHERE Uniprot_Ac IN ('%s')" % "','".join(unstored_u_acs)
                try:
                    cursor.execute(sql)
                    results = cursor.fetchall()
                    db.commit()
                except:
                    config.errorlog.add_error("Database error in IdMapping")

                for row in results:
                    u_ac = row[0]
                    u_id = row[1]

                    if u_ac in proteins:
                        proteins[u_ac].u_id = u_id
                    if u_ac in ac_iso_map:
                        for iso in ac_iso_map[u_ac]:
                            proteins['%s-%s' % (u_ac,iso)].u_id = u_id

        else:
            ac_id_map = getUniprotIds(id_search,'ACC',target_type="ID")
            for u_ac in ac_id_map:
                u_id = ac_id_map[u_ac]
                if u_ac in proteins:
                    proteins[u_ac].u_id = u_id
                if u_ac in ac_iso_map:
                    for iso in ac_iso_map[u_ac]:
                        proteins['%s-%s' % (u_ac,iso)].u_id = u_id

    if db != None:
        sql = "SELECT Uniprot_Ac,Refseq FROM AC_Refseq WHERE Uniprot_Ac IN ('%s')" % "','".join(list(proteins.keys()))
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            config.errorlog.add_error("Database error in IdMapping")

        for row in results:
            u_ac = row[0]
            ref = row[1]
            proteins[u_ac].add_ref_id(ref)
    else:
        ac_np_map = getUniprotIds(list(proteins.keys()),'ACC',target_type="P_REFSEQ_AC")
        for u_ac in ac_np_map:
            ref = ac_np_map[u_ac]
            proteins[u_ac].add_ref_id(ref)

    for pdb_tuple in pdb_map:
        protein = sdsc.Protein(pdb_id=pdb_tuple,positions = pdb_map[pdb_tuple])
        proteins[pdb_tuple] = protein

    if db != None:
        db.close()

    return proteins,indel_map


#called by serializedPipeline
def getUniprotIds(query_ids,querytype,target_type="ID"):
    #print query_ids
    if len(query_ids) == 0:
        return {}
    query = ' '.join(query_ids)
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
    'from':'%s' % (querytype),
    'to':'%s' % (target_type),
    'format':'tab',
    'query':'%s' % (query)
    }
    #print params
    data = urllib.parse.urlencode(params).encode('utf-8')
    request = urllib.request.Request(url, data)
    contact = "agress@mpi-inf.mpg.de" # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    try:
        response = urllib.request.urlopen(request)
    except:
        print("ERROR: Uniprot did not answer")
        return {}
    page = response.read(2000000).decode('utf-8')
    uniprot_ids = {}
    try:
        error = False
        lines = page.split("\n")
        for line in lines[1:]:
            words = line.split()
            if len(words) > 1:
                quer = words[0]
                target = words[1]
                if not quer in query_ids:
                    error = True
                    break
                if quer == target:
                    target = words[3].split(';')[0]
                uniprot_ids[quer] = target
        if error:

            if len(query_ids) == 1:
                return {query_ids.pop():None}
            #try a divide and conqer solution
            set_A = set()
            set_B = set()
            while len(query_ids) > 0:
                set_A.add(query_ids.pop())
                if len(query_ids) > 0:
                    set_B.add(query_ids.pop())
            map_A = getUniprotIds(set_A,querytype,target_type=target_type)
            map_B = getUniprotIds(set_B,querytype,target_type=target_type)
            return map_A.update(map_B)

    #If unusable result, try without version number
    except:
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
    return uniprot_ids

#called by serializedPipeline
def getSequencesPlain(u_acs,config,max_seq_len=None,debug=False,filtering_db=None):
    gene_sequence_map = {}
    filtered_set = set()

    missing_set = set()

    try:
        if config.mapping_db != None:
            db = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,config.mapping_db)
            cursor = db.cursor()
        else:
            db = None
    except:
        db = None
        [e,f,g] = sys.exc_info()
        config.errorlog.add_warning('Connection to mapping DB failed\n%s' % f)

    if filtering_db != None:
        filtering_db_path,dbs = filtering_db

    if db != None:
        if len(u_acs) == 0:
            return {}
        if debug:
            t0 = time.time()
        sql = "SELECT Uniprot_Ac,Sequence,Disorder_Scores,Disorder_Regions FROM Sequences WHERE Uniprot_Ac IN ('%s')" % "','".join(u_acs)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            config.errorlog.add_warning("Error selecting in getSequencesPlain: %s,%s" % (sql,f))

        db.close()

        if debug:
            t1 = time.time()
            print("getSequences Part 1: ",str(t1-t0))

        n = 0
        for row in results:
            u_ac = row[0]
            if not u_ac in u_acs:
                continue
            seq = row[1]
            disorder_scores = row[2]
            disorder_regions = row[3]
            if max_seq_len != None:
                if len(seq) > max_seq_len:
                    filtered_set.add(u_ac)
                    n += 1
                    continue
            if disorder_scores != None:
                disorder_score_list = [float(x) for x in disorder_scores[1:-1].split(',')]
                disorder_scores = {}
                for pos,score in enumerate(disorder_score_list):
                    disorder_scores[pos+1] = score
                if len(disorder_regions) == 2: #this happens for completely globular proteins
                    disorder_regions_datastruct = []
                else:
                    disorder_regions = disorder_regions[:-2].replace('[','').replace(' ','')
                    disorder_regions_datastruct = []
                    for region in disorder_regions.split('],'):
                        lower_bound,upper_bound,region_type = region.split(',')
                        lower_bound = int(lower_bound)
                        upper_bound = int(upper_bound)
                        disorder_regions_datastruct.append((lower_bound,upper_bound,region_type))
            else:
                disorder_scores = None
                disorder_regions_datastruct = None

            gene_sequence_map[u_ac] = seq,disorder_scores,disorder_regions_datastruct

        if n > 0 and debug:
            print('Filtered ',n,' Sequences due to max length: ',max_seq_len)

        if debug:
            t2 = time.time()
            print("getSequences Part 2: ",str(t2-t1))

    elif debug:
        t2 = time.time()
        missing_set = set(u_acs)

    in_db = set()
    for u_ac in u_acs:
        if filtering_db != None:
            inside_all = True
            for db_name in dbs:
                folder_key = u_ac.split('-')[0][-2:]
                filename = '%s/%s/%s_%s_gpw.fasta.gz' % (filtering_db_path,folder_key,u_ac,db_name)
                if not os.path.isfile(filename):
                     inside_all = False
            if inside_all:
                in_db.add(u_ac)
                continue

        if u_ac in filtered_set:
            continue

        if not u_ac in gene_sequence_map:
            missing_set.add(u_ac)
            #print u_ac

    if debug:
        t3 = time.time()
        print("getSequences Part 3: ",str(t3-t2))

    if debug:
        print('Size of missing set: ',len(missing_set))
        print(missing_set)

    for u_ac in missing_set:
        seq_out = getSequence(u_ac,config)
        if seq_out == None:
            config.errorlog.add_warning('getSequence output is None for %s' % u_ac)
            gene_sequence_map[u_ac] = 0,None,None
            continue
        seq,refseqs,go_terms,pathways = seq_out
        gene_sequence_map[u_ac] = seq,None,None

    if debug:
        t4 = time.time()
        print("getSequences Part 4: ",str(t4-t3))

    if filtering_db != None:
        return gene_sequence_map,in_db
    return gene_sequence_map

#called by serializePipeline
def getSequences(proteins,config):

    t0 = time.time()
    u_acs = set()
    iso_map = {}
    protein_map = proteins.get_protein_map()
    for u_ac in protein_map:
        if protein_map[u_ac].sequence == None:
            u_acs.add(u_ac)

        if u_ac.count('-') == 1:
            [base,iso] = u_ac.split('-')
        else:
            base = u_ac
            iso = 'c'
        if not base in iso_map:
            iso_map[base] = [iso]
        else:
            iso_map[base].append(iso)

    gene_sequence_map = getSequencesPlain(u_acs,config)

    for u_ac in u_acs:
        protein_map[u_ac].sequence = gene_sequence_map[u_ac][0]
        proteins.set_disorder_scores(u_ac,gene_sequence_map[u_ac][1])
        proteins.set_disorder_regions(u_ac,gene_sequence_map[u_ac][2])

    info_map_path = '%s/human_info_map.tab' % config.human_id_mapping_path

    f = open(info_map_path,'r')
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
                    [reac_id,pathway] = path_double.split(':',1)
                    pathways[reac_id] = pathway
            for iso in iso_map[u_ac]:
                if iso == 'c':
                    iso_u_ac = u_ac
                else:
                    iso_u_ac = '%s-%s' % (u_ac,iso)
                protein_map[iso_u_ac].go_terms = go_terms
                protein_map[iso_u_ac].pathways = pathways

    return

def getSequence(uniprot_ac,config,tries=0,return_id=False):
    if config.verbosity >= 2:
        print('uniprot.getSequence for ',uniprot_ac)

    if sdsc.is_mutant_ac(uniprot_ac):
        config.errorlog.add_error('Cannot call getSequence with a mutant protein: %s' % uniprot_ac)
        return None

    #new part just for the sequence
    time.sleep(2.0**tries-1.0)
    if len(uniprot_ac) < 2:
        return None

    if not uniprot_ac[0:3] == 'UPI':
        url = 'https://www.uniprot.org/uniprot/%s.fasta' %uniprot_ac
    else:
        url = 'https://www.uniprot.org/uniparc/%s.fasta' %uniprot_ac
    try:
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request)
        page = response.read(9000000).decode('utf-8')
    except:
        print('Error trying to reach: ',url)
        return None

    lines = page.split("\n")

    wildtype_sequences = []

    for line in lines:
        if line == '':
            continue
        if line[0] == '>':
            continue
        wildtype_sequences.append(line)

    wildtype_sequence = ("".join(wildtype_sequences)).replace(" ","").replace("\n","")

    if uniprot_ac[0:3] == 'UPI':
        if return_id:
            return (wildtype_sequence,{},{},{},None)
        return (wildtype_sequence,{},{},{})

    #old part, now just for refseqs,go and reactome
    time.sleep(2.0**tries-1.0)
    if len(uniprot_ac) < 2:
        return (0,"",{},{})

    url = 'https://www.uniprot.org/uniprot/%s.txt' %uniprot_ac
    try:
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request)
        page = response.read(9000000).decode('utf-8')
    except:
        if tries < 4:
            return getSequence(uniprot_ac,config,tries=tries+1)

        else:
            #print uniprot_ac
            return (1,"",{},{})
    lines = page.split("\n")
    #print(lines)
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
            return (1,"",{},{})
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
                    reac_id = split[1].replace(" ","")
                    pathway = split[2][1:-1]
                    pathways[reac_id] = pathway
                #TODO filter out isoforms
                if words[1] == "RefSeq;":
                    if words[-1].count('[') > 0:
                        u_ac_iso = words[-1][1:-1]
                        refs = words[2:-1]
                    else:
                        u_ac_iso = uniprot_ac.split('-')[0]
                        refs = words[2:]
                    refs = [x[:-1] for x in refs]
                    if not u_ac_iso in refseqs:
                        refseqs[u_ac_iso] = []
                    refseqs[u_ac_iso] += refs
    if uniprot_ac in refseqs:
        refseqs = ",".join(refseqs[uniprot_ac])
    else:
        refseqs = ''
    if return_id:
        return (wildtype_sequence,refseqs,go_terms,pathways,u_id)
    return (wildtype_sequence,refseqs,go_terms,pathways)
