import urllib,urllib2
import time
import MySQLdb
import sys

def getUniprotId(query,querytype):
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
    'from':'%s' % (querytype),
    'to':'ID',
    'format':'tab',
    'query':'%s' % (query)
    } 
    #print params
    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = "" # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    try:
        response = urllib2.urlopen(request)
    except:
        return "-"
    page = response.read(200000)
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

#called by serializedPipeline
def IdMapping(ac_map,id_map,np_map,db,cursor,tag_map,species_map):
    genes = {}
    new_tag_map = {}
    new_species_map = {}
    for ac in ac_map:
        genes[ac] = ['',set([]),ac_map[ac]]
        if ac in tag_map:
            new_tag_map[ac] = tag_map[ac]
        if ac in species_map:
            new_species_map[ac] = species_map[ac]

    #Step one: map everything to uniprot-ac
    if len(id_map) > 0:
        if db != None:
            sql = "SELECT Uniprot_Ac,Uniprot_Id FROM AC_ID WHERE Uniprot_Id IN ('%s')" % "','".join(id_map.keys())
            try:
                cursor.execute(sql)
                results = cursor.fetchall()
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in mutationCheck: %s,%s" % (sql,f))

            for row in results:
                u_ac = row[0]
                u_id = row[1]
                if not u_ac in genes:
                    genes[u_ac] = [u_id,set([]),id_map[u_id]]
                    if u_id in tag_map:
                        new_tag_map[u_ac] = tag_map[u_id]

                    if u_id in species_map:
                        new_species_map[u_ac] = species_map[u_id]

                else:
                    genes[u_ac][0] = u_id
                    genes[u_ac][2] = genes[u_ac][2]|id_map[u_id]
                    new_tag_map[u_ac] = tag_map[u_id]

                    if u_id in tag_map:
                        if u_ac in new_tag_map:
                            new_tag_map[u_ac].update(tag_map[u_id])
                        else:
                            new_tag_map[u_ac] = tag_map[u_id]

                    if u_id in species_map:
                        new_species_map[u_ac] = species_map[u_id]

        else:
            id_ac_map = getUniprotIds(id_map.keys(),'ID',target_type="ACC")
            for u_id in id_ac_map:
                u_ac = id_ac_map[u_id]
                if not u_ac in genes:
                    genes[u_ac] = [u_id,set([]),id_map[u_id]]
                    if u_id in tag_map:
                        new_tag_map[u_ac] = tag_map[u_id]

                    if u_id in species_map:
                        new_species_map[u_ac] = species_map[u_id]

                else:
                    genes[u_ac][0] = u_id
                    genes[u_ac][2] = genes[u_ac][2]|id_map[u_id]
                    new_tag_map[u_ac] = tag_map[u_id]

                    if u_id in tag_map:
                        if u_ac in new_tag_map:
                            new_tag_map[u_ac].update(tag_map[u_id])
                        else:
                            new_tag_map[u_ac] = tag_map[u_id]

                    if u_id in species_map:
                        new_species_map[u_ac] = species_map[u_id]

    if len(np_map) > 0:
        if db != None:
            sql = "SELECT Uniprot_Ac,Refseq FROM AC_Refseq WHERE Refseq IN ('%s')" % "','".join(np_map.keys())
            try:
                cursor.execute(sql)
                results = cursor.fetchall()
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in mutationCheck: %s,%s" % (sql,f))

            for row in results:
                u_ac = row[0]
                ref = row[1]
                if not u_ac in genes:
                    genes[u_ac] = ['',set([ref]),np_map[ref]]
                    if ref in tag_map:
                        new_tag_map[u_ac] = tag_map[ref]

                    if ref in species_map:
                        new_species_map[u_ac] = species_map[ref]
                else:
                    genes[u_ac][1].add(ref)
                    genes[u_ac][2] = genes[u_ac][2]|np_map[ref]

                    if ref in tag_map:
                        if u_ac in new_tag_map:
                            new_tag_map[u_ac].update(tag_map[ref])
                        else:
                            new_tag_map[u_ac] = tag_map[ref]

                    if ref in species_map:
                        new_species_map[u_ac] = species_map[ref]

        else:
            np_ac_map = getUniprotIds(np_map.keys(),'P_REFSEQ_AC',target_type="ACC")
            for ref in np_ac_map:
                u_ac = np_ac_map[ref]
                if not u_ac in genes:
                    genes[u_ac] = ['',set([ref]),np_map[ref]]
                    if ref in tag_map:
                        new_tag_map[u_ac] = tag_map[ref]

                    if ref in species_map:
                        new_species_map[u_ac] = species_map[ref]
                else:
                    genes[u_ac][1].add(ref)
                    genes[u_ac][2] = genes[u_ac][2]|np_map[ref]

                    if ref in tag_map:
                        if u_ac in new_tag_map:
                            new_tag_map[u_ac].update(tag_map[ref])
                        else:
                            new_tag_map[u_ac] = tag_map[ref]

                    if ref in species_map:
                        new_species_map[u_ac] = species_map[ref]

    #Step two: get uniprot-id and refseqs from uniprot-ac

    ac_iso_map = {}
    id_search = []
    for u_ac in genes:
        if genes[u_ac][0] != '':
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
            sql = "SELECT Uniprot_Ac,Uniprot_Id FROM AC_ID WHERE Uniprot_Ac IN ('%s')" % "','".join(id_search)
            try:
                cursor.execute(sql)
                results = cursor.fetchall()
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in mutationCheck: %s,%s" % (sql,f))

            for row in results:
                u_ac = row[0]
                u_id = row[1]
                if u_ac in genes:
                    genes[u_ac][0] = u_id
                if u_ac in ac_iso_map:
                    for iso in ac_iso_map[u_ac]:
                        genes['%s-%s' % (u_ac,iso)][0] = u_id
        else:
            ac_id_map = getUniprotIds(id_search,'ACC',target_type="ID")
            for u_ac in ac_id_map:
                u_id = ac_id_map[u_ac]
                if u_ac in genes:
                    genes[u_ac][0] = u_id
                if u_ac in ac_iso_map:
                    for iso in ac_iso_map[u_ac]:
                        genes['%s-%s' % (u_ac,iso)][0] = u_id

    if db != None:
        sql = "SELECT Uniprot_Ac,Refseq FROM AC_Refseq WHERE Uniprot_Ac IN ('%s')" % "','".join(genes.keys())
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in mutationCheck: %s,%s" % (sql,f))

        for row in results:
            u_ac = row[0]
            ref = row[1]
            genes[u_ac][1].add(ref)
    else:
        ac_np_map = getUniprotIds(genes.keys(),'ACC',target_type="P_REFSEQ_AC")
        for u_ac in ac_np_map:
            ref = ac_np_map[u_ac]
            genes[u_ac][1].add(ref)

    return genes,new_tag_map,new_species_map


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
    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = "agress@mpi-inf.mpg.de" # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    try:
        response = urllib2.urlopen(request)
    except:
        print "ERROR: Uniprot did not answer" 
        return {}
    page = response.read(2000000)
    #print page
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
def getSequencesPlain(u_acs,db,cursor):
    gene_sequence_map = {}
    if db != None:
        if len(u_acs) == 0:
            return {}
        t0 = time.time()
        sql = "SELECT Uniprot_Ac,Sequence FROM Sequences WHERE Uniprot_Ac IN ('%s')" % "','".join(u_acs)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in mutationCheck: %s,%s" % (sql,f))

        t1 = time.time()
        print "getSequences Part 1: ",str(t1-t0)

        
        for row in results:
            u_ac = row[0]
            if not u_ac in u_acs:
                continue
            seq = row[1]
            gene_sequence_map[u_ac] = seq

        t2 = time.time()
        print "getSequences Part 2: ",str(t2-t1)

    else:
        t2 = time.time()

    missing_set = set()

    for u_ac in u_acs:
        if not u_ac in gene_sequence_map:
            missing_set.add(u_ac)
            #print u_ac

    t3 = time.time()
    print "getSequences Part 3: ",str(t3-t2)

    #print len(missing_set)
    #sys.exit()
    for u_ac in missing_set:
        seq,refseqs,go_terms,pathways = getSequence(u_ac)
        gene_sequence_map[u_ac] = seq

    t4 = time.time()
    print "getSequences Part 4: ",str(t4-t3)

    return gene_sequence_map

#called by serializePipeline
def getSequences(u_acs,info_map_path,db,cursor,pdb_dict):
    #print u_acs
    t0 = time.time()

    gene_sequence_map = getSequencesPlain(u_acs,db,cursor)

    
    f = open(info_map_path,'r')
    lines = f.readlines()
    f.close()

    gene_info_map = {}
    iso_map = {}
    for u_ac in u_acs:
        #if u_ac in missing_set:
        #    continue
        if u_ac.count('-') == 1:
            [base,iso] = u_ac.split('-')
        else:
            base = u_ac
            iso = 'c'
        if not base in iso_map:
            iso_map[base] = [iso]
        else:
            iso_map[base].append(iso)
        gene_info_map[u_acs[u_ac]] = ({},{},gene_sequence_map[u_ac])

    t4 = time.time()
    print "Time for extracting from Fasta Part 4: %s" % str(t4-t0)

    
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
                    go_id = ':'.join(go_double.split(':')[:2])
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
                gene_info_map[u_acs[iso_u_ac]] = (go_terms,pathways,gene_sequence_map[iso_u_ac])

    t5 = time.time()
    print "Time for extracting from Fasta Part 5: %s" % str(t5-t4)

    """
    #print missing_set
    #sys.exit()
    for u_ac in missing_set:
        seq,refseqs,go_terms,pathways = getSequence(u_ac)
        gene_sequence_map[u_ac] = seq
        gene_info_map[u_acs[u_ac]] = (refseqs,go_terms,pathways,seq)

    t6 = time.time()
    print "Time for extracting from Fasta Part 6: %s" % str(t6-t5)
    """

    #print gene_sequence_map

    return gene_sequence_map,gene_info_map

#called by serializedPipeline
def getSequence(uniprot_ac,tries=0):
    #new part just for the sequence
    time.sleep(2.0**tries-1.0)
    if len(uniprot_ac) < 2:
        return (0,"",{},{})
    url = 'https://www.uniprot.org/uniprot/%s.fasta' %uniprot_ac
    try:
        request = urllib2.Request(url)
        response = urllib2.urlopen(request)
        page = response.read(9000000)
    except:
        if tries < 3:
            return getSequence(uniprot_ac,tries=tries+1)

        else:
            #print uniprot_ac
            return (1,"",{},{})
    #print uniprot_ac,page
    lines = page.split("\n")

    wildtype_sequences = []

    for line in lines:
        if line == '':
            continue
        if line[0] == '>':
            continue
        wildtype_sequences.append(line)

    wildtype_sequence = ("".join(wildtype_sequences)).replace(" ","").replace("\n","")

    #old part, now just for refseqs,go and reactome
    time.sleep(2.0**tries-1.0)
    if len(uniprot_ac) < 2:
        return (0,"",{},{})
    url = 'https://www.uniprot.org/uniprot/%s.txt' %uniprot_ac
    try:
        request = urllib2.Request(url)
        response = urllib2.urlopen(request)
        page = response.read(9000000)
    except:
        if tries < 4:
            return getSequence(uniprot_ac,tries=tries+1)

        else:
            #print uniprot_ac
            return (1,"",{},{})
    lines = page.split("\n")
    #print(lines)
    refseqs = {}
    go_terms = {}
    pathways = {}
    for line in lines:
        if line == '':
            continue
        words = line.split()
        if len(words) == 0:
            print uniprot_ac
            return (1,"",{},{})
        if words[0] == "DR":
            if len(words) > 1:
                if words[1] == "GO;":
                    words_ = line.split(";")
                    go_id = words_[1].strip(";")
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

    return (wildtype_sequence,refseqs,go_terms,pathways)
