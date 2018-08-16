import urllib,urllib2
import time

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
def humanIdMapping(u_acs,ac_map,map_folder,obsolete_check=True):

    #print u_acs
    #print ac_map.keys()

    iso_map = {}
    for u_ac in u_acs:
        if u_ac.count('-') == 1:
            [u_ac,iso] = u_ac.split('-')
        else:
            iso = 'c'
        if not u_ac in iso_map:
            iso_map[u_ac] = [iso]
        else:
            iso_map[u_ac].append(iso)

    u_ac_id_map = {}

    u_ac_refseq_map = {}

    if obsolete_check:
        obs_set = set([])
        upd_set = {}
        obs_map = {}
        f = open('%s/human_obsolete_acs.tab' % map_folder,'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            [obs_ac,upd_ac] = line[:-1].split('\t')
            if obs_ac in iso_map:
                obs_set.add(obs_ac)
                upd_set[upd_ac] = iso_map[obs_ac]
                obs_map[upd_ac] = ac_map[obs_ac]
        for obs_ac in obs_set:
            del ac_map[obs_ac]
            del iso_map[obs_ac]
        iso_map.update(upd_set)
        ac_map.update(obs_map)
    
    f = open('%s/human_ac_mapping.tab' % map_folder,'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        words = line[:-1].split('\t')
        if words[0] in iso_map:
            ref_map = {}
            if words[2] != '':
                refs = words[2][:-1].split(';')
                for ref in refs:
                    [iso,refdouble] = ref.split(':')
                    if not iso in ref_map:
                        ref_map[iso] = [refdouble]
                    else:
                        ref_map[iso].append(refdouble)
            for iso in iso_map[words[0]]:
                if iso == 'c':
                    iso_ac = words[0]
                else:
                    iso_ac = '%s-%s' % (words[0],iso)
                u_ac_id_map[iso_ac] = words[1]
                if iso in ref_map:
                    ref_str = ','.join(ref_map[iso])
                    u_ac_refseq_map[iso_ac] = ref_str

    #print ac_map.keys()
    #print u_ac_id_map

    return u_ac_id_map,u_ac_refseq_map,ac_map

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
                if not words[0] in query_ids:
                    error = True
                    break
                uniprot_ids[words[0]] = words[1]
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

#called by serializePipeline
def getHumanSequences(u_acs,u_ac_refseq_map,info_map_path,human_proteome_path):
    f = open(human_proteome_path,'r')
    lines = f.readlines()
    f.close()

    canon_set = set([])
    control_set = set([])
    for u_ac in u_acs:
        canon_set.add(u_ac.split('-')[0])
        control_set.add(u_ac)
    canon_map = {}
    canon = False

    seq_go = False
    seq = ''
    gene_sequence_map = {}
    for line in lines:
        if line[0] == '>':
            if seq_go:
                if canon:
                    canon_map[u_ac] = seq
                else:
                    gene_sequence_map[u_ac] = seq
                seq_go = False
                seq = ''
            u_ac = line.split('|')[1]
            if u_ac.count('-') == 0:
                canon = True
                if u_ac in canon_set:
                    seq_go = True
            else:
                canon = False
                if u_ac in u_acs:
                    seq_go = True
                    control_set.remove(u_ac)
        elif seq_go:
            seq += line[:-1]
    if seq_go:
        if canon:
            canon_map[u_ac] = seq
        else:
            gene_sequence_map[u_ac] = seq
        gene_sequence_map[u_ac] = seq

    for u_ac in control_set:
        gene_sequence_map[u_ac] = canon_map[u_ac.split('-')[0]]

    #print u_acs
    #print gene_sequence_map.keys()

    f = open(info_map_path,'r')
    lines = f.readlines()
    f.close()

    iso_map = {}
    for u_ac in u_acs:
        if u_ac.count('-') == 1:
            [u_ac,iso] = u_ac.split('-')
        else:
            iso = 'c'
        if not u_ac in iso_map:
            iso_map[u_ac] = [iso]
        else:
            iso_map[u_ac].append(iso)

    gene_info_map = {}
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
                if iso_u_ac in u_ac_refseq_map:
                    refseqs = u_ac_refseq_map[iso_u_ac]
                else:
                    refseqs = ''
                gene_info_map[u_acs[iso_u_ac]] = (refseqs,go_terms,pathways,gene_sequence_map[iso_u_ac])

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
        if tries < 4:
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
