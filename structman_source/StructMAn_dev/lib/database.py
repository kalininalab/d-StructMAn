#TODO replace all sql query, similiar to first big select in mindistout
import MySQLdb
import pdbParser as pdb
import sys
import os
import time
import templateFiltering
from operator import itemgetter
from multiprocessing import Process, Queue, Manager
from Bio.PDB import *
from Bio.SubsMat import MatrixInfo

chem_dist_matrix =  {'A': {'A': 0.0, 'C': 0.996, 'E': 0.544, 'D': 0.637, 'G': 0.663, 'F': 0.663, 'I': 0.626, 'H': 0.711, 'K': 0.768, 'M': 0.616, 'L': 0.435, 'N': 0.708, 'Q': 0.641, 'P': 0.949, 'S': 0.634, 'R': 0.977, 'T': 0.657, 'W': 0.985, 'V': 0.646, 'Y': 0.973}, 'C': {'A': 0.996, 'C': 0.0, 'E': 1.051, 'D': 0.804, 'G': 0.901, 'F': 0.859, 'I': 0.856, 'H': 0.595, 'K': 1.153, 'M': 0.651, 'L': 1.093, 'N': 0.687, 'Q': 0.753, 'P': 1.184, 'S': 0.744, 'R': 1.028, 'T': 0.737, 'W': 0.97, 'V': 0.82, 'Y': 0.835}, 'E': {'A': 0.544, 'C': 1.051, 'E': 0.0, 'D': 0.416, 'G': 0.923, 'F': 0.743, 'I': 0.892, 'H': 0.545, 'K': 0.647, 'M': 0.578, 'L': 0.724, 'N': 0.647, 'Q': 0.532, 'P': 0.994, 'S': 0.811, 'R': 0.897, 'T': 0.844, 'W': 0.882, 'V': 0.96, 'Y': 0.973}, 'D': {'A': 0.637, 'C': 0.804, 'E': 0.416, 'D': 0.0, 'G': 0.66, 'F': 0.68, 'I': 0.835, 'H': 0.361, 'K': 0.648, 'M': 0.58, 'L': 0.803, 'N': 0.291, 'Q': 0.385, 'P': 0.747, 'S': 0.555, 'R': 0.793, 'T': 0.64, 'W': 0.76, 'V': 0.886, 'Y': 0.744}, 'G': {'A': 0.663, 'C': 0.901, 'E': 0.923, 'D': 0.66, 'G': 0.0, 'F': 0.813, 'I': 0.814, 'H': 0.82, 'K': 0.974, 'M': 0.902, 'L': 0.827, 'N': 0.576, 'Q': 0.78, 'P': 0.629, 'S': 0.452, 'R': 1.081, 'T': 0.601, 'W': 1.017, 'V': 0.812, 'Y': 0.875}, 'F': {'A': 0.663, 'C': 0.859, 'E': 0.743, 'D': 0.68, 'G': 0.813, 'F': 0.0, 'I': 0.414, 'H': 0.53, 'K': 0.775, 'M': 0.442, 'L': 0.439, 'N': 0.656, 'Q': 0.607, 'P': 0.732, 'S': 0.723, 'R': 0.871, 'T': 0.666, 'W': 0.379, 'V': 0.625, 'Y': 0.509}, 'I': {'A': 0.626, 'C': 0.856, 'E': 0.892, 'D': 0.835, 'G': 0.814, 'F': 0.414, 'I': 0.0, 'H': 0.673, 'K': 0.741, 'M': 0.602, 'L': 0.382, 'N': 0.717, 'Q': 0.611, 'P': 0.9, 'S': 0.602, 'R': 0.754, 'T': 0.469, 'W': 0.733, 'V': 0.239, 'Y': 0.578}, 'H': {'A': 0.711, 'C': 0.595, 'E': 0.545, 'D': 0.361, 'G': 0.82, 'F': 0.53, 'I': 0.673, 'H': 0.0, 'K': 0.669, 'M': 0.346, 'L': 0.758, 'N': 0.365, 'Q': 0.299, 'P': 0.883, 'S': 0.598, 'R': 0.684, 'T': 0.586, 'W': 0.602, 'V': 0.736, 'Y': 0.579}, 'K': {'A': 0.768, 'C': 1.153, 'E': 0.647, 'D': 0.648, 'G': 0.974, 'F': 0.775, 'I': 0.741, 'H': 0.669, 'K': 0.0, 'M': 0.844, 'L': 0.702, 'N': 0.604, 'Q': 0.412, 'P': 0.883, 'S': 0.656, 'R': 0.383, 'T': 0.605, 'W': 0.879, 'V': 0.777, 'Y': 0.71}, 'M': {'A': 0.616, 'C': 0.651, 'E': 0.578, 'D': 0.58, 'G': 0.902, 'F': 0.442, 'I': 0.602, 'H': 0.346, 'K': 0.844, 'M': 0.0, 'L': 0.639, 'N': 0.639, 'Q': 0.534, 'P': 1.024, 'S': 0.762, 'R': 0.903, 'T': 0.725, 'W': 0.637, 'V': 0.698, 'Y': 0.745}, 'L': {'A': 0.435, 'C': 1.093, 'E': 0.724, 'D': 0.803, 'G': 0.827, 'F': 0.439, 'I': 0.382, 'H': 0.758, 'K': 0.702, 'M': 0.639, 'L': 0.0, 'N': 0.8, 'Q': 0.682, 'P': 0.867, 'S': 0.729, 'R': 0.894, 'T': 0.665, 'W': 0.778, 'V': 0.53, 'Y': 0.786}, 'N': {'A': 0.708, 'C': 0.687, 'E': 0.647, 'D': 0.291, 'G': 0.576, 'F': 0.656, 'I': 0.717, 'H': 0.365, 'K': 0.604, 'M': 0.639, 'L': 0.8, 'N': 0.0, 'Q': 0.304, 'P': 0.675, 'S': 0.339, 'R': 0.635, 'T': 0.418, 'W': 0.744, 'V': 0.735, 'Y': 0.555}, 'Q': {'A': 0.641, 'C': 0.753, 'E': 0.532, 'D': 0.385, 'G': 0.78, 'F': 0.607, 'I': 0.611, 'H': 0.299, 'K': 0.412, 'M': 0.534, 'L': 0.682, 'N': 0.304, 'Q': 0.0, 'P': 0.849, 'S': 0.446, 'R': 0.447, 'T': 0.413, 'W': 0.737, 'V': 0.628, 'Y': 0.57}, 'P': {'A': 0.949, 'C': 1.184, 'E': 0.994, 'D': 0.747, 'G': 0.629, 'F': 0.732, 'I': 0.9, 'H': 0.883, 'K': 0.883, 'M': 1.024, 'L': 0.867, 'N': 0.675, 'Q': 0.849, 'P': 0.0, 'S': 0.734, 'R': 1.034, 'T': 0.805, 'W': 0.734, 'V': 1.021, 'Y': 0.676}, 'S': {'A': 0.634, 'C': 0.744, 'E': 0.811, 'D': 0.555, 'G': 0.452, 'F': 0.723, 'I': 0.602, 'H': 0.598, 'K': 0.656, 'M': 0.762, 'L': 0.729, 'N': 0.339, 'Q': 0.446, 'P': 0.734, 'S': 0.0, 'R': 0.662, 'T': 0.189, 'W': 0.924, 'V': 0.539, 'Y': 0.639}, 'R': {'A': 0.977, 'C': 1.028, 'E': 0.897, 'D': 0.793, 'G': 1.081, 'F': 0.871, 'I': 0.754, 'H': 0.684, 'K': 0.383, 'M': 0.903, 'L': 0.894, 'N': 0.635, 'Q': 0.447, 'P': 1.034, 'S': 0.662, 'R': 0.0, 'T': 0.555, 'W': 0.939, 'V': 0.735, 'Y': 0.626}, 'T': {'A': 0.657, 'C': 0.737, 'E': 0.844, 'D': 0.64, 'G': 0.601, 'F': 0.666, 'I': 0.469, 'H': 0.586, 'K': 0.605, 'M': 0.725, 'L': 0.665, 'N': 0.418, 'Q': 0.413, 'P': 0.805, 'S': 0.189, 'R': 0.555, 'T': 0.0, 'W': 0.883, 'V': 0.389, 'Y': 0.56}, 'W': {'A': 0.985, 'C': 0.97, 'E': 0.882, 'D': 0.76, 'G': 1.017, 'F': 0.379, 'I': 0.733, 'H': 0.602, 'K': 0.879, 'M': 0.637, 'L': 0.778, 'N': 0.744, 'Q': 0.737, 'P': 0.734, 'S': 0.924, 'R': 0.939, 'T': 0.883, 'W': 0.0, 'V': 0.932, 'Y': 0.474}, 'V': {'A': 0.646, 'C': 0.82, 'E': 0.96, 'D': 0.886, 'G': 0.812, 'F': 0.625, 'I': 0.239, 'H': 0.736, 'K': 0.777, 'M': 0.698, 'L': 0.53, 'N': 0.735, 'Q': 0.628, 'P': 1.021, 'S': 0.539, 'R': 0.735, 'T': 0.389, 'W': 0.932, 'V': 0.0, 'Y': 0.695}, 'Y': {'A': 0.973, 'C': 0.835, 'E': 0.973, 'D': 0.744, 'G': 0.875, 'F': 0.509, 'I': 0.578, 'H': 0.579, 'K': 0.71, 'M': 0.745, 'L': 0.786, 'N': 0.555, 'Q': 0.57, 'P': 0.676, 'S': 0.639, 'R': 0.626, 'T': 0.56, 'W': 0.474, 'V': 0.695, 'Y': 0.0}}

def getGeneScoreDict(gene_id_list,session_id,db,cursor):
    t0 = time.time()
    max_g = max(gene_id_list)
    min_g = min(gene_id_list)
    print max_g,min_g
    sql = "SELECT Gene,Gene_Score,Session FROM RS_Gene_Session WHERE Session = '%s' AND Gene <= %s AND Gene >= %s" % (str(session_id),str(max_g),str(min_g))
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("NameError in getGeneScoreDict: %s" % sql)
    score_dict = {}

    for row in results:
        if row[2] != session_id:
            continue
        g_id = row[0]
        if g_id in gene_id_list:
            score_dict[g_id] = row[1]
    t1 = time.time()
    sql = "SELECT Gene_Id,Uniprot_Ac,Genbank_Protein_Accession_Number,Uniprot_Id,Error_Code,Error,Species FROM Gene WHERE Gene_Id <= %s AND Gene_Id >= %s" % (str(max_g),str(min_g))
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("NameError in getGeneScoreDict: %s" % sql)
    gene_score_dict = {}
    for row in results:
        g_id = row[0]
        if g_id in gene_id_list and g_id in score_dict:
            gene_score_dict[g_id] = (row[1],row[2],row[3],row[4],row[5],row[6],score_dict[g_id])
    t2 = time.time()

    print "Time for getGeneScoreDict1: ",t1-t0
    print "Time for getGeneScoreDict2: ",t2-t1

    return gene_score_dict


#called from serializedPipeline
def geneScan(genes_aac_list,db,cursor):
    #Scan the database for stored Genes
    if len(genes_aac_list) == 0:
        return {},{}
    sql = "SELECT Uniprot_Ac FROM Gene"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in geneCheck: %s,\n%s" % (sql,f))
    
    t1 = time.time()
    stored_genes = {}
    new_genes = {}

    for row in results:
        u_id = row[0]
        if u_id in genes_aac_list:
            stored_genes[u_id] = genes_aac_list[u_id]

    for u_id in genes_aac_list:
        if not u_id in stored_genes:
            new_genes[u_id] = genes_aac_list[u_id]

    return stored_genes,new_genes


#called from serializedPipeline
def geneCheck(genes_aac_list,ac_id_map,species_map,session_id,db,cursor,seq_thresh,res_thresh,cov_thresh,GS_db):
    t0 = time.time()
    #Scan the database for stored Genes
    if len(genes_aac_list) == 0:
        return {},{}
    sql = "SELECT Gene_Id,Uniprot_Ac FROM Gene"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in geneCheck: %s,\n%s" % (sql,f))
    
    t1 = time.time()
    stored_genes = {}
    
    gene_id_list = set([])

    max_g_id = 0
    for row in results:
        gene_id = row[0]
        u_id = row[1]
        if not u_id in genes_aac_list:
            continue
        stored_genes[u_id] = [gene_id,None]
        gene_id_list.add(gene_id)
        if gene_id > max_g_id:
            max_g_id = gene_id

    t2 = time.time()
    gene_dict = {}
    #Get the corresponding Session-Ids
    if len(gene_id_list) > 0:
        sql = "SELECT Session,Gene FROM RS_Gene_Session"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in geneCheck: %s,\n%s" % (sql,f))

        session_id_list = set([])        
        
        for row in results:
            session_id = row[0]
            gene_id = row[1]
            if not gene_id in gene_id_list:
                continue

            session_id_list.add(session_id)

            if gene_id not in gene_dict:
                gene_dict[gene_id] = set()
                gene_dict[gene_id].add(session_id)
            else:
                gene_dict[gene_id].add(session_id)

        if len(session_id_list) == 0:
            print gene_id_list
            raise NameError("Error in geneCheck: Found Gene without Session")

        #Retrieve informations about the sessions
        sql = "SELECT Session_Id,Sequence_Identity_Threshold,Coverage_Threshold,Resolution_Threshold FROM Session"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in geneCheck: %s,\n%s" % (sql,f))

    t3 = time.time()
    #For the already stored genes, we have to figure out, if there was a more restrictive session
    restrictive_sessions = {}

    for gene in stored_genes:
        for row in results:
            session_id = row[0]
            if not session_id in session_id_list:
                continue
            if session_id in gene_dict[stored_genes[gene][0]]:
                if row[1] == None:
                    seq_id = 0.0
                else:
                    seq_id = float(row[1])
                if row[2] == None:
                    cov = 0.0
                else:
                    cov = float(row[2])
                if row[3] == None:
                    res = 100000.0
                else:
                    res = float(row[3])

                if float(seq_thresh) >= seq_id and float(res_thresh) <= res and float(cov_thresh) >= cov:
                    if stored_genes[gene][1] == None:
                        stored_genes[gene][1] = True                       
        if stored_genes[gene][1] == None:
            stored_genes[gene][1] = False
    t4 = time.time()
    #Insert the new genes into the database
    new_genes = set()
    for gene in genes_aac_list:
        if not gene in stored_genes:
            new_genes.add(gene)
    t5 = time.time()
    if len(new_genes) > 0:
        value_strs = []
        for gene in new_genes:
            if gene in species_map:
                value_strs.append("('%s','%s','%s','%s')" % (gene.replace("'","\\'"),ac_id_map[gene],str(session_id),species_map[gene]))
            else:
                value_strs.append("('%s','%s','%s',NULL)" % (gene.replace("'","\\'"),ac_id_map[gene],str(session_id)))

        size = len(value_strs)
        print size
        if size > 0:
            bound = 50000
            if size > bound:
                m = size/bound
                rest = size%bound
                if rest == 0:
                    part_size = size/m
                else:
                    m += 1
                    part_size = size/m
                parts = []
                for i in range(0,m-1):
                    parts.append(value_strs[(i*part_size):((i+1)*part_size)])
                parts.append(value_strs[(m-1)*part_size:])
            else:
                parts = [value_strs]
            for part in parts:

                sql = ("""INSERT INTO Gene(Uniprot_Ac,
                        Uniprot_Id,Original_Session,Species)
                        VALUES %s""") % ','.join(part)
                try:
                    cursor.execute(sql)
                    db.commit()
                except:
                    [e,f,g] = sys.exc_info()
                    raise NameError("Couldn't insert geneCheck: %s,%s" % (sql,f))


            
        #Retrieve the Gene-Ids from the new added genes
        sql = "SELECT Gene_Id,Uniprot_Ac FROM Gene WHERE Gene_Id > %s" % str(max_g_id)# WHERE Uniprot_Ac IN ('%s')" % "','".join(new_genes)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in geneCheck: %s,\n%s" % (sql,f))
    else:
        results = ()
    t6 = time.time()
    new_gene_map = {}
    new_gene_ids = set()
    for row in results:
        if row[1] in new_genes:
            new_gene_map[row[1]] = row[0]
            new_gene_ids.add(row[0])
    t7 = time.time()
    #Insert the Gene-Session-Connections into the database
    value_strs = []
    for gene in new_gene_map:
        value_strs.append("('%s','%s')" % (str(new_gene_map[gene]),str(session_id)))

    for gene in stored_genes:
        value_strs.append("('%s','%s')" % (str(stored_genes[gene][0]),str(session_id)))

    process = Process(target=backgroundInsertGS,args=(value_strs,GS_db))
    process.start()

    #Return the mappings uniprot-id:gene_id(and 'more restrictive?')
    t8 = time.time()

    print "Time for geneCheck Part 1: %s" % str(t1-t0)
    print "Time for geneCheck Part 2: %s" % str(t2-t1)
    print "Time for geneCheck Part 3: %s" % str(t3-t2)
    print "Time for geneCheck Part 4: %s" % str(t4-t3)
    print "Time for geneCheck Part 5: %s" % str(t5-t4)
    print "Time for geneCheck Part 6: %s" % str(t6-t5)
    print "Time for geneCheck Part 7: %s" % str(t7-t6)
    print "Time for geneCheck Part 8: %s" % str(t8-t7)
    return stored_genes,new_gene_map,gene_id_list,new_gene_ids,process

def backgroundInsertGS(value_strs,db):
    cursor = db.cursor()
    sql = ("""INSERT IGNORE INTO RS_Gene_Session(
                Gene, Session)
                VALUES %s""") % ','.join(value_strs)
    try:
        cursor.execute(sql)
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Couldn't insert geneCheck: %s,%s" % (sql,f))

def getLastAutoInc(db,cursor):
    sql = "SELECT LAST_INSERT_ID()"
         
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        for row in results:
            auto_id = row[0]
        db.commit()
    except:
        print "NameError in getLastAutoInc"
        db.rollback()
    return auto_id

def getUniprotFromId(gene_id,db,cursor):
    sql = "SELECT Uniprot_Id FROM Gene WHERE Gene_Id = '%s'" % (str(gene_id))
         
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

#called by babel
def getUniprotAcFromId(gene_id,db,cursor):
    sql = "SELECT Uniprot_Ac FROM Gene WHERE Gene_Id = '%s'" % (str(gene_id))
         
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


def getLigandDict(template_id_list,db,cursor):
    t0 = time.time()
    max_t = max(template_id_list)
    min_t = min(template_id_list)
    sql = "SELECT Ligand,Template,Chain,Residue FROM RS_Ligand_Template WHERE Template <= %s AND Template >= %s" % (str(max_t),str(min_t))
         
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getLigands: %s" % sql)
        db.rollback()

    if results == ():
        return []

    template_dict = {}
    ligand_id_list = set()
    for row in results:
        t_id = row[1]
        if t_id in template_id_list:
            if not t_id in template_dict:
                template_dict[t_id] = {row[0]:["Ligand","",row[3],row[2]]}
            elif not row[0] in template_dict[t_id]:
                template_dict[t_id][row[0]] = ["Ligand","",row[3],row[2]]
            ligand_id_list.add(row[0])

    t1 = time.time()


    max_l = max(ligand_id_list)
    min_l = min(ligand_id_list)
    print max_l,min_l
    if len(ligand_id_list) > 0:
        sql = "SELECT Ligand_Id,Name FROM Ligand WHERE Ligand_Id <= %s AND Ligand_Id >= %s" % (str(max_l),str(min_l))
             
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in getLigands: %s" % sql)
            db.rollback()
        
        lig_map = {}
        for row in results:
            lig_id = row[0]
            if lig_id in ligand_id_list:
                lig_map[lig_id] = row[1]


        for template_id in template_dict:
            for lig_id in template_dict[template_id]:
                template_dict[template_id][lig_id][1] = lig_map[lig_id]
    t2 = time.time()

    print "Time for getLigandDict1: ",t1-t0
    print "Time for getLigandDict2: ",t2-t1

    return template_dict

#called by babel
def getAAChange(mutation_id,db,cursor):
    sql = "SELECT Amino_Acid_Change FROM Mutation WHERE Mutation_Id = '%s'" % str(mutation_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getAAChange: %s" % sql)
        db.rollback()
    #print(results)
    if results == ():
        return ""
    for row in results:
       aac = row[0]
       return aac
    return ""

#called by serializedPipeline
def mutationCheck(gene_aaclist_map,stored_genes,stored_gene_ids,new_genes,new_gene_ids,database_session,tag_map,db,cursor,MS_db):
    #structure of stored_genes: {Uniprot-Id:(gene_id,more_restrictive)}
    #structure of new_genes: {Uniprot-Id:gene_id}
    gene_mut_map_new = {}
    gene_mut_map_pseudo = {}
    gene_mut_map_twins = {}
    t0 = time.time()

    #print stored_genes
    #print stored_gene_ids

    #stored_gene_ids = set([x[0] for x in stored_genes.values()])
    stored_mutation_map = {}
    update_map = {}
    max_g_id = None
    min_g_id = None
    binsize = 10000
    bins = set([])
    if len(stored_gene_ids) > 0:

        for g_id in stored_gene_ids:
            bin_number = g_id//binsize
            bins.add(bin_number)

            if max_g_id == None or int(g_id) > max_g_id:
                max_g_id = int(g_id)
            if min_g_id == None or int(g_id) < min_g_id:
                min_g_id = int(g_id)

        min_max_tuples = []
        if max_g_id - min_g_id > binsize:
            for bin_number in bins:
                min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize))
        else:
            min_max_tuples = [(min_g_id,max_g_id)]


        total_results = []
        for (run_min_g_id,run_max_g_id) in min_max_tuples:

            sql = "SELECT Amino_Acid_Change,Gene,Mutation_Id FROM Mutation WHERE Gene BETWEEN %s and %s" % (str(run_min_g_id),str(run_max_g_id))

            try:
                cursor.execute(sql)
                results = cursor.fetchall()
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in mutationCheck: %s,%s" % (sql,f))
            total_results += results


        for row in total_results:
            gene_id = row[1]
            if not gene_id in stored_gene_ids:
                continue
            aacs = row[0].split(",")
            aacbase = aacs[0][:-1]
            mutation_id = row[2]

            if not gene_id in stored_mutation_map:
                stored_mutation_map[gene_id] = set()
                stored_mutation_map[gene_id].add(aacbase)
            else:
                stored_mutation_map[gene_id].add(aacbase)

            if not gene_id in update_map:
                update_map[gene_id] = {aacbase:(mutation_id,row[0])}
            else:
                update_map[gene_id][aacbase] = (mutation_id,row[0])

    t1 = time.time()

    for gene in gene_aaclist_map:
        twin_set = set()
        aaclist = gene_aaclist_map[gene]
        if gene in stored_genes:
            for aac in aaclist:
                aacbase = aac[:-1]
                if aacbase in twin_set:
                    if gene not in gene_mut_map_twins:
                        gene_mut_map_twins[gene] = set()
                        gene_mut_map_twins[gene].add(aac)
                    else:
                        gene_mut_map_twins[gene].add(aac)
                else:
                    twin_set.add(aacbase)
                    if aacbase in stored_mutation_map[stored_genes[gene][0]]:
                        if gene not in gene_mut_map_pseudo:
                            gene_mut_map_pseudo[gene] = set()
                            gene_mut_map_pseudo[gene].add(aac)
                        else:
                            gene_mut_map_pseudo[gene].add(aac)
                    else:
                        if gene not in gene_mut_map_new:
                            gene_mut_map_new[gene] = set()
                            gene_mut_map_new[gene].add(aac)
                        else:
                            gene_mut_map_new[gene].add(aac)
        else:
            for aac in aaclist:
                aacbase = aac[:-1]
                if aacbase in twin_set:
                    if gene not in gene_mut_map_twins:
                        gene_mut_map_twins[gene] = set()
                        gene_mut_map_twins[gene].add(aac)
                    else:
                        gene_mut_map_twins[gene].add(aac)
                else:
                    twin_set.add(aacbase)
                    if gene not in gene_mut_map_new:
                        gene_mut_map_new[gene] = set()
                        gene_mut_map_new[gene].add(aac)
                    else:
                        gene_mut_map_new[gene].add(aac)

    #print gene_mut_map_pseudo

    t2 = time.time()

    value_strs = []
    for gene in gene_mut_map_new:
        if gene in stored_genes:
            gene_id = stored_genes[gene][0]
        else:
            gene_id = new_genes[gene]
        base_map = {}
        for aac in gene_mut_map_new[gene]:
            aacb = aac[:-1]
            if not aacb in base_map:
                base_map[aacb] = [aac[-1]]                       
            else:
                base_map[aacb].append(aac[-1])
        if gene in gene_mut_map_twins:
            for aac in gene_mut_map_twins[gene]:
                aacb = aac[:-1]
                if aacb in base_map:
                    base_map[aacb].append(aac[-1])
        for base in base_map:
            value_strs.append("('%s','%s%s','-')" % (str(gene_id),base,','.join(base_map[base])))

    t3 = time.time()

    if len(value_strs) > 0:

        size = len(value_strs)
        if size > 0:
            bound = 50000
            if size > bound:
                m = size/bound
                rest = size%bound
                if rest == 0:
                    part_size = size/m
                else:
                    m += 1
                    part_size = size/m
                parts = []
                for i in range(0,m-1):
                    parts.append(value_strs[(i*part_size):((i+1)*part_size)])
                parts.append(value_strs[(m-1)*part_size:])
            else:
                parts = [value_strs]
            for part in parts:

                sql = "INSERT IGNORE INTO Mutation(Gene,Amino_Acid_Change,Polyphen2_Score) VALUES %s" % ','.join(part)
                try:
                    cursor.execute(sql)
                    db.commit()
                except:
                    [e,f,g] = sys.exc_info()
                    raise NameError("Error in mutationCheck: %s,%s" % (sql,f))
    t4 = time.time()
    value_strs = []
    mut_ids = set()
    for gene in gene_mut_map_pseudo:
        gene_id = stored_genes[gene][0]
        aacb_mut_id_aacs_map = update_map[gene_id]
        base_map = {}
        for aac in gene_mut_map_pseudo[gene]:
            aacb = aac[:-1]
            if not aacb in base_map:
                base_map[aacb] = [aac[-1]]                       
            else:
                base_map[aacb].append(aac[-1])
        if gene in gene_mut_map_twins:
            for aac in gene_mut_map_twins[gene]:
                aacb = aac[:-1]
                if aacb in base_map:
                    base_map[aacb].append(aac[-1])
        for base in base_map:
            (mut_id,aacs) =  aacb_mut_id_aacs_map[base]
            new_aa = []
            for aa in base_map[base]:
                if aacs[1:].count(aa) < 1:
                    new_aa.append(aa)
            """
            #This clause controls the update of the update of the type of the substituted amino acid, since this information can be found in the mutation-session table and is not currently used in the annotation, we can comment it, since the updates (as in this state) can be very expensive in a large mutation table
            #FIX: add WHERE Mutation_Id BETWEEN min/max clause add the end of the update query
            #Alternativ: Think about changing Amino_Acid_Change in the Mutation table to aac_base and save the substitution type in the mutation_session table
            if len(new_aa) > 0:
                value_strs.append("WHEN '%s' THEN '%s,%s'" % (str(mut_id),aacs,','.join(new_aa)))
                mut_ids.add(str(mut_id))
            """
    t5 = time.time()
    if len(value_strs) > 0:

        size = len(value_strs)
        if size > 0:
            bound = 50000
            if size > bound:
                m = size/bound
                rest = size%bound
                if rest == 0:
                    part_size = size/m
                else:
                    m += 1
                    part_size = size/m
                parts = []
                for i in range(0,m-1):
                    parts.append(value_strs[(i*part_size):((i+1)*part_size)])
                parts.append(value_strs[(m-1)*part_size:])
            else:
                parts = [value_strs]
            for part in parts:
                sql = "UPDATE Mutation SET Amino_Acid_Change = CASE Mutation_Id %s ELSE Amino_Acid_Change END" % (" ".join(part))
                try:
                    cursor.execute(sql)
                    db.commit()
                except:
                    [e,f,g] = sys.exc_info()
                    raise NameError("Error in mutationCheck: %s,%s" % (sql,f))
    t6 = time.time()

    for g_id in new_gene_ids: #Note: the max and min values are indeed the ones from the beginning of the method!
        bin_number = g_id//binsize
        bins.add(bin_number)

        if max_g_id == None or int(g_id) > max_g_id:
            max_g_id = int(g_id)
        if min_g_id == None or int(g_id) < min_g_id:
            min_g_id = int(g_id)

    min_max_tuples = []
    if max_g_id - min_g_id > binsize:
        for bin_number in bins:
            min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize))
    else:
        min_max_tuples = [(min_g_id,max_g_id)]


    total_results = []
    for (run_min_g_id,run_max_g_id) in min_max_tuples:

        sql = "SELECT Amino_Acid_Change,Mutation_Id,Gene FROM Mutation WHERE Gene BETWEEN %s and %s" % (str(run_min_g_id),str(run_max_g_id))

        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in mutationCheck: %s,%s" % (sql,f))
        total_results += results

    t7 = time.time()


    update_map = {}
    for row in total_results:
        gene_id = row[2]
        if not (gene_id in stored_gene_ids or gene_id in new_gene_ids):
            continue
        aacs = row[0].split(",")
        aacbase = aacs[0][:-1]
        mutation_id = row[1]
        if not gene_id in update_map:
            update_map[gene_id] = {aacbase:mutation_id}
        else:
            update_map[gene_id][aacbase] = mutation_id

    #print update_map

    t8 = time.time()

    mut_new_aa_map = {}
    gene_mut_map_new_id = {}
    gene_mut_map_pseudo_id = {}
    tag_update_map = {}

    for gene in gene_aaclist_map:
        aaclist = gene_aaclist_map[gene]
        if gene in stored_genes:
            gene_id = stored_genes[gene][0]
        else:
            gene_id = new_genes[gene]       
        for aac in aaclist:
            aac_base = aac.split(',')[0][:-1]
            mut_id = update_map[gene_id][aac[:-1]]
            if not mut_id in mut_new_aa_map:
                mut_new_aa_map[mut_id] = set()
                mut_new_aa_map[mut_id].add(aac[-1])
            else:
                mut_new_aa_map[mut_id].add(aac[-1])
            if gene in tag_map:
                if aac in tag_map[gene]:
                    tag_update_map[mut_id] = tag_map[gene][aac]
            else:
                tag_update_map[mut_id] = None
            if gene in gene_mut_map_new:
                if not gene in gene_mut_map_new_id:
                    gene_mut_map_new_id[gene] = (gene_id,{aac_base:mut_id})
                else:
                    gene_mut_map_new_id[gene][1][aac_base] = mut_id
            elif gene in gene_mut_map_pseudo:
                if not gene in gene_mut_map_pseudo_id:
                    gene_mut_map_pseudo_id[gene] = (gene_id,{aac_base:mut_id})
                else:
                    gene_mut_map_pseudo_id[gene][1][aac_base] = mut_id

    t9 = time.time()

    value_strs = []
    for mut_id in mut_new_aa_map:
        if tag_update_map[mut_id] == None:
            value_strs.append("('%s','%s','%s',NULL)" % (str(database_session),str(mut_id),','.join(mut_new_aa_map[mut_id])))
        else:
            value_strs.append("('%s','%s','%s','%s')" % (str(database_session),str(mut_id),','.join(mut_new_aa_map[mut_id]),tag_update_map[mut_id]))

    process = Process(target = backgroundInsertMS,args=(value_strs,MS_db))
    process.start()

    t10 = time.time()

    print "Time for mutationCheck Part 1: %s" % (str(t1-t0))
    print "Time for mutationCheck Part 2: %s" % (str(t2-t1))
    print "Time for mutationCheck Part 3: %s" % (str(t3-t2))
    print "Time for mutationCheck Part 4: %s" % (str(t4-t3))
    print "Time for mutationCheck Part 5: %s" % (str(t5-t4))
    print "Time for mutationCheck Part 6: %s" % (str(t6-t5))
    print "Time for mutationCheck Part 7: %s" % (str(t7-t6))
    print "Time for mutationCheck Part 8: %s" % (str(t8-t7))
    print "Time for mutationCheck Part 9: %s" % (str(t9-t8))
    print "Time for mutationCheck Part 10: %s" % (str(t10-t9))

    return gene_mut_map_new_id,gene_mut_map_pseudo_id,process

def backgroundInsertMS(value_strs,db):
    cursor = db.cursor()
    size = len(value_strs)
    if size > 0:
        bound = 50000
        if size > bound:
            m = size/bound
            rest = size%bound
            if rest == 0:
                part_size = size/m
            else:
                m += 1
                part_size = size/m
            parts = []
            for i in range(0,m-1):
                parts.append(value_strs[(i*part_size):((i+1)*part_size)])
            parts.append(value_strs[(m-1)*part_size:])
        else:
            parts = [value_strs]
        for part in parts:

            sql = "INSERT IGNORE INTO RS_Mutation_Session(Session,Mutation,New_AA,Tag) VALUES %s" % ','.join(part)
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in mutationCheck: %s,%s" % (sql,f))

#called by serializedPipeline
def getTemplates(gene_id_list,db,cursor):
    gene_template_map = {}
    if len(gene_id_list) > 0:
        sql = "SELECT Template_Id,Name,Sequence_Identity,Alignment_Length,Target_Chain,Resolution,R_Value,Quality_Score,Original_Target_Chain,Original_Chains,Chains,Homooligomer,Gene  FROM Template"
             
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in getTemplates: %s,%s" % (sql,f))

        if results == ():
            return gene_template_map
        for row in results:
            template_id = row[0]
            pdb_id = row[1]
            seq_id = row[2]
            aln_length = row[3]
            target_chain = row[4]
            resolution = row[5]
            r_value = row[6]
            score = row[7]
            original_target_chain = row[8]
            original_chains = row[9]
            chains = row[10]
            oligos = row[11]
            gene_id = row[12]
            if not gene_id in gene_id_list:
                continue
            (sub_lig_dist,sub_ano) = (-0.1,0.0)
            template = [pdb_id,seq_id,target_chain,aln_length,resolution,[],r_value,score,sub_lig_dist,sub_ano,original_target_chain,original_chains]
            oligo_set = set([])
            for oligo in oligos:
                oligo_set.add(oligo)
            if not gene_id in gene_template_map:
                gene_template_map[gene_id] = ({template_id:template},{pdb_id:oligo_set})
            else:
                gene_template_map[gene_id][0][template_id] = template
                gene_template_map[gene_id][1][pdb_id] = oligo_set
    return gene_template_map

def getTemplatesFromIdList(template_id_list,db,cursor):
    t0 = time.time()
    template_ligand_dict = getLigandDict(template_id_list,db,cursor)
    t1 = time.time()
    
    sql = "SELECT Template_Id,Name,Sequence_Identity,Gene,Alignment_Length,Target_Chain,Resolution,R_Value,Quality_Score,Chains FROM Template WHERE Template_Id BETWEEN %s AND %s" % (str(min(template_id_list)),str(max(template_id_list)))
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("NameError in getTemplatesFromIdList: %s" % sql)
    template_results = {}
    t2 = time.time()
    for row in results:
        template_id = row[0]
        if template_id in template_id_list:
            pdb_id = row[1]
            seq_id = row[2]
            gene_id = row[3]
            aln_length = row[4]
            target_chain = row[5]
            resolution = row[6]
            r_value = row[7]
            score = row[8]
            chains = row[9]
            if template_id in template_ligand_dict:
                ligands = template_ligand_dict[template_id].values()
            else:
                ligands = []
            iaps = chains.split(";")
            for iap in iaps:
                iap_info = iap.split(":")
                if len(iap_info) > 1:
                    ligands.append([iap_info[1],iap_info[0]])

            template = [pdb_id,seq_id,target_chain,aln_length,resolution,ligands,r_value,score]
            template_results[template_id] = (template,gene_id,chains)
    t3 = time.time()
    time_1 = t1 - t0
    time_2 = t2 - t1
    time_3 = t3 - t2
    print("getTemplatesFromId time 1: " + str(time_1))
    print("getTemplatesFromId time 2: " + str(time_2))
    print("getTemplatesFromId time 3: " + str(time_3))
    return template_results

#called by babel
def getPDB(template_id,db,cursor):
    sql = "SELECT Name FROM Template WHERE Template_Id = '%s'" % str(template_id)
         
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        print("NameError in getPDB")
        db.rollback()

    return results[0][0]

#called by babel
def checkMutationSession(mutation_id,session_id,db,cursor):
    sql = "SELECT Mutation FROM RS_Mutation_Session WHERE Mutation = '%s' AND Session = '%s'" % (str(mutation_id),str(session_id))
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

#called by babel
def getGeneFromTemplate(template_id,db,cursor):
    sql = "SELECT Gene FROM Template \
       WHERE Template_Id = '%s'" % template_id
         
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getGeneFromTemplate: %s" % sql)
        db.rollback()
    #print(pdb_id)
    #print(results)
    if results == ():
        return 0
    return results[0][0]

#called by babel
def getAnnoIds(db,cursor,template_id=0,mutation_id=0):
    if template_id == 0:
        if mutation_id == 0:
            raise NameError("Invalid Input in getAnnoIds")
        sql = ("""SELECT Template FROM RS_Mutation_Template WHERE Mutation = '%s'""" % str(mutation_id))
    elif mutation_id == 0:
        sql = ("""SELECT Mutation FROM RS_Mutation_Template WHERE Template = '%s'""" % str(template_id))
    else:
        raise NameError("Invalid Input in getAnnoIds")
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getAnnoIds: %s" % sql)
        db.rollback()
    ids = []
    for row in results:
        ids.append(row[0])
    return ids

def getMutationDict(mutation_id_list,db,cursor):

    min_mut_id = min(mutation_id_list)
    max_mut_id = max(mutation_id_list)

    sql = "SELECT Mutation_Id,Amino_Acid_Change,Gene,Polyphen2_Score FROM Mutation WHERE Mutation_Id BETWEEN %s AND %s" % (str(min_mut_id),str(max_mut_id))
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("NameError in getMutationDict: %s" % sql)
    mutation_dict= {}
    for row in results:
        mut_id = row[0]
        if mut_id in mutation_id_list:
            mutation_dict[mut_id] = (row[1],row[2],row[3])
    return mutation_dict

#called by serializedPipeline
def updateSession(session_id,time,db,cursor):
    sql = "UPDATE IGNORE Session SET End = '%s' WHERE Session_Id = '%s'" % (str(time),str(session_id))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Commit your changes in the database
        db.commit()
    except:
        print("Couldn't update Session '%s'") % (str(session_id))
        # Rollback in case there is any NameError
        db.rollback()

#called by serializedPipeline
def insertSession(time,ori_filename,seq_thresh,res_thresh,cov_thresh,db,cursor):
    sql = "INSERT IGNORE INTO Session(Input_File,Start,Sequence_Identity_Threshold,Coverage_Threshold,Resolution_Threshold) VALUES ('%s','%s','%s','%s','%s')" % (ori_filename,str(time),str(seq_thresh),str(cov_thresh),str(res_thresh))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Couldn't insert Session '%s'") % (sql)
        # Rollback in case there is any NameError
        db.rollback()
    session_id = getLastAutoInc(db,cursor)
    return session_id

#called by output
#called by serializedPipeline
#called by structman
def getSessionId(infile,db,cursor):
    sql = "SELECT Session_Id FROM Session WHERE Input_File = '%s'" % infile
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getSessionId: %s" % sql)
        # Rollback in case there is any NameError
    if results == ():
        return None
    try:
        sid = results[0][0]
    except:
        raise NameError("Input_File not in Session: %s" % sql)
    return results[0][0]

#called by babel
def getSessionFile(session_id,db,cursor):
    sql = "SELECT Input_File FROM Session WHERE Session_Id = '%s'" % str(session_id)
    #print sql
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getSessionFile: '%s'" % sql)
        # Rollback in case there is any NameError
    #print results
    if results == ():
        return None
    return results[0][0]

#called by serializedPipeline
def addGeneInfos(gene_info_map,db,cursor):
    #structure of gene_info_map: {gene_id:(refseqs,{go_term_id:go_term_name},{reactome_id:pathway_name}),sequence}
    ref_value_strs = []
    go_term_ids = set()
    reac_ids = set()
    seq_value_strs = []
    for gene_id in gene_info_map:
        (refseqs,go_terms,pathways,sequence) = gene_info_map[gene_id]
        for go_id in go_terms:
            go_term_ids.add(go_id)
        for reac_id in pathways:
            reac_ids.add(reac_id)
        ref_value_strs.append("WHEN '%s' THEN '%s'" % (str(gene_id),refseqs))
        seq_value_strs.append("WHEN '%s' THEN '%s'" % (str(gene_id),sequence))

    stored_go_terms = {}
    if len(go_term_ids) > 0:
        sql = "SELECT GO_Term_Id,Id FROM GO_Term WHERE Id IN ('%s')" % "','".join(go_term_ids)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in addGeneInfos: %s,%s" % (sql,f))
        for row in results:
            stored_go_terms[str(row[1])] = row[0]

    stored_pathways = {}
    if len(reac_ids) > 0:
        sql = "SELECT Pathway_Id,Reactome_Id FROM Pathway WHERE Reactome_Id IN ('%s')" % "','".join(reac_ids)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in addGeneInfos: %s,%s" % (sql,f))
        for row in results:
            stored_pathways[str(row[1])] = row[0]

    new_go_terms = {}
    go_value_strs = []
    new_pathways = {}
    pathway_value_strs = []
    for gene_id in gene_info_map:
        (refseqs,go_terms,pathways,sequence) = gene_info_map[gene_id]
        for go_id in go_terms:
            if go_id not in stored_go_terms:
                new_go_terms[go_id] = go_terms[go_id]
                go_value_strs.append("('%s','%s')" % (go_terms[go_id].replace("'","''"),go_id))
        for reac_id in pathways:
            if reac_id not in stored_pathways:
                new_pathways[reac_id] = pathways[reac_id]
                pathway_value_strs.append("('%s','%s')" % (pathways[reac_id].replace("'","''"),reac_id))

    if len(go_value_strs) > 0:
        sql = "INSERT IGNORE INTO GO_Term(Name,Id) VALUES %s" % ','.join(go_value_strs)
        try:
            cursor.execute(sql)
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in addGeneInfos: %s,%s" % (sql,f))

    if len(pathway_value_strs) > 0:
        sql = "INSERT IGNORE INTO Pathway(Name,Reactome_Id) VALUES %s" % ','.join(pathway_value_strs)
        try:
            cursor.execute(sql)
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in addGeneInfos: %s,%s" % (sql,f))
    
    stored_go_terms = {}
    if len(go_term_ids) > 0:
        sql = "SELECT GO_Term_Id,Id FROM GO_Term WHERE Id IN ('%s')" % "','".join(go_term_ids)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in addGeneInfos: %s,%s" % (sql,f))
        for row in results:
            stored_go_terms[str(row[1])] = row[0]

    stored_pathways = {}
    if len(reac_ids) > 0:
        sql = "SELECT Pathway_Id,Reactome_Id FROM Pathway WHERE Reactome_Id IN ('%s')" % "','".join(reac_ids)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in addGeneInfos: %s,%s" % (sql,f))
        for row in results:
            stored_pathways[str(row[1])] = row[0]
    
    go_value_strs = []
    pathway_value_strs = []
    for gene_id in gene_info_map:
        (refseqs,go_terms,pathways,sequence) = gene_info_map[gene_id]
        for go_id in go_terms:
            go_value_strs.append("('%s','%s')" % (str(gene_id),str(stored_go_terms[go_id])))
        for reac_id in pathways:
            pathway_value_strs.append("('%s','%s')" % (str(gene_id),str(stored_pathways[reac_id])))
    
    if len(go_value_strs) > 0:
        sql = "INSERT IGNORE INTO RS_Gene_GO_Term(Gene,GO_Term) VALUES %s" % ','.join(go_value_strs)
        try:
            cursor.execute(sql)
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in addGeneInfos: %s,%s" % (sql,f))

    if len(pathway_value_strs) > 0:
        sql = "INSERT IGNORE INTO RS_Gene_Pathway(Gene,Pathway) VALUES %s" % ','.join(pathway_value_strs)
        try:
            cursor.execute(sql)
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in addGeneInfos: %s,%s" % (sql,f))

    gene_id_strs = [str(x) for x in gene_info_map.keys()]

    if len(ref_value_strs) > 0:
        sql = "UPDATE Gene SET Genbank_Protein_Accession_Number = CASE Gene_Id %s ELSE Genbank_Protein_Accession_Number END WHERE Gene_Id IN (%s)" % (" ".join(ref_value_strs),','.join(gene_id_strs))
        try:
            cursor.execute(sql)
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in addGeneInfos: %s,%s" % (sql,f))

    if len(seq_value_strs) > 0:
        sql = "UPDATE Gene SET Sequence = CASE Gene_Id %s ELSE Sequence END WHERE Gene_Id IN (%s)" % (" ".join(seq_value_strs),','.join(gene_id_strs))
        try:
            cursor.execute(sql)
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in addGeneInfos: %s,%s" % (sql,f))

    return

#called by serializedPipeline
def getSequences(stored_genes,sequence_map,db,cursor):
    gene_id_gene_map = {}
    stored_gene_ids = set([])
    for gene in stored_genes:
        stored_gene_ids.add(stored_genes[gene][0])
        gene_id_gene_map[stored_genes[gene][0]] = gene
    gene_id_sequence_map = {}
    if len(stored_gene_ids) > 0:
        sql = "SELECT Gene_Id,Sequence FROM Gene"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in getSequences: %s,%s" % (sql,f))
        for row in results:
            gene_id = row[0]
            if not gene_id in stored_gene_ids:
                continue
            gene_id_sequence_map[gene_id] = row[1]

    for gene_id in gene_id_sequence_map:
        sequence_map[gene_id_gene_map[gene_id]] = gene_id_sequence_map[gene_id]

    return sequence_map

#called by serializedPipeline
def addErrorCodeToGene(gene_error_map,db,cursor):
    #structure of gene_error_map: {gene_id:error-id}
    error = {0:"Illegal uniprot ID",1:"Illegal uniprot ID or uniprot connection error",2:"Received empty wt-sequence",3:"Found no usable template",4:"Unknown Error"}
    value_strs_c = []
    value_strs_e = []
    for gene_id in gene_error_map:
        error_id = gene_error_map[gene_id]
        value_strs_c.append("WHEN %s THEN %d" % (str(gene_id),error_id))
        value_strs_e.append("WHEN %s THEN '%s'" % (str(gene_id),error[error_id]))
    if len(value_strs_c) > 0:
        gene_id_string = ','.join([str(x) for x in gene_error_map.keys()])
        sql = "UPDATE IGNORE Gene SET Error_Code = CASE Gene_Id %s ELSE Error_Code END , ERROR = CASE Gene_Id %s ELSE Error END WHERE Gene_Id IN (%s)" % (" ".join(value_strs_c)," ".join(value_strs_e),gene_id_string)
        try:
            cursor.execute(sql)
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Couldn't update Gene-Error: %s,%s" % (sql,f))
    return

def backgroundInsertAS(value_strings,max_m,min_m,db,anno_anno_calculation):
    cursor = db.cursor()
    size = len(value_strings)
    if size > 0:
        bound = 2000
        if size > bound:
            m = size/bound
            rest = size%bound
            if rest == 0:
                part_size = size/m
            else:
                m += 1
                part_size = size/m
            parts = []
            for i in range(0,m-1):
                parts.append(value_strings[(i*part_size):((i+1)*part_size)])
            parts.append(value_strings[(m-1)*part_size:])
        else:
            parts = [value_strings]
        for part in parts:
            sql = "INSERT IGNORE INTO RS_Annotation_Session(Session,Mutation,Template) VALUES %s;" % ",".join(part)
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in MultipleAnnotationSession: %s,%s" % (sql,f))
    else:
        return {}

    if anno_anno_calculation:

        sql = 'SELECT Template_1,Template_2,Mutation_1,Mutation_2,Chain_1,Chain_2,Distance,Atompair FROM RS_Annotation_Annotation WHERE Mutation_1 <= %s AND Mutation_1 >= %s' % (str(max_m),str(min_m))
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in backgroundInsertAS: %s\n%s" % (sql,f))

        anno_anno_value_strs = []
        dupl_set = set()
        for row in results:
            t_1 = row[0]
            t_2 = row[1]
            m_1 = row[2]
            m_2 = row[3]
            c_1 = row[4]
            c_2 = row[5]
            d = str(row[6])
            ap = row[7]
            if t_1 in t_m_map and t_2 in t_m_map:
                if m_1 in t_m_map[t_1] and m_2 in t_m_map[t_2]:
                    if not (m_1,m_2) in dupl_set:
                        anno_anno_value_strs.append("('%s','%s','%s','%s','%s','%s','%s','%s','%s')" % (str(t_1),str(t_2),str(m_1),str(m_2),c_1,c_2,str(session_id),d,ap))
                        dupl_set.add((m_1,m_2))
        size = len(anno_anno_value_strs)
        if size > 0:
            bound = 2000
            if size > bound:
                m = size/bound
                rest = size%bound
                if rest == 0:
                    part_size = size/m
                else:
                    m += 1
                    part_size = size/m
                parts = []
                for i in range(0,m-1):
                    parts.append(anno_anno_value_strs[(i*part_size):((i+1)*part_size)])
                parts.append(anno_anno_value_strs[(m-1)*part_size:])
            else:
                parts = [anno_anno_value_strs]
            for part in parts:
                sql = "INSERT IGNORE INTO RS_Annotation_Annotation(Template_1,Template_2,Mutation_1,Mutation_2,Chain_1,Chain_2,Session,Distance,Atompair) VALUES %s" % ','.join(part)
                try:
                    cursor.execute(sql)
                    db.commit()
                except:
                    [e,f,g] = sys.exc_info()
                    raise NameError("Error in MultipleAnnotationSession: %s,%s" % (sql,f))

#called by serializedPipeline
def insertMultipleAnnotationSession(session_id,pseudo_template_map,gene_mut_map_pseudo,db,cursor,AS_db,anno_session_mapping = True,anno_anno_calculation=False):
    #structure of pseudo_template_map: {Uniprot-Id:(template-list,{template-id:stored-template},oligo_map)}
    #structure of gene_mut_map_pseudo: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    #print pseudo_template_map
    #print gene_mut_map_pseudo
    value_strings = []
    stored_annotations = {}
    t_m_map = {}
    max_m = None
    min_m = None
    for gene in pseudo_template_map:
        if gene in gene_mut_map_pseudo:
            (gene_id,aac_map) = gene_mut_map_pseudo[gene]
            pseudo_templates = pseudo_template_map[gene][1]

            for aac_base in aac_map:
                mut_id = aac_map[aac_base]
                if max_m == None or int(mut_id) > max_m:
                    max_m = int(mut_id)
                if min_m == None or int(mut_id) < min_m:
                    min_m = int(mut_id)
                for template_id in pseudo_templates:
                    value_strings.append("('%s','%s','%s')" % (str(session_id),str(mut_id),str(template_id)))
                    if not template_id in stored_annotations:
                        stored_annotations[template_id] = {}
                        t_m_map[template_id] = set()
                    stored_annotations[template_id][aac_base] = mut_id
                    t_m_map[template_id].add(mut_id)

    if not anno_session_mapping:
        return stored_annotations,None

    proc = Process(target=backgroundInsertAS, args=(value_strings,max_m,min_m,AS_db,anno_anno_calculation))
    proc.start()

    return stored_annotations,proc

#called by serializedPipeline
def insertTemplates(template_list,session_id,db,cursor,cwd,pdb_path,smiles_path,inchi_path):
    t0 = time.time()
    #structure of template_list: [(gene_id,template,oligo_counter,alignment)]
    value_strs = []
    gene_id_list = set()
    ligand_map = {}
    ligands = set()
    for (gene_id,template,oligos,alignment) in template_list:
        chains = []
        for iap in template[5]:
            if iap[0] != "Ligand":
                #print template
                #print iap
                chains.append("%s:%s" % (iap[1],iap[0]))
        chains = ";".join(chains)
        oligos = ''.join(oligos)

        value_strs.append("('%s','%d','%1.3f','%1.2f','%s','%1.4f','%1.4f','%f','%s','%s','%s','%s','%s')" % (template[0],int(gene_id),float(template[1]),float(template[3]),template[2],float(template[4]),float(template[6]),float(template[7]),template[10],template[11],chains,alignment.replace("'","\\'"),oligos))
        gene_id_list.add(gene_id)
        if not gene_id in ligand_map:
            ligand_map[gene_id] = {}
        for iap in template[5]:
            ia_type = iap[0]
            if ia_type == "Ligand":
                if not template[0] in ligand_map[gene_id]:
                    ligand_map[gene_id][template[0]] = [iap]
                else:
                    ligand_map[gene_id][template[0]].append(iap)
                ligands.add(iap[1])
    t1 = time.time()
    if len(value_strs) > 0:
        #For large amounts of templates, a single insert can destroy the connection to the database => Solution: divide large inserts into several smaller ones
        value_strs_list = []
        i = 0
        part_size = 500
        while i*part_size < len(value_strs):
            a = i*part_size
            b = min(((i+1)*part_size,len(value_strs)))
            i += 1
            sql = "INSERT IGNORE INTO Template(Name,Gene,Sequence_Identity,Alignment_Length,Target_Chain,Resolution,R_Value,Quality_Score,Original_Target_Chain,Original_Chains,Chains,Alignment,Homooligomer) VALUES %s" % ','.join(value_strs[a:b])
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in insertTemplates: %s" % (f))

    t2 = time.time()
    template_id_map = {}
    if len(gene_id_list) > 0:
        sql = "SELECT Template_Id,Name,Gene FROM Template"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in insertTemplates: %s\n%s" % (sql,f))

        for row in results:
            gene_id = row[2]
            if not gene_id in gene_id_list:
                continue
            template_id = row[0]
            pdb_id = row[1]
            if not gene_id in template_id_map:
                template_id_map[gene_id] = {pdb_id:template_id}
            else:
                template_id_map[gene_id][pdb_id] = template_id
    t3 = time.time()
    stored_ligands = {}        
    if len(ligands) > 0:
        sql = "SELECT Ligand_Id,Name FROM Ligand"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in insertTemplates: %s\n%s" % (sql,f))
        for row in results:
            ligand_name = row[1]
            if ligand_name not in ligands:
                continue
            stored_ligands[ligand_name] = row[0]


    t4 = time.time()
    value_strs = []
    new_ligands = set()

    ligand_db = pdb.parseLigandDB(smiles_path,inchi_path)
    #print ligand_db

    #print ligand_db
    t5 = time.time()

    #print ligand_map
    for gene_id in ligand_map:
        for pdb_id in ligand_map[gene_id]:
            iaps = ligand_map[gene_id][pdb_id]
            for iap in iaps:
                name = iap[1]
                if not (name in stored_ligands or name in new_ligands):
                    if not (name == "UNK" or name == "UNX"):
                        new_ligands.add(name)
                        res = iap[2]
                        chain = iap[3]
                        if name in ligand_db:
                            (smiles,inchi) = ligand_db[name]
                        elif len(name) < 3:
                            if len(name) == 1:
                                smiles = "[%s]" % name
                            else:
                                smiles = "[%s%s]" % (name[0],name[1].lower())
                            inchi = "InChI=1S/%s" % smiles
                        else:
                            (smiles,inchi) = pdb.getSI(pdb_id,name,res,chain,cwd,pdb_path)
                        value_strs.append("('%s','%s','%s')" % (name,smiles,inchi))


    t6 = time.time()

    if len(value_strs) > 0:
        sql = "INSERT IGNORE INTO Ligand(Name,Smiles,Inchi) VALUES %s " % ','.join(value_strs)
        try:
            cursor.execute(sql)
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in insertTemplates: %s\n%s" % (sql,f))

    t7 = time.time()

    new_ligand_map = {}
    if len(new_ligands) > 0:
        sql = "SELECT Ligand_Id,Name FROM Ligand"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in insertTemplates: %s\n%s" % (sql,f))
        for row in results:
            ligand_name = row[1]
            if not ligand_name in new_ligands:
                continue
            new_ligand_map[ligand_name] = row[0]

    t8 = time.time()

    value_strs = []

    for (gene_id,template,oligo_counter,alignment_pir) in template_list:
        pdb_id = template[0]
        template_id = template_id_map[gene_id][pdb_id]
        for iap in template[5]:
            ia_type = iap[0]
            if ia_type == "Ligand":
                name = iap[1]
                if name in stored_ligands:
                    ligand_id = stored_ligands[name]
                elif name in new_ligand_map:
                    ligand_id = new_ligand_map[name]
                else:
                    continue
                res = iap[2]
                chain = iap[3]
                value_strs.append("('%s','%s','%s','%s')" % (str(ligand_id),str(template_id),str(chain),str(res)))

    t9 = time.time()

    if len(value_strs) > 0:
        sql = "INSERT IGNORE INTO RS_Ligand_Template(Ligand,Template,Chain,Residue) VALUES %s" % ','.join(value_strs)
        try:
            cursor.execute(sql)
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in insertTemplates: %s\n%s" % (sql,f))

    t10 = time.time()

    """
    print "insertTemplates part 1: ",t1-t0
    print "insertTemplates part 2: ",t2-t1
    print "insertTemplates part 3: ",t3-t2
    print "insertTemplates part 4: ",t4-t3
    print "insertTemplates part 5: ",t5-t4
    print "insertTemplates part 6: ",t6-t5
    print "insertTemplates part 7: ",t7-t6
    print "insertTemplates part 8: ",t8-t7
    print "insertTemplates part 9: ",t9-t8
    print "insertTemplates part 10: ",t10-t9
    """
    return template_id_map

#called by serializedPipeline
def insertAnnotations(template_annotation_map,anno_anno_dict,template_mutation_map,lig_wf,chain_wf,session,db,cursor,full_anno_anno_calculation=False,error_annotations_into_db=True,anno_session_mapping = True):
    #structure of template_mutation_map: {pdb_id:[sub_info_map,{template_id:{aac_base:mutation_id}},template]}
    #structure of anno_anno_dict: {(smaller template_id,bigger template_id,aac_base of first template,aac_base of second template):(min_dist,atompair)}
    #print template_annotation_map
    #print template_mutation_map
    error = {0:"Mutation mapped to a gap in the target-template alignment",1:"Mutation position outside of the chain (This should not happen)"}
    value_strs = []
    session_value_strs = []
    anno_anno_value_strs = []
    for template_id in template_annotation_map:
        [pdb_id,annotations] = template_annotation_map[template_id]
        sub_infos = template_mutation_map[pdb_id][0][template_id][0]
        template = template_mutation_map[pdb_id][2]
        template_id_aac_map = template_mutation_map[pdb_id][1]
        aac_bases = set()
        for aac_base in template_id_aac_map[template_id]:
            if not aac_base in aac_bases:
                aac_bases.add(aac_base)
                #print template_id,pdb_id
                if not aac_base in annotations:
                    continue
                annotation = annotations[aac_base]
                mutation_id = template_id_aac_map[template_id][aac_base]
                if annotation == 0:
                    if error_annotations_into_db:
                        value_strs.append("('%s','%s',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'0','%s',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)" % (str(mutation_id),str(template_id),error[0]))
                        session_value_strs.append("('%s','%s','%s')" % (str(mutation_id),str(template_id),str(session)))
                elif annotation == 1:
                    if error_annotations_into_db:
                        value_strs.append("('%s','%s',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,'1','%s',NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)" % (str(mutation_id),str(template_id),error[1]))
                        session_value_strs.append("('%s','%s','%s')" % (str(mutation_id),str(template_id),str(session)))
                else:
                    lig_dists = annotation[0]
                    chain_dists = annotation[1]
                    rsa = annotation[2]
                    res_id = str(annotation[3])
                    ssa = annotation[4]
                    homomer_map = annotation[5]
                    profile = annotation[6]
                    if profile == None:
                        profile_str = "NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL"
                    else:
                        profile_str = "'%s','%1.4f','%s','%1.4f','%s','%1.4f','%s','%1.4f','%s','%1.4f'" % (profile['ligand'][0],profile['ligand'][1],profile['interchain'][0],profile['interchain'][1],profile['neighbor'][0],profile['neighbor'][1],profile['short'][0],profile['short'][1],profile['long'][0],profile['long'][1])

                    sub_info = sub_infos[aac_base]
                    substitution_pos_template = sub_info[0]
                    substitution_res_template = sub_info[1]
                    lig_dist_strings = []
                    dists = []
                    for lig_id in lig_dists:
                        (dist,atom_pair) = lig_dists[lig_id]
                        if not dist == None:
                            dists.append(dist)
                            lig_dist_strings.append("%s:%s:(%s-%s)" % (str(lig_id),dist,str(atom_pair[0]),str(atom_pair[1])))
                    if len(dists) > 0:
                        min_lig_dist = min(dists)
                    else:
                        min_lig_dist = -1
                    dists = []
                    lig_dist_string = ",".join(lig_dist_strings)
                    chain_dist_strings = []
                    for chain_id in chain_dists:
                        (dist,atom_pair,min_resi) = chain_dists[chain_id]
                        if not dist == None:
                            dists.append(dist)
                            chain_dist_strings.append("%s.%s:%s:(%s-%s)" % (chain_id,min_resi,dist,str(atom_pair[0]),str(atom_pair[1])))
                    if len(dists) > 0:
                        min_chain_dist = min(dists)
                    else:
                        min_chain_dist = -1
                    chain_dist_string = ",".join(chain_dist_strings)
                    candidate = templateFiltering.candidateScore(min_lig_dist,min_chain_dist,lig_wf=lig_wf,chain_wf=chain_wf)
                    if len(lig_dists) == 0:
                        lig_dist_string = "-"
                    if len(chain_dists) == 0:
                        chain_dist_string = "-"
                    res_id = res_id.replace("'","\\'")
                    if rsa == None:
                        rsa_str = "NULL"
                    else:
                        rsa_str = "'%1.2f'" % float(rsa)
                    if ssa == None:
                        ssa = "NULL"
                    else:
                        ssa = "'%s'" % ssa
                    homo_strs = []
                    for homo_chain in homomer_map:
                        (min_d,atom,atom2) = homomer_map[homo_chain]
                        homo_strs.append('%s:%1.2f(%s-%s)' % (homo_chain,min_d,atom,atom2))
                    homo_str = ','.join(homo_strs)
                    value_strs.append("('%s','%s','%s','%s',%s,%s,'%s','%s','%1.2f','%s',NULL,NULL,%s)" % (str(mutation_id),str(template_id),lig_dist_string,chain_dist_string,rsa_str,ssa,str(res_id),substitution_res_template,float(candidate),homo_str,profile_str))
                    session_value_strs.append("('%s','%s','%s')" % (str(mutation_id),str(template_id),str(session)))
                       
    for (m_id,m_id_2) in anno_anno_dict:
        if not full_anno_anno_calculation:
            (template_id,template_id_2,min_dist,atom,atom2,chain,chain_2) = anno_anno_dict[(m_id,m_id_2)]

            anno_anno_value_strs.append("('%s','%s','%s','%s','%s','%s','%s','%s','%s-%s')" % (str(template_id),str(template_id_2),str(m_id),str(m_id_2),chain,chain_2,str(session),str(min_dist),atom,atom2))
        else:
            for (template_id,template_id_2,min_dist,atom,atom2,chain,chain_2) in anno_anno_dict[(m_id,m_id_2)]:
                anno_anno_value_strs.append("('%s','%s','%s','%s','%s','%s','%s','%s','%s-%s')" % (str(template_id),str(template_id_2),str(m_id),str(m_id_2),chain,chain_2,str(session),str(min_dist),atom,atom2))
                                        

    size = len(value_strs)
    #print size
    if size > 0:
        bound = 2000
        if size > bound:
            m = size/bound
            rest = size%bound
            if rest == 0:
                part_size = size/m
            else:
                m += 1
                part_size = size/m
            parts = []
            for i in range(0,m-1):
                parts.append(value_strs[(i*part_size):((i+1)*part_size)])
            parts.append(value_strs[(m-1)*part_size:])
        else:
            parts = [value_strs]
        for part in parts:
            sql = """INSERT IGNORE INTO RS_Mutation_Template
            (Mutation,
            Template,
            Sub_Lig_Dist,
            Sub_Chain_Distances,
            Relative_Surface_Access,
            Secondary_Structure_Assignment,
            Residue_Id,Amino_Acid,
            Candidate_Score,
            Homomer_Distances,
            Error_Code,
            Error,
            Ligand_Interaction_Degree,
            Ligand_Interaction_Score,
            Chain_Interaction_Degree,
            Chain_Interaction_Score,
            Short_Interaction_Degree,
            Short_Interaction_Score,
            Medium_Interaction_Degree,
            Medium_Interaction_Score,
            Long_Interaction_Degree,
            Long_Interaction_Score
            ) VALUES %s""" % ','.join(part)
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in insertAnnotations: %s" % (f))

    if anno_session_mapping:
        size = len(session_value_strs)
        if size > 0:
            bound = 100000
            if size > bound:
                m = size/bound
                rest = size%bound
                if rest == 0:
                    part_size = size/m
                else:
                    m += 1
                    part_size = size/m
                parts = []
                for i in range(0,m-1):
                    parts.append(session_value_strs[(i*part_size):((i+1)*part_size)])
                parts.append(session_value_strs[(m-1)*part_size:])
            else:
                parts = [session_value_strs]
            for part in parts:
                sql = "INSERT IGNORE INTO RS_Annotation_Session(Mutation,Template,Session) VALUES %s" % ','.join(part)
                try:
                    cursor.execute(sql)
                    db.commit()
                except:
                    [e,f,g] = sys.exc_info()
                    raise NameError("Error in insertAnnotations: %s,\n%s" % (sql,f))

    size = len(anno_anno_value_strs)
    if size > 0:
        bound = 100000
        if size > bound:
            m = size/bound
            rest = size%bound
            if rest == 0:
                part_size = size/m
            else:
                m += 1
                part_size = size/m
            parts = []
            for i in range(0,m-1):
                parts.append(anno_anno_value_strs[(i*part_size):((i+1)*part_size)])
            parts.append(anno_anno_value_strs[(m-1)*part_size:])
        else:
            parts = [anno_anno_value_strs]
        for part in parts:
            sql = "INSERT IGNORE INTO RS_Annotation_Annotation(Template_1,Template_2,Mutation_1,Mutation_2,Chain_1,Chain_2,Session,Distance,Atompair) VALUES %s" % ','.join(part)
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in insertAnnotations: %s,%s" % (sql,f))


#called by serializedPipeline
def getAlignments(g_t_a_map,db,cursor):
    t0 = time.time()

    template_id_list = set()
    min_t = None
    max_t = None
    for gene_id in g_t_a_map:
        for t_id in g_t_a_map[gene_id].keys():
            t = int(t_id)
            if min_t == None or t < min_t:
                min_t = t
            if max_t == None or t > max_t:
                max_t = t
            template_id_list.add(t_id) 

    gene_template_alignment_map = {}

    pdb_map = {}
    
    max_diff = 20000
    if len(template_id_list) > 0:
        total_max = max_t
        while max_t - min_t > max_diff:            
            max_t = max_t - max_diff
 
        all_done = False
        while not all_done:
            sql = "SELECT Template_Id,Name,Sequence_Identity,Gene,Alignment_Length,Alignment,Target_Chain FROM Template WHERE Template_Id BETWEEN %s AND %s" % (str(min_t),str(max_t))
            try:
                cursor.execute(sql)
                results = cursor.fetchall()
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in getAlignments: %s,\n%s" % (sql,f))

            #print "Time for Select Alignments from Database: %s" % (str(t2-t))
            
            if max_t == total_max:
                all_done = True
            else:
                min_t = max_t + 1
                max_t = min(max_t + max_diff,total_max)

            for row in results:

                template_id = row[0]
                if not template_id in template_id_list:
                    continue
                pdb_id =row[1]
                seq_id = row[2]
                gene_id = row[3]
                aln_length = row[4]
                alignment = row[5]
                chain = row[6]

                lines = alignment.split("\n")
                #print(lines)
                nextline = False
                target_start = False
                template_start = False
                target_lines = []
                template_lines = []
                target_name = ""
                for line in lines:
                    if len(line) > 0:
                        if line[0] == ">":
                            ids = line.split(";")
                            if target_name == "":
                                target_name = ids[1]
                                target_name = target_name.replace(" ","")
                                target_name = target_name.replace("\n","")
                            nextline = True
                        elif nextline:
                            if not target_start:
                                target_start = True
                            else:
                                target_start = False
                                template_start = True
                                words = line.split(":")
                                startres = words[2]
                                endres = words[4]
                                chain = words[3]
                            nextline = False
                        elif line[0] == "\n":
                            template_start = False
                        elif target_start:
                            target_lines.append(line)
                        elif template_start:
                            template_lines.append(line)

                target_seq = "".join(target_lines)
                target_seq = target_seq.replace("*","")
                template_seq = "".join(template_lines)
                template_seq = template_seq.replace("*","")
                aaclist = g_t_a_map[gene_id][template_id]

                if not pdb_id in pdb_map:
                    pdb_map[pdb_id] = {}
                if not chain in pdb_map[pdb_id]:
                    pdb_map[pdb_id][chain] = []
                pdb_map[pdb_id][chain].append((target_seq,template_seq,aaclist,gene_id,aln_length,seq_id,template_id))

    t1 = time.time()
    print "Time for part 1 in getAlignments: %s" % (str(t1-t0))
    return pdb_map

def checkLigand(name,db,cursor):
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

#called by Babel
def getLigandTemplates(name,db,cursor):
    ligand_id = checkLigand(name,db,cursor)
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

#called by babel
def createLigandDB(outfile,session_id,db,cursor):
    sql = "SELECT Template,Mutation FROM RS_Annotation_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in createLigandDB: %s,\n%s" % (sql,f))
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
            [e,f,g] = sys.exc_info()
            raise NameError("Error in createLigandDB: %s,\n%s" % (sql,f))

        
        for row in results:
            if not row[1] in template_ids:
                continue
            if row[0] in mutation_ids:
                filtered_template_ids.add(row[1])

    #print filtered_template_ids

    ligand_ids = set()
    if len(filtered_template_ids) > 0:
        sql = "SELECT Ligand,Template FROM RS_Ligand_Template"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in createLigandDB: %s,\n%s" % (sql,f))

        for row in results:
            if not row[1] in filtered_template_ids:
                continue
            ligand_ids.add(row[0])

    #print ligand_ids

    lines = []
    if len(ligand_ids) > 0:
        sql = "SELECT Name,Smiles,Ligand_Id FROM Ligand"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in createLigandDB: %s,\n%s" % (sql,f))

        for row in results:
            if not row[2] in ligand_ids:
                continue
            lines.append("%s\t%s" % (row[1],row[0]))

    page = "\n".join(lines)
    #print page
    f = open(outfile, "wb")
    f.write(page)
    f.close()                

def getClassWT(lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc,surface_t_id,chain_t_id,lig_t_id):
    surface_threshold = 0.16
    if rel_sur_acc != None:
        if rel_sur_acc > 1.0:
            rel_sur_acc = 0.01*rel_sur_acc
    if lig_sub_dist == None:
        lig_sub_dist = "NONE"
    if chain_sub_dist == None:
        chain_sub_dist = "NONE"
    t_id = None
    try:
        mut_class = ""
        if (lig_sub_dist == "NONE" and chain_sub_dist == "NONE"):
            if rel_sur_acc == None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > surface_threshold:
                mut_class = "Surface isolated chain"
                t_id = surface_t_id
            else:
                mut_class = "Core isolated chain"
                t_id = surface_t_id
        elif lig_sub_dist == "NONE":
            if float(chain_sub_dist) > 5.0:
                if rel_sur_acc == None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > surface_threshold:
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
                if rel_sur_acc == None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > surface_threshold:
                    mut_class = "Surface"
                    t_id = surface_t_id
                else:
                    mut_class = "Core"
                    t_id = surface_t_id
            else:
                mut_class = "Contact Ligand"
                t_id = lig_t_id
                
                
        elif (float(lig_sub_dist) > 5.0 and float(chain_sub_dist) > 5.0):
            if rel_sur_acc == None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > surface_threshold:
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
        print "Error in class definition"
        print lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc
        print sys.exc_info()
    if mut_class == "":
        print lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc
        
    return mut_class,t_id

def getClass(lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc):
    surface_threshold = 0.16
    if rel_sur_acc != None:
        if rel_sur_acc > 1.0:
            rel_sur_acc = 0.01*rel_sur_acc
    if lig_sub_dist == None:
        lig_sub_dist = "NONE"
    if chain_sub_dist == None:
        chain_sub_dist = "NONE"
    try:
        mut_class = ""
        if (lig_sub_dist == "NONE" and chain_sub_dist == "NONE"):
            if rel_sur_acc == None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > surface_threshold:
                mut_class = "Surface isolated chain"
            else:
                mut_class = "Core isolated chain"
        elif lig_sub_dist == "NONE":
            if float(chain_sub_dist) > 5.0:
                if rel_sur_acc == None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > surface_threshold:
                    mut_class = "Surface"
                else:
                    mut_class = "Core"
            else:
                mut_class = "Contact %s" % chain_type
        elif chain_sub_dist == "NONE":
            if float(lig_sub_dist) > 5.0:
                if rel_sur_acc == None:
                    mut_class = "Unknown"
                elif float(rel_sur_acc) > surface_threshold:
                    mut_class = "Surface"
                else:
                    mut_class = "Core"
            else:
                mut_class = "Contact Ligand"
                
                
        elif (float(lig_sub_dist) > 5.0 and float(chain_sub_dist) > 5.0):
            if rel_sur_acc == None:
                mut_class = "Unknown"
            elif float(rel_sur_acc) > surface_threshold:
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
        print "Error in class definition"
        print lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc
        print sys.exc_info()
    if mut_class == "":
        print lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc
        
    return mut_class

def majority_vote(secs):
    class_dict = {'helix':0.0,'sheet':0.0,'disorder':0.0}
    sec_class = {'H':'helix','B':'disorder','E':'sheet','G':'helix','I':'helix','T':'disorder','S':'disorder'}
    for (sec,qual) in secs:
        if sec in sec_class:
            class_dict[sec_class[sec]] += qual

    if class_dict['disorder'] >= class_dict['helix']:
        if class_dict['disorder'] >= class_dict['sheet']:
            return 'disorder'
        else:
            return 'sheet'
    elif class_dict['helix'] >= class_dict['sheet']:
        return 'helix'
    else:
        return 'sheet'   

def getChemicalDistance(aac):
    #print aac
    try:
        if aac.count(',') < 1:
            chemical_distance = chem_dist_matrix[aac[0]][aac[-1]]
        else:
            aa1 = aac[0]
            aa2s = aac.split(",")
            aa2s[0] = aa2s[0][-1]
            chem_dists = []
            for aa2 in aa2s:
                chem_dists.append(chem_dist_matrix[aa1][aa2])
            chemical_distance = float(sum(chem_dists))/float(len(chem_dists))
    except:
        return None
    return chemical_distance

def getBlosumValue(aac):
    if aac.count(',') < 1:
        try:
            try:
                blosum_value = MatrixInfo.blosum62[(aac[0],aac[-1])]
            except:
                blosum_value = MatrixInfo.blosum62[(aac[-1],aac[0])]
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
                    blosum_value = MatrixInfo.blosum62[(aa1,aa2)]
                except:
                    blosum_value = MatrixInfo.blosum62[(aa2,aa1)]
            except:
                blosum_value = 0.0
            bvs.append(blosum_value)
        blosum_value = float(sum(bvs))/float(len(bvs))
    return blosum_value

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
        [e,f,g] = sys.exc_info()
        raise NameError("Error in getLigandAnnotation: %s\n%s" % (sql,f))
    for row in results:
        ligand_name = row[1]
        if not ligand_name in chosen_ones:
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

    print "TESTPRINT: ",len(template_mutation_map)

    min_t_id = min(template_ids)
    max_t_id = max(template_ids)

    t = 50000
    old_max = None
    if max_t_id - min_t_id > t:
        old_max = max_t_id
        max_t_id = min_t_id + t

    while True:
        sql = "SELECT Mutation,Template,Sub_Lig_Dist FROM RS_Mutation_Template WHERE Template BETWEEN %s AND %s" % (str(min_t_id),str(max_t_id))
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

        if old_max == None:
            break
        else:
            min_t_id = max_t_id + 1
            max_t_id = old_max
            if max_t_id - min_t_id > t:
                old_max = max_t_id
                max_t_id = min_t_id + t
            else:
                old_max = None

    #print "TESTPRINT: ",len(template_mutation_map)

    for lig in chosen_ones:
        if not lig in ligand_map:
            continue
        lig_id = ligand_map[lig]
        filtered_annotation_ids = {}
        if not lig_id in ligand_template_map:
            continue
        for template_id in ligand_template_map[lig_id]:
            if not template_id in template_mutation_map:
                continue
            for mutation_id in template_mutation_map[template_id]:
                
                sub_lig_distances = template_mutation_map[template_id][mutation_id]
                if sub_lig_distances != None:
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
                        good_dists = sorted(good_dists,key=lambda x:float(x.rsplit(":")[1]))
                        if template_id in filtered_annotation_ids:
                            filtered_annotation_ids[template_id].append((mutation_id,good_dists))
                        else:
                            filtered_annotation_ids[template_id] = [(mutation_id,good_dists)]
        if len(filtered_annotation_ids) > 0:
            for template_id in filtered_annotation_ids:
                filtered_annotation_ids[template_id] = sorted(filtered_annotation_ids[template_id],key=lambda x:float(x[1][0].split(":")[1]))
            anno_dict[lig] =  (chosen_ones[lig],filtered_annotation_ids)
    return anno_dict

#called by output
def minDistOut(outfolder,session_name,session_id,db,cursor,ligand_filter=None,intertable=False,templ_qual_separation=False,overwrite=False):
    outfile = '%s/%s' % (outfolder,session_name)

    if os.path.isfile("%s.classification.tsv" % outfile):
        print "A classification file already exists:\n%s.classification.tsv\n" % outfile
        if overwrite:
            print "Overwrite flag is active, classification gets recalculated and old files will be overwritten."
        else:
            print "The classification gets skipped, if you wish to reclassify, use the --overwrite flag."
            files = os.listdir(outfolder)
            class_files = []
            for fname in files:
                if fname.count('.classification') > 0:
                    class_files.append('%s/%s' % (outfolder,fname))
            return class_files

    if ligand_filter != None:
        f = open(ligand_filter,'r')
        lines = f.readlines()
        f.close()
        ligand_filter = set()
        for line in lines:
            ligand_filter.add(line.replace('\n','').replace(' ',''))
    global chem_dist_matrix

    sql = "SELECT Mutation,New_AA,Tag FROM RS_Mutation_Session WHERE Session = '%s'" % str(session_id)
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in minDistOut2: %s,\n%s" % (sql,f))
    
    new_aa_map = {}
    tag_map = {}
    for row in results:
        new_aa_map[row[0]] = row[1]
        tag_map[row[0]] = row[2]

    sql = "SELECT Mutation FROM RS_Mutation_Session WHERE Session = '%s'" % str(session_id)
    
    cursor = db.cursor()
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in minDistOut: %s,\n%s" % (sql,f))

    if results == ():
        return []




    mutation_id_list = set()
    max_mut_id = 0
    min_mut_id = None

    binsize = 10000
    bins = set([])

    for row in results:
        mut_id = row[0]

        bin_number = mut_id//binsize
        bins.add(bin_number)

        if mut_id > max_mut_id:
            max_mut_id = mut_id
        if min_mut_id == None or mut_id < min_mut_id:
            min_mut_id = mut_id
        mutation_id_list.add(mut_id)

    t2 = time.time()

    min_max_tuples = []
    if max_mut_id - min_mut_id > binsize:
        for bin_number in bins:
            min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize))
    else:
        min_max_tuples = [(min_mut_id,max_mut_id)]

    total_results = []
    for (min_mut_id,max_mut_id) in min_max_tuples:

        sql = ("SELECT Mutation,Sub_Lig_Dist,Sub_Chain_Distances,Template,Relative_Surface_Access,Secondary_Structure_Assignment,Ligand_Interaction_Degree,Ligand_Interaction_Score,Chain_Interaction_Degree,Chain_Interaction_Score,Short_Interaction_Degree,Short_Interaction_Score,Medium_Interaction_Degree,	Medium_Interaction_Score,Long_Interaction_Degree,Long_Interaction_Score FROM RS_Mutation_Template WHERE Error IS NULL AND Mutation BETWEEN %s AND  %s") % (str(min_mut_id),str(max_mut_id))
        try:
            # Execute the SQL command
            cursor.execute("SELECT Mutation,Sub_Lig_Dist,Sub_Chain_Distances,Template,Relative_Surface_Access,Secondary_Structure_Assignment,Ligand_Interaction_Degree,Ligand_Interaction_Score,Chain_Interaction_Degree,Chain_Interaction_Score,Short_Interaction_Degree,Short_Interaction_Score,Medium_Interaction_Degree,	Medium_Interaction_Score,Long_Interaction_Degree,Long_Interaction_Score FROM RS_Mutation_Template WHERE Error IS NULL AND Mutation BETWEEN %s AND  %s", (min_mut_id, max_mut_id))
            anno_results = cursor.fetchall()
            # Commit your changes in the database
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in minDistOut: %s,\n%s" % (sql,f))

        for row in anno_results:
            mut_id = row[0]
            if mut_id in mutation_id_list:
                total_results.append(row)


    template_id_list = set()
    for row in total_results:
        template_id_list.add(row[3]) 
    #print template_id_list
    template_dict = {}
    sql = "SELECT Template_Id,Chains,Sequence_Identity,Name,Original_Chains,Quality_Score,Alignment_Length FROM Template"
    try:
        cursor.execute(sql)
        temp_results = cursor.fetchall()
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in minDistOut: %s,\n%s" % (sql,f))
    
    seq_ids = set()
    ba_ids = set()
    seq_ba_ids = set()
    template_map = {}
    for row in temp_results:
        t_id = row[0]
        if t_id in template_id_list:
            chains = str(row[1])
            seq_id = float(row[2])
            pdb_id = row[3]
            ori_chains = row[4]
            qual = float(row[5])
            cov = float(row[6])
            cis = chains.split(";")
            template_dict[t_id] = {}
            for ci in cis:
                template_dict[t_id][ci[0]] = ci[2:]
            if seq_id >= 90.0:
                seq_ids.add(t_id)
            if not pdb_id.count('_na') > 0:
                ba_ids.add(t_id)
                if seq_id >= 90.0:
                    seq_ba_ids.add(t_id)
            template_map[t_id] = (pdb_id,ori_chains,qual,cov)

    temp_results_seq = []
    temp_results_ba = []
    temp_results_seq_ba = []
    for row in total_results:
        t_id = row[3]
        if t_id in seq_ids:
            temp_results_seq.append(row)
        if t_id in ba_ids:
            temp_results_ba.append(row)
        if t_id in seq_ba_ids:
            temp_results_seq_ba.append(row)

    mutation_dict = getMutationDict(mutation_id_list,db,cursor)
    gene_id_list = set()
    for m in mutation_dict:
        gene_id_list.add(mutation_dict[m][1])
    sql = "SELECT Gene_Id,Uniprot_Ac,Uniprot_Id,Species FROM Gene"
    try:
        cursor.execute(sql)
        gene_results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("NameError in minDistOut: %s" % sql)
    g_u_dict = {}
    for row in gene_results:
        gene_id = row[0]
        if gene_id in gene_id_list:
            g_u_dict[gene_id] = (row[1],row[2],row[3])

    if templ_qual_separation:
        options_map = {"":total_results,".bioassembly":temp_results_ba,".seqid90":temp_results_seq,".seqid90.bioassembly":temp_results_seq_ba}
    else:
        options_map = {"":total_results}

    class_files = []
    for option in options_map:
        temp_results = options_map[option]
        min_l_dict,min_m_dict,min_c_dict,min_r_dict,min_d_dict,mutation_surface_dict,mutation_sec_dict,mutation_inter_dict = createTemplateDicts(temp_results,template_dict,template_map,ligand_filter)
        class_dict = createClassDict(min_l_dict,min_m_dict,min_c_dict,min_r_dict,min_d_dict,mutation_surface_dict,mutation_sec_dict)
        if intertable:
            inter_dict = createInterDict(mutation_inter_dict)

        class_file = "%s.classification%s.tsv" % (outfile,option)
        class_files.append(class_file)
        writeClassFile(class_file,mutation_surface_dict,mutation_sec_dict,mutation_dict,g_u_dict,class_dict,template_map,new_aa_map,tag_map)
        if intertable:
            writeInterFile("%s.interaction_profiles%s.tsv" % (outfile,option),inter_dict,mutation_dict,g_u_dict,template_map,new_aa_map,tag_map)

    return class_files

#method for comparing/developing the new classifiction
def diffSurfs(mutation_surface_dict,g_u_dict,mutation_dict,outfile):
    surface_threshold = 0.16
    class_dict = {}
    for m in mutation_surface_dict:
        (u_ac,u_id) = g_u_dict[mutation_dict[m][1]]
        aac = mutation_dict[m][0]
        DSC = 0.0 #decision sum core
        DSS = 0.0 #decision sum surface
        n = 0.0
        min_surf = 2.0
        lines = ["Coverage\trASA"]
        for (surface,qual,cov) in mutation_surface_dict[m]:
            if surface < surface_threshold:
                DSC += qual*(cov**5)
            else:
                DSS += qual*(cov**5)
            lines.append("%s\t%s" % (str(cov),str(surface)))
            n += 1.0

            if surface < min_surf:
                min_surf = surface
                
        weighted_surface_value = DSS-2*DSC
        if weighted_surface_value > 0:
            weighted_sc = "Surface"
            if min_surf < surface_threshold:
                f = open("%s.%s_%s.tsv" % (outfile,u_ac,aac),'w')
                f.write("\n".join(lines))
                f.close()
        else:
            weighted_sc = "Core"

        conf_sc = (1.0-1.0/(n+1.0))*abs(DSS-2*DSC)/(DSS+2*DSC)

def createTemplateDicts(template_result_list,template_dict,template_map,ligand_filter=None):
    ions = set(["LI","NA","K","MG","CA","RB","CS","BE","SR","BA","SC","TI","V","CR","MN","FE","CO","NI","CU","ZN","F","CL","SI","AL","Y",
        "ZR","NB","MO","TC","RU","PD","AG","CD","IN","SN","GA","AS","SB","TE","I","BR","AT","PT","AU","HG","TL","PB"])

    min_l_dict = {}
    min_m_dict = {}
    min_c_dict = {}
    min_r_dict = {}
    min_d_dict = {}

    mutation_surface_dict = {}

    mutation_sec_dict = {}

    mutation_inter_dict = {}

    for row in template_result_list:
        #print row
        m_id = row[0]
        sld = str(row[1])
        scd = str(row[2])
        t_id = row[3]
        (pdb_id,ori_chains,qual,cov) = template_map[t_id]
        try:
            rel_sur_acc = float(row[4])
            if rel_sur_acc > 1.0:
                rel_sur_acc = 0.01*rel_sur_acc
        except:
            rel_sur_acc = None
        try:
            sec_str_ass = str(row[5])
        except:
            sec_str_ass = None

        #print row

        Ligand_Interaction_Degree = row[6]
        Ligand_Interaction_Score = row[7]
        Chain_Interaction_Degree = row[8]
        Chain_Interaction_Score = row[9]
        Short_Interaction_Degree = row[10]
        Short_Interaction_Score = row[11]
        Medium_Interaction_Degree = row[12]
        Medium_Interaction_Score = row[13]
        Long_Interaction_Degree = row[14]
        Long_Interaction_Score = row[15]

        if not m_id in mutation_inter_dict:
            mutation_inter_dict[m_id] = [(Ligand_Interaction_Degree,
                                    Ligand_Interaction_Score,
                                    Chain_Interaction_Degree,
                                    Chain_Interaction_Score,
                                    Short_Interaction_Degree,
                                    Short_Interaction_Score,
                                    Medium_Interaction_Degree,
                                    Medium_Interaction_Score,
                                    Long_Interaction_Degree,
                                    Long_Interaction_Score,
                                    qual)]
        else:
            mutation_inter_dict[m_id].append((Ligand_Interaction_Degree,
                                    Ligand_Interaction_Score,
                                    Chain_Interaction_Degree,
                                    Chain_Interaction_Score,
                                    Short_Interaction_Degree,
                                    Short_Interaction_Score,
                                    Medium_Interaction_Degree,
                                    Medium_Interaction_Score,
                                    Long_Interaction_Degree,
                                    Long_Interaction_Score,
                                    qual))

        min_ld = None
        min_md = None
        min_cd = None
        min_dd = None
        min_rd = None

        slds = sld.split(",")
        ldists = []
        mdists = []
        for ld in slds:
            ldi = ld.split(":")
            if len(ldi) > 1:
                lig_name = ldi[0].split('_')[0]
                if ligand_filter == None:
                    if lig_name in ions:
                        mdists.append(float(ldi[1]))
                    else:
                        ldists.append(float(ldi[1]))
                elif not lig_name in ligand_filter:
                    if lig_name in ions:
                        mdists.append(float(ldi[1]))
                    else:
                        ldists.append(float(ldi[1]))
        if len(ldists) > 0:
            min_ld = min(ldists)

        if len(mdists) > 0:
            min_md = min(mdists)

        scds = scd.split(",")
        cdists = []
        ddists = []
        rdists = []        
        for sd in scds:
            sdi = sd.split(":")
            if len(sdi) > 1:
                chain = sdi[0][0]
                mc_d = float(sdi[1])
                try:
                    chaintype = template_dict[t_id][chain]     
                    if chaintype == "Protein":
                        cdists.append(mc_d)
                    elif chaintype == "RNA":
                        rdists.append(mc_d)
                    elif chaintype == "DNA":
                        ddists.append(mc_d)
                except:
                    print "Strange Error: ",t_id,m_id,chain

        if len(cdists) > 0:
            min_cd = min(cdists)
        if len(rdists) > 0:
            min_rd = min(rdists)
        if len(ddists) > 0:
            min_dd = min(ddists)

        minimal_distances = []
        if min_cd != None:
            minimal_distances.append(min_cd)
        if min_dd != None:
            minimal_distances.append(min_dd)
        if min_rd != None:
            minimal_distances.append(min_rd)
        if min_ld != None:
            minimal_distances.append(min_ld)
        if min_md != None:
            minimal_distances.append(min_md)

        if len(minimal_distances) == 0:
            min_minimal_distances = 2.0
        else:
            min_minimal_distances = min(minimal_distances)

        if not min_minimal_distances < 1.2:
            if rel_sur_acc != None and cov > 0.0:
                if m_id not in mutation_surface_dict:
                    mutation_surface_dict[m_id] = [(rel_sur_acc,qual,cov)]
                else:
                    mutation_surface_dict[m_id].append((rel_sur_acc,qual,cov))

            if not m_id in mutation_sec_dict:
                mutation_sec_dict[m_id] = [(sec_str_ass,qual)]
            else:
                mutation_sec_dict[m_id].append((sec_str_ass,qual))

            if not min_ld == None:
                if not m_id in min_l_dict:
                    min_l_dict[m_id] = [(min_ld,qual)]
                else:
                    min_l_dict[m_id].append((min_ld,qual))
            if not min_md == None:
                if not m_id in min_m_dict:
                    min_m_dict[m_id] = [(min_md,qual)]
                else:
                    min_m_dict[m_id].append((min_md,qual))
            if not min_cd == None:
                if not m_id in min_c_dict:
                    min_c_dict[m_id] = [(min_cd,qual)]
                else:
                    min_c_dict[m_id].append((min_cd,qual))
            if not min_rd == None:
                if not m_id in min_r_dict:
                    min_r_dict[m_id] = [(min_rd,qual)]
                else:
                    min_r_dict[m_id].append((min_rd,qual))
            if not min_dd == None:
                if not m_id in min_d_dict:
                    min_d_dict[m_id] = [(min_dd,qual)]
                else:
                    min_d_dict[m_id].append((min_dd,qual))

    return min_l_dict,min_m_dict,min_c_dict,min_r_dict,min_d_dict,mutation_surface_dict,mutation_sec_dict,mutation_inter_dict

def createInterDict(mutation_inter_dict):

    inter_dict = {}

    for m_id in mutation_inter_dict:
        profiles = mutation_inter_dict[m_id]
        total_qual = 0.0
        average_profile = [0.0]*10
        for profile in profiles:
            if None in profile:
                continue
            (Ligand_Interaction_Degree,
            Ligand_Interaction_Score,
            Chain_Interaction_Degree,
            Chain_Interaction_Score,
            Short_Interaction_Degree,
            Short_Interaction_Score,
            Medium_Interaction_Degree,
            Medium_Interaction_Score,
            Long_Interaction_Degree,
            Long_Interaction_Score,
            qual) = profile
            total_qual += qual

            average_profile[0] += Ligand_Interaction_Degree*qual
            average_profile[1] += Ligand_Interaction_Score*qual
            average_profile[2] += Chain_Interaction_Degree*qual
            average_profile[3] += Chain_Interaction_Score*qual
            average_profile[4] += Short_Interaction_Degree*qual
            average_profile[5] += Short_Interaction_Score*qual
            average_profile[6] += Medium_Interaction_Degree*qual
            average_profile[7] += Medium_Interaction_Score*qual
            average_profile[8] += Long_Interaction_Degree*qual
            average_profile[9] += Long_Interaction_Score*qual
        if total_qual > 0.0:
            average_profile = [str(x/total_qual) for x in average_profile]
        inter_dict[m_id] = average_profile

    return inter_dict

def median(l):
    #print l
    n = len(l)
    l = sorted(l)
    if n % 2 == 0:
        med = (l[(n/2)-1]+l[n/2])/2.0
    else:
        med = l[(n-1)/2]
    return med

def getWeightedClass(weighted_sc,conf_sc,weighted_c,conf_c,weighted_d,conf_d,weighted_r,conf_r,weighted_l,conf_l,weighted_m,conf_m):
    dt = 5.0
    edt = 8.0
    #Decide on contact macromolecule
    if weighted_c < dt and weighted_d < dt and weighted_r < dt: #If there is an Protein,DNA and RNA present (artificial scenario)
        if conf_c >= conf_d and conf_c >= conf_r:
            contact_macro = "Protein"
            conf_macro = conf_c
        elif conf_d >= conf_r:
            contact_macro = "DNA"
            conf_macro = conf_d
        else:
            contact_macro = "RNA"
            conf_macro = conf_r
    elif weighted_c < dt and weighted_d < dt:
        if conf_c >= conf_d:
            contact_macro = "Protein"
            conf_macro = conf_c
        else:
            contact_macro = "DNA"
            conf_macro = conf_d
    elif weighted_c < dt and weighted_r < dt:
        if conf_c >= conf_r:
            contact_macro = "Protein"
            conf_macro = conf_c
        else:
            contact_macro = "RNA"
            conf_macro = conf_r
    elif weighted_d < dt and weighted_r < dt:
        if conf_d >= conf_r:
            contact_macro = "DNA"
            conf_macro = conf_d
        else:
            contact_macro = "RNA"
            conf_macro = conf_r
    elif weighted_c < dt:
        contact_macro = "Protein"
        conf_macro = conf_c
    elif weighted_d < dt:
        contact_macro = "DNA"
        conf_macro = conf_d
    elif weighted_r < dt:
        contact_macro = "RNA"
        conf_macro = conf_r

    elif weighted_c < edt and weighted_d < edt and weighted_r < edt:
        if conf_c >= conf_d and conf_c >= conf_r:
            contact_macro = "Protein far"
            conf_macro = conf_c
        elif conf_d >= conf_r:
            contact_macro = "DNA far"
            conf_macro = conf_d
        else:
            contact_macro = "RNA far"
            conf_macro = conf_r
    elif weighted_c < edt and weighted_d < edt:
        if conf_c >= conf_d:
            contact_macro = "Protein far"
            conf_macro = conf_c
        else:
            contact_macro = "DNA far"
            conf_macro = conf_d
    elif weighted_c < edt and weighted_r < edt:
        if conf_c >= conf_r:
            contact_macro = "Protein far"
            conf_macro = conf_c
        else:
            contact_macro = "RNA far"
            conf_macro = conf_r
    elif weighted_d < edt and weighted_r < edt:
        if conf_d >= conf_r:
            contact_macro = "DNA far"
            conf_macro = conf_d
        else:
            contact_macro = "RNA far"
            conf_macro = conf_r
    elif weighted_c < edt:
        contact_macro = "Protein far"
        conf_macro = conf_c
    elif weighted_d < edt:
        contact_macro = "DNA far"
        conf_macro = conf_d
    elif weighted_r < edt:
        contact_macro = "RNA far"
        conf_macro = conf_r

    else:
        contact_macro = None
        conf_macro = None

    #Decide on contact low weight molecule
    if weighted_l < dt and weighted_m < dt:
        if conf_l >= conf_m:
            contact_ligand = "Ligand"
            conf_ligand = conf_l
        else:
            contact_ligand = "Metal"
            conf_ligand = conf_m
    elif weighted_l < dt:
        contact_ligand = "Ligand"
        conf_ligand = conf_l
    elif weighted_m < dt:
        contact_ligand = "Metal"
        conf_ligand = conf_m

    elif weighted_l < edt and weighted_m < edt:
        if conf_l >= conf_m:
            contact_ligand = "Ligand far"
            conf_ligand = conf_l
        else:
            contact_ligand = "Metal far"
            conf_ligand = conf_m
    elif weighted_l < edt:
        contact_ligand = "Ligand far"
        conf_ligand = conf_l
    elif weighted_m < edt:
        contact_ligand = "Metal far"
        conf_ligand = conf_m

    else:
        contact_ligand = None
        conf_ligand = None

    #Classification Decision Tree
    if contact_macro != None and contact_ligand != None:
        clas = "Double Interaction: %s and %s" % (contact_macro,contact_ligand)
        conf = (conf_ligand+conf_macro)/2.0
    elif contact_macro != None:
        clas = "%s Interaction" % contact_macro
        conf = conf_macro
    elif contact_ligand != None:
        clas = "%s Interaction" % contact_ligand
        conf = conf_ligand
    else:
        clas = weighted_sc
        conf = conf_sc

    return clas,conf

def createClassDict(min_l_dict,min_m_dict,min_c_dict,min_r_dict,min_d_dict,mutation_surface_dict,mutation_sec_dict):
    distance_threshold = 8.0
    surface_threshold = 0.16
    class_dict = {}
    for m in mutation_surface_dict:
        DSC = 0.0 #decision sum core
        DSS = 0.0 #decision sum surface
        n = 0.0
        qs = []
        for (surface,qual,cov) in mutation_surface_dict[m]:
            if surface < surface_threshold:
                DSC += qual*(cov**5)
            else:
                DSS += qual*(cov**5)
            qs.append(qual)
            n += 1.0

        weighted_surface_value = DSS-2*DSC
        if weighted_surface_value > 0:
            weighted_sc = "Surface"
        else:
            weighted_sc = "Core"
        conf_sc = (1.0-1.0/(n+1.0))*(abs(weighted_surface_value)/(DSS+2*DSC))*median(qs)

        if m in min_c_dict:
            nom = 0.0
            denom = 0.0
            n = 0.0
            qs = []
            for (min_c,qual) in min_c_dict[m]:
                if min_c < distance_threshold:
                    nom += min_c*qual
                    denom += qual
                    n += 1.0
                    qs.append(qual)
            if denom > 0:
                weighted_c = nom/denom
                conf_c = (1.0-1.0/(n+1.0))*(max(qs)+median(qs))/2
            else:
                weighted_c = 10.0
                conf_c = 0.0
        else:
            weighted_c = 10.0
            conf_c = 0.0

        if m in min_d_dict:
            nom = 0.0
            denom = 0.0
            n = 0.0
            qs = []
            for (min_d,qual) in min_d_dict[m]:
                if min_d < distance_threshold:
                    nom += min_d*qual
                    denom += qual
                    n += 1.0
                    qs.append(qual)
            if denom > 0.0:
                weighted_d = nom/denom
                conf_d = (1.0-1.0/(n+1.0))*(max(qs)+median(qs))/2
            else:
                weighted_d = 10.0
                conf_d = 0.0
        else:
            weighted_d = 10.0
            conf_d = 0.0

        if m in min_r_dict:
            nom = 0.0
            denom = 0.0
            n = 0.0
            qs = []
            for (min_r,qual) in min_r_dict[m]:
                if min_r < distance_threshold:
                    nom += min_r*qual
                    denom += qual
                    n += 1.0
                    qs.append(qual)
            if denom > 0.0:
                weighted_r = nom/denom
                conf_r = (1.0-1.0/(n+1.0))*(max(qs)+median(qs))/2
            else:
                weighted_r = 10.0
                conf_r = 0.0
        else:
            weighted_r = 10.0
            conf_r = 0.0

        if m in min_l_dict:
            nom = 0.0
            denom = 0.0
            n = 0.0
            qs = []
            for (min_l,qual) in min_l_dict[m]:
                if min_l < distance_threshold:
                    nom += min_l*qual
                    denom += qual
                    n += 1.0
                    qs.append(qual)
            if denom > 0.0:
                weighted_l = nom/denom
                conf_l = (1.0-1.0/(n+1.0))*(max(qs)+median(qs))/2
            else:
                weighted_l = 10.0
                conf_l = 0.0
        else:
            weighted_l = 10.0
            conf_l = 0.0

        if m in min_m_dict:
            nom = 0.0
            denom = 0.0
            n = 0.0
            qs = []
            for (min_m,qual) in min_m_dict[m]:
                if min_m < distance_threshold:
                    #print min_m,qual
                    nom += min_m*qual
                    denom += qual
                    n += 1.0
                    qs.append(qual)
            if denom > 0.0:
                weighted_m = nom/denom
                conf_m = (1.0-1.0/(n+1.0))*(max(qs)+median(qs))/2
            else:
                weighted_m = 10.0
                conf_m = 0.0
        else:
            weighted_m = 10.0
            conf_m = 0.0
        
        mut_class,conf = getWeightedClass(weighted_sc,conf_sc,weighted_c,conf_c,weighted_d,conf_d,weighted_r,conf_r,weighted_l,conf_l,weighted_m,conf_m)
        class_dict[m] = mut_class,conf,weighted_sc
    return class_dict

def simplifyClass(c):
    if c == "Surface" or c == "Core":
        return c
    if c.count("Double") > 0:
        return "Ligand Interaction"
    if c.count("Ligand") > 0 or c.count("Metal") > 0:
        return "Ligand Interaction"
    if c.count("DNA") > 0:
        return "DNA Interaction"
    if c.count("RNA") > 0:
        return "RNA Interaction"
    if c.count("Protein") > 0:
        return "Protein Interaction"
    else:
        print c

def writeInterFile(outfile,inter_dict,mutation_dict,g_u_dict,template_map,new_aa_map,tag_map):
    startline = "Uniprot\tAAC\tSpecie\tTag\tLigand_Interaction_Degree\tLigand_Interaction_Score\tChain_Interaction_Degree\tChain_Interaction_Score\tShort_Interaction_Degree\tShort_Interaction_Score\tMedium_Interaction_Degree\tMedium_Interaction_Score\tLong_Interaction_Degree\tLong_Interaction_Score"
    lines = [startline]
    for m in inter_dict:
        aac = mutation_dict[m][0]

        new_aa = new_aa_map[m]
        aac = "%s%s" % (aac.split(',')[0][:-1],new_aa)
        (u_ac,u_id,species) = g_u_dict[mutation_dict[m][1]]

        if m in inter_dict:
            interstr = '\t'.join([str(x) for x in inter_dict[m]])
        else:
            interstr = '\t'.join((['None']*10))

        line = "%s\t%s\t%s\t%s\t%s" % (u_ac,aac,species,tag_map[m],interstr)
        lines.append(line)

    f = open(outfile,'w')
    f.write("\n".join(lines))
    f.close()

def writeClassFile(outfile,mutation_surface_dict,mutation_sec_dict,mutation_dict,g_u_dict,class_dict,template_map,new_aa_map,tag_map):
    startline = "Uniprot\tAAC\tSpecies\tTag\tWeighted Surface/Core\tClass\tSimple Class\tConfidence Value\tSecondary Structure\tChemical Distance\tBlosum62 Value\tUniprot Id"

    lines = [startline]
    for m in mutation_surface_dict:
        if m in mutation_sec_dict:
            mv_sec_ass = majority_vote(mutation_sec_dict[m])
        else:
            mv_sec_ass = None

        Class,conf,weighted_sc = class_dict[m]
        simple_class = simplifyClass(Class)

        aac = mutation_dict[m][0]

        new_aa = new_aa_map[m]
        aac = "%s%s" % (aac.split(',')[0][:-1],new_aa)

        chemical_distance = getChemicalDistance(aac)
        blosum_value = getBlosumValue(aac)
        (u_ac,u_id,species) = g_u_dict[mutation_dict[m][1]]
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (u_ac,mutation_dict[m][0],species,tag_map[m],weighted_sc,Class,simple_class,str(conf),mv_sec_ass,str(chemical_distance),str(blosum_value),str(u_id)))
    f = open(outfile,'w')
    f.write("\n".join(lines))
    f.close()



#called by output
def prodAnoOut(output_file,session_id,db_name,db_adress,db_user_name,db_password,top=None,ligand_filter=None,proteome=False):
    if ligand_filter != None:
        f = open(ligand_filter,'r')
        lines = f.readlines()
        f.close()
        ligand_filter = set()
        for line in lines:
            ligand_filter.add(line.replace('\n','').replace(' ',''))    

    t0 = time.time()
    db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    cursor = db.cursor()

    t03 = time.time()

    mutation_id_list = set()
    template_id_list = set()

    print db_name,session_id

    t035 = time.time()
    sql = "SELECT Mutation,Template FROM RS_Annotation_Session WHERE Session = '%s'" % (str(session_id))
    if proteome:
        sql = "SELECT Mutation,Template FROM RS_Mutation_Template"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in prodAnoOut: %s" % sql)
    t04 = time.time()
    sql_time = t04 - t035
    print("SQL time: " + str(sql_time))
    if not results == ():

        binsize = 10000
        bins = set([])
        max_mut_id = 0
        min_mut_id = None
        for row in results:
            mut_id = row[0]

            bin_number = mut_id//binsize
            bins.add(bin_number)

            if mut_id > max_mut_id:
                max_mut_id = mut_id
            if min_mut_id == None or mut_id < min_mut_id:
                min_mut_id = mut_id
            mutation_id_list.add(mut_id)
            template_id_list.add(row[1])
        t045 = time.time()
        p_time = t045 - t04
        print("Python set time: " + str(p_time))


    sql = "SELECT Mutation,New_AA,Tag FROM RS_Mutation_Session WHERE Session = '%s'" % (str(session_id))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in prodAnoOut: %s,\n%s" % (sql,f))
    
    new_aa_map = {}
    tag_map = {}
    for row in results:
        new_aa_map[row[0]] = row[1]
        tag_map[row[0]] = row[2]

    #print len(template_id_list)
    total_results = []
    if not len(template_id_list) == 0:

        t2 = time.time()

        min_max_tuples = []
        if max_mut_id - min_mut_id > binsize:
            for bin_number in bins:
                min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize))
        else:
            min_max_tuples = [(min_mut_id,max_mut_id)]

        
        for (min_mut_id,max_mut_id) in min_max_tuples:
            results = ()
            sql = "SELECT * FROM RS_Mutation_Template WHERE Error IS NULL AND Mutation BETWEEN %s AND %s" % (str(min_mut_id),str(max_mut_id))
            try:
                # Execute the SQL command
                cursor.execute(sql)
                results = cursor.fetchall()
                # Commit your changes in the database
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in prodAnoOut: %s,\n%s" % (sql,f))

            if not results == ():
                for row in results:
                    template_id = row[1]
                    if template_id in template_id_list:
                        mutation_id = row[0]
                        if mutation_id in mutation_id_list:
                            if len(row) < 5:
                                print row,min_mut_id,max_mut_id
                            total_results.append(row)

        if len(total_results) == 0:
            t06 = time.time()  
            t07 = time.time()
            template_dict = {}
            mutation_dict = {}
            reduced = False
            return (output_file,reduced)

        t06 = time.time()  
          
        if top != None:
            len_before = len(total_results)
            total_results = sorted(total_results,key = lambda x: float(x[7]),reverse=True)
            total_results = total_results[:top]
            len_after = len(total_results)
            mutation_id_list = []
            template_id_list = []
            if len_before != len_after:
                reduced = True
            else:
                reduced = False
        else:
            reduced = False

        mutation_id_list = set()
        template_id_list = set()
        for row in total_results:
            mutation_id_list.add(row[0])
            template_id_list.add(row[1])

        mutation_dict = getMutationDict(mutation_id_list,db,cursor)            
        template_dict = getTemplatesFromIdList(template_id_list,db,cursor)

        t07 = time.time()
            
            
    else:
        t05 = time.time()
        t06 = time.time()  
        t07 = time.time()
        template_dict = {}
        mutation_dict = {}
        reduced = False

    t08 = time.time()

    #print "template_id_lst length: ",len(template_id_list)

    gene_id_list = set([x[1] for x in template_dict.values()])
    #print gene_id_list

    t09 = time.time()

    if len(gene_id_list) == 0:
        gene_score_dict = {}
    else:
        gene_score_dict = getGeneScoreDict(gene_id_list,session_id,db,cursor)
    cursor.close()
    db.close()


    t1 = time.time()
    
    lines = produceAnotationOutput(session_id,total_results,template_dict,gene_score_dict,mutation_dict,ligand_filter,new_aa_map,tag_map)

    startline = "\t".join(("PDB","Gene",'Species',"SP ID","REFSEQS","AA change","Residue_id(Template)","Sequence Id","Alignment Lenght","Resolution","R-Factor","Target Chain","Ligands(Name,Residue Nr)","Sub-Lig-Dists","Min Sub-Lig-Dist","Chain Type","Sub Chain Distances","Min Sub-Chain-Dist","Chaintypes","Relative Surface Access","Secondary Structure Assignment","Quality Score","Candidate Score","Combined Score","Gene_Score","PP2 Score","Class","Chemical Distance","Blosum62 Value","Tag"))
    print len(lines)
    lines = sorted(lines,key= lambda x: float(x[22]),reverse=True)
    lines2 = [startline]
    for line in lines:
        try:
            line2 = "\t".join(line)
            lines2.append(line2)
        except:
            print line
    page = "\n".join(lines2)
    f = open(output_file, "wb")
    f.write(page)
    f.close()

    print "%s done" % output_file
    t2 = time.time()
    totaltime = t2 - t0
    database_prep_time = t03 - t0

    rest2_prep_time = t07 - t06
    rest3_prep_time = t08 - t07
    rest4_prep_time = t09 - t08
    rest5_prep_time = t1 - t09
    para_time = t2 - t1
    print("Database Prep time: " + str(database_prep_time))

    print("Rest2 Prep time: " + str(rest2_prep_time))
    print("Rest3 Prep time: " + str(rest3_prep_time))
    print("Rest4 Prep time: " + str(rest4_prep_time))
    print("Rest5 Prep time: " + str(rest5_prep_time))
    print("Para time: " + str(para_time))
    print("Total time: " + str(totaltime))
    return (output_file,reduced)

def produceAnotationOutput(session_id,input_queue,template_dict,gene_score_dict,mutation_dict,ligand_filter,new_aa_map,tag_map):
    global chem_dist_matrix
    outlines = []
    for row in input_queue:
        if len(row) < 9:
            print "Error: too less words in row: ",row
            continue
        mutation_id = row[0]
        template_id = row[1]
        sub_lig_dist = row[2]
        sub_chain_dists = row[3]
        rel_sur_acc = row[4]
        sec_str_ass = row[5]
        residue_id = row[6]
        candidate = row[8]
        (template,gene_id,chains) = template_dict[template_id]
        if not gene_id in gene_score_dict:
            print "Possible error for: ",mutation_id,template_id,gene_id
            continue
        (gene_name,gpan,sp_id,error_code,error,species,gene_score) = gene_score_dict[gene_id]

        # template: [pdb_id,seq_id,target_chain,aln_length,resolution,ligands,r_value,score,sub_lig_dist,sub_ano,original_target_chain,original_chains]
        pdb_id = template[0]
        seq_id = template[1]
        aln_length = template[3]
        chain = template[2]
        resolution = template[4]
        ligands = template[5]
        r_value = template[6]
        quality = template[7]
        
        tag = tag_map[mutation_id]

        (aachange,gene_id,pp2) = mutation_dict[mutation_id]
        new_aa = new_aa_map[mutation_id]
        aachange = "%s%s" % (aachange.split(',')[0][:-1],new_aa)
        #if sp_id == 'RET_HUMAN':
        #    if mutation_id == 843:
        #        print aachange,mutation_id      

        ds = []
        if sub_lig_dist != None:
            lig_dists = sub_lig_dist.split(",")
            for lig_dist in lig_dists:
                ldp = lig_dist.split(":")
                if len(ldp) > 1:
                    lig_name = ldp[0].split('_')[0]
                    if ligand_filter == None:
                        dfl = list(ldp[1])
                        if len(dfl) > 1:
                            ds.append(float(ldp[1]))
                    elif not lig_name in ligand_filter:
                        dfl = list(ldp[1])
                        if len(dfl) > 1:
                            ds.append(float(ldp[1]))

        if len(ds) > 0:
            lig_sub_dist = min(ds)
        else:
            lig_sub_dist = "NONE"

        ds = []
        if sub_chain_dists != None:
            chain_dists = sub_chain_dists.split(",")
            for chain_dist in chain_dists:
                cdp = chain_dist.split(":")
                if len(cdp) > 1:
                    try:
                        ds.append([float(cdp[1]),cdp[0]])
                    except:
                        continue
        if len(ds) > 0:
            ds = sorted(ds,key = lambda x: x[0])
            chain_sub_dist = ds[0][0]
            chain_type = ds[0][1]
        else:
            chain_sub_dist = "NONE"
            chain_type = "NONE"

        if candidate != None:
            combined = candidate*quality
            candidate_str = "%.3f" % candidate
            combined_str = "%.3f" % combined
        else:
            combined = None
            candidate_str = '0.0'
            combined_str = '0.0'
        

        #What does this?
        ligand_texts = []        
        for ligand in ligands:
            if ligand[0] == "Ligand":
                t = "%s,%s" %(ligand[1],ligand[2])
                ligand_texts.append(t)
        ligtext = ";".join(ligand_texts)

        chain_type_dict = {}
        for c_pair in chains.split(";"):
            chain_type_dict[c_pair.split(":")[0]] = c_pair.split(":")[1]
        if chain_type == "NONE":
            class_chain_type = chain_type
        else:
            try:
                class_chain_type = chain_type_dict[chain_type[0]]
            except:
                chain_type = "NONE"
                class_chain_type = "NONE"

        mut_class = getClass(lig_sub_dist,chain_sub_dist,class_chain_type,rel_sur_acc)

        chemical_distance = getChemicalDistance(aachange)
        blosum_value = getBlosumValue(aachange)        

        line = (str(pdb_id),str(gene_name),species,str(sp_id),str(gpan),str(aachange),str(residue_id),str(seq_id),str(aln_length),str(resolution),str(r_value),str(chain),str(ligtext),str(sub_lig_dist),str(lig_sub_dist),str(chain_type),str(sub_chain_dists),str(chain_sub_dist),str(chains),str(rel_sur_acc),str(sec_str_ass),"%.3f" % quality,candidate_str,combined_str,str(gene_score),str(pp2),mut_class,str(chemical_distance),str(blosum_value),tag)
        outlines.append(line)

    return outlines

#called by output
def updateGeneScores(session_id,db,cursor):
    sql = "SELECT Gene,Gene_Score FROM RS_Gene_Session WHERE Session = '%s'" % session_id
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error: %s,\n%s" % (sql,f))
    
    gene_mut_map = {}    
    gene_ids = set()
    for row in results:
        if row != None:
            #No need for updating, if already in database
            return
        gene_ids.add(row[0])
        gene_mut_map[row[0]] = set()

    min_g = min(gene_ids)
    max_g = max(gene_ids)

    sql = "SELECT Mutation FROM RS_Mutation_Session WHERE Session = '%s'" % session_id
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error: %s,\n%s" % (sql,f))

    legal_muts = set()
    min_m = None
    max_m = None
    for row in results:
        m_id = row[0]
        legal_muts.add(m_id)
        if min_m == None or m_id < min_m:
            min_m = m_id
        if max_m == None or m_id > max_m:
            max_m = m_id

    if len(gene_ids) > 0:
        sql = "SELECT Mutation_Id,Gene FROM Mutation WHERE Gene BETWEEN %s AND %s" % (str(min_g),str(max_g))
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error: %s,\n%s" % (sql,f))

        for row in results:
            if not row[1] in gene_ids:
                continue
            if row[0] in legal_muts:
                gene_mut_map[row[1]].add(row[0])

    t0 = time.time()
    geneScores = computeGeneScores(gene_mut_map,min_m,max_m,legal_muts,db,cursor)
    t1 = time.time()
    value_strs = []

    for gene in gene_mut_map:
        if gene in geneScores:
            gs = geneScores[gene]
        else:
            gs = 0.0
        value_strs.append("WHEN '%s' THEN '%s'" % (gene,str(gs)))
    sql = "UPDATE IGNORE RS_Gene_Session SET Gene_Score = CASE GENE %s ELSE Gene_Score END WHERE Session = '%s' AND Gene BETWEEN %s AND %s" % (' '.join(value_strs),str(session_id),str(min_g),str(max_g))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Commit your changes in the database
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error: %s" % (f))
    
    t2 = time.time()

    print "Time for calculating Gene Scores: ",t1-t0
    print "Time for updating Gene Scores into Database: ",t2-t1

def computeGeneScores(gene_mut_map,min_m,max_m,mut_ids,db,cursor):
    t0 = time.time()

    t1 = time.time()
    sql = "SELECT Candidate_Score,Template,Mutation FROM RS_Mutation_Template WHERE Error IS NULL AND Mutation BETWEEN %s AND %s" % (str(min_m),str(max_m))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in computeGeneScores: %s" % sql)
    t2 = time.time()

    template_map = {}
    
    for row in results:
        m_id = row[2]
        if not m_id in mut_ids:
            continue
        t_id = row[1]
        cs = float(row[0])
        if not t_id in template_map:
            template_map[t_id] = cs
        else:
            template_map[t_id] += cs
    if len(template_map) == 0:
        return {}

    template_ids = set(template_map.keys())

    max_t = max(template_ids)
    min_t = min(template_ids)

    t3 = time.time()
    sql = "SELECT Template_Id,Gene,Quality_Score FROM Template WHERE Template_Id <= %s AND Template_Id >= %s" % (str(max_t),str(min_t))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in computeGeneScores: %s" % sql)
    t4 = time.time()

    geneScores = {}
    for row in results:
        t_id = row[0]
        if not t_id in template_ids:
            continue
        g_id = row[1]
        template_score = float(row[2])*template_map[t_id]
        if not g_id in geneScores:
            geneScores[g_id] = template_score
        elif template_score > geneScores[g_id]:
            geneScores[g_id] = template_score
    t5 = time.time()

    #print "computeGeneScores detailed times: ",t1-t0,t2-t1,t3-t2,t4-t3,t5-t4

    return geneScores        

#called by output
def sortGenes(session_id,outfile,db,cursor):

    sql = "SELECT Gene,Gene_Score FROM RS_Gene_Session WHERE Session = '%s'" % session_id

    try:
        # Execute the SQL command
        cursor.execute(sql)
        gene_scores = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("NameError sortGenes %s" % sql)

    gene_scores = sorted(gene_scores,key = lambda x: x[1],reverse = True)

    lines = ["sp_id\tgene_id\tScore"]
    for gene_score in gene_scores:
        sp_id = getUniprotFromId(str(gene_score[0]),db,cursor)
        line = "%s\t%s\t%s" % (sp_id,str(gene_score[0]),str(gene_score[1]))
        lines.append(line)
        #updateGeneScore(gene_score[0],session_id,gene_score[2],db,cursor)
    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()

#called by output, reScore
def updateGeneScore(gene_id,session_id,score,db,cursor):
    sql = "UPDATE IGNORE RS_Gene_Session SET Gene_Score = '%s' WHERE Session = '%s' AND Gene = '%s'" % (str(score),str(session_id),str(gene_id))
    try:
        # Execute the SQL command
        cursor.execute(sql)
        # Commit your changes in the database
        db.commit()
    except:
        raise NameError("Error in updateGeneScore: %s" % sql)
        # Rollback in case there is any NameError
        db.rollback()  

#called by output    
def goTermAnalysis(session_id,outfile,db,cursor):
    #print  session_id
    sql = "SELECT Gene,Gene_Score FROM RS_Gene_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)
    datasetsize = len(results)
    go_dict = {}
    gene_dict = {}
    gene_id_list = set([])
    for row in results:
        gene_score = row[1]
        if gene_score == None:
            gene_score = 0.0
        gene_dict[row[0]] = gene_score
        gene_id_list.add(row[0])

    sql = "SELECT GO_Term,Gene FROM RS_Gene_GO_Term"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in goTermAnalyis: %s" % sql)
    go_term_id_list = set([])
    for row in results:
        if not row[1] in gene_id_list:
            continue
        gene_score = gene_dict[row[1]]
        if not row[0] in go_dict:
            go_term_id_list.add(row[0])
            go_dict[row[0]] = ["","",gene_score,1.0,gene_score]
        else:
            i = go_dict[row[0]][3] + 1
            avg_score = go_dict[row[0]][2]
            go_dict[row[0]][2] = ((i-1)/i)*avg_score + (1/i)*gene_score
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

    go_list = sorted(go_dict.values(),key = lambda x: x[4],reverse = True)
    lines = ["GO-Term\tGO-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
    for go in go_list:
        normalized_score = go[4]/datasetsize
        lines.append("%s\t%s\t%s\t%s\t%s\t%s" % (str(go[0]),str(go[1]),str(go[4]),str(go[2]),str(int(go[3])),str(normalized_score)))
    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()

#called by output
def pathwayAnalysis(session_id,outfile,db,cursor):
    sql = "SELECT Gene,Gene_Score FROM RS_Gene_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)
    datasetsize = len(results)
    path_dict = {}
    gene_dict = {}
    gene_id_list = set([])
    for row in results:
        gene_score = row[1]
        if gene_score == None:
            gene_score = 0.0
        gene_dict[row[0]] = gene_score
        gene_id_list.add(row[0])

    sql = "SELECT Pathway,Gene FROM RS_Gene_Pathway"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in pathwayAnalyis: %s" % sql)
    pathway_id_list = set([])
    for row in results:
        if not row[1] in gene_id_list:
            continue
        gene_score = gene_dict[row[1]]
        if not row[0] in path_dict:
            pathway_id_list.add(row[0])
            path_dict[row[0]] = ["","",gene_score,1.0,gene_score]
        else:
            i = path_dict[row[0]][3] + 1
            avg_score = path_dict[row[0]][2]
            path_dict[row[0]][2] = ((i-1)/i)*avg_score + (1/i)*gene_score
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

    path_list = sorted(path_dict.values(),key = lambda x: x[4],reverse = True)
    lines = ["Pathway\tReactome-ID\tTotal Score\tAVG-Score\tGene_Amount\tNormalized Score"]
    for path in path_list:
        normalized_score = path[4]/datasetsize
        lines.append("%s\t%s\t%s\t%s\t%s\t%s" % (str(path[0]),str(path[1]),str(path[4]),str(path[2]),str(int(path[3])),str(normalized_score)))
    page = "\n".join(lines)
    f = open(outfile, "wb")
    f.write(page)
    f.close()

#called by postAnnoAnno
def updateAnnoAnno(anno_anno_map,session_id,db,cursor):
    valuestrs = []
    for (m_id1,m_id2) in anno_anno_map:
        (min_d,atom,atom2,t_id1,t_id2,chain1,chain2) = anno_anno_map[(m_id1,m_id2)]
        valuestrs.append("('%s','%s','%s','%s','%s','%s','%s','%s','%s %s')" % (str(t_id1),str(t_id2),str(m_id1),str(m_id2),chain1,chain2,str(session_id),str(min_d),str(atom),str(atom2)))
    sql = "INSERT IGNORE INTO RS_Annotation_Annotation(Template_1,Template_2,Mutation_1,Mutation_2,Chain_1,Chain_2,Session,Distance,Atompair) VALUES %s" % ','.join(valuestrs)
    try:
        cursor.execute(sql)
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in updateAnnoAnno: %s" % (f))

#called by calcAnnoRate
def calculateAnnotationRate(db,cursor,session_id):
    sql = "SELECT Mutation FROM RS_Mutation_Session WHERE Session = %s" % session_id
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
        sql = "SELECT Mutation_Id,Gene FROM Mutation WHERE Mutation_Id in (%s)" % ",".join(mut_id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        gene_map = {}
        gene_id_list = set()
        for row in results:
            gene_id = str(row[1])
            mut_id = str(row[0])
            gene_map[mut_id] = gene_id
            gene_id_list.add(gene_id)

        sql = "SELECT Gene_Id,Error_Code FROM Gene WHERE Gene_id in (%s)" % ",".join(gene_id_list)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in calculateAnnotationRate: %s" % sql)

        error_map = {}
        for row in results:
            error_map[str(row[0])] = str(row[1])

        #print error_map        

        sql = "SELECT Mutation,Template,Error_Code FROM RS_Mutation_Template WHERE Mutation in (%s)" % ",".join(mut_id_list)
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
            if row[2] == None:
                anno_mut_id_list.add(mut_id)
                all_error_muts.discard(mut_id)
                fine_templates.add(row[1])
            elif mut_id not in anno_mut_id_list:
                all_error_muts.add(mut_id)

        sql = "SELECT Gene,Sequence_Identity,Template_Id FROM Template WHERE Gene in (%s)" % ",".join(gene_id_list)
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
            
        #print "Total Mutations: ",len(mut_id_list)
        #print "Mutations with Annotation: ",len(muts_with_anno)
        #print "Mutations without template: ",len(muts_without_template)
        #print "Mutations with unknown gene error: ",len(muts_with_unknown_gene_error)
        #print "Mutations, which are mapped to gaps in all templates: ",len(all_error_muts)
        #print "Annotation-Rate: ",float(len(muts_with_anno))/float(len(mut_id_list))

        return (len(mut_id_list),len(muts_with_anno),len(muts_without_template),len(muts_with_unknown_gene_error),len(all_error_muts),float(len(muts_with_anno))/float(len(mut_id_list)),number_of_templates,len(greater90),len(greater95),len(greater98),len(greater99),muts_with_temp_98)
    else:
        return (0,0,0,0,0,0,0,0,0,0,0,0)
