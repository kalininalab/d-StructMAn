#TODO replace all sql query, similiar to first big select in mindistout
import pymysql as MySQLdb
import pdbParser as pdb
import sys
import os
import time
import templateFiltering
import rin
from operator import itemgetter
import multiprocessing
from Bio.PDB import *
from Bio.SubsMat import MatrixInfo
import numpy as np
import re

short_distance_threshold = 5.0
long_distance_threshold = 8.0
surface_threshold = 0.16

# PDB Three-letter-codes (spaces stripped)
# Generated from http://ligand-expo.rcsb.org/dictionaries/Components-smiles-oe.smi
# Using: grep -P "^\[?(Li|Be|Na|Mg|Al|K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|Fr|Ra|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs|Cn)[+-]?[0-9]*\]?\t" Components-smiles-oe.smi
# Does not include metalloids
metal_atoms = set(['0BE','3CO','3NI','4MO','4PU','4TI','6MO','AG','AL',
                    'AM','AU','AU3','BA','BS3','CA','CD','CE','CF',
                    'CO','CR','CS','CU','CU1','CU3','DY','ER3','EU','EU3',
                    'FE','FE2','GA','GD','GD3','HG','HO','HO3','IN','IR','IR3','K',
                    'LA','LI','LU','MG','MN','MN3','MO','NA','NI','OS','OS4',
                    'PB','PD','PR','PT','PT4','RB','RE','RH','RH3','RU','SM','SR',
                    'TA0','TB','TH','TL','U1','V','W','Y1','YB','YB2','YT3','ZCM',
                    'ZN','ZN2','ZR'])

# PDB Three-letter-codes (spaces stripped)
# Generated from http://ligand-expo.rcsb.org/dictionaries/Components-smiles-oe.smi
# Using: grep -P "^\[?(As|Si|F|I|Cl|Sb|Te|At|Br)[+-][0-9]*\]?\t" Components-smiles-oe.smi
# Includes (unused/obsolete) -group entries (ending in O) in addition to ions
ion_atoms = set(['BR', 'BRO', 'CL', 'CLO', 'F', 'FLO', 'IDO', 'IOD', 'SB'])

chem_dist_matrix =  {'A': {'A': 0.0, 'C': 0.996, 'E': 0.544, 'D': 0.637, 'G': 0.663, 'F': 0.663, 'I': 0.626, 'H': 0.711, 'K': 0.768, 'M': 0.616, 'L': 0.435, 'N': 0.708, 'Q': 0.641, 'P': 0.949, 'S': 0.634, 'R': 0.977, 'T': 0.657, 'W': 0.985, 'V': 0.646, 'Y': 0.973}, 'C': {'A': 0.996, 'C': 0.0, 'E': 1.051, 'D': 0.804, 'G': 0.901, 'F': 0.859, 'I': 0.856, 'H': 0.595, 'K': 1.153, 'M': 0.651, 'L': 1.093, 'N': 0.687, 'Q': 0.753, 'P': 1.184, 'S': 0.744, 'R': 1.028, 'T': 0.737, 'W': 0.97, 'V': 0.82, 'Y': 0.835}, 'E': {'A': 0.544, 'C': 1.051, 'E': 0.0, 'D': 0.416, 'G': 0.923, 'F': 0.743, 'I': 0.892, 'H': 0.545, 'K': 0.647, 'M': 0.578, 'L': 0.724, 'N': 0.647, 'Q': 0.532, 'P': 0.994, 'S': 0.811, 'R': 0.897, 'T': 0.844, 'W': 0.882, 'V': 0.96, 'Y': 0.973}, 'D': {'A': 0.637, 'C': 0.804, 'E': 0.416, 'D': 0.0, 'G': 0.66, 'F': 0.68, 'I': 0.835, 'H': 0.361, 'K': 0.648, 'M': 0.58, 'L': 0.803, 'N': 0.291, 'Q': 0.385, 'P': 0.747, 'S': 0.555, 'R': 0.793, 'T': 0.64, 'W': 0.76, 'V': 0.886, 'Y': 0.744}, 'G': {'A': 0.663, 'C': 0.901, 'E': 0.923, 'D': 0.66, 'G': 0.0, 'F': 0.813, 'I': 0.814, 'H': 0.82, 'K': 0.974, 'M': 0.902, 'L': 0.827, 'N': 0.576, 'Q': 0.78, 'P': 0.629, 'S': 0.452, 'R': 1.081, 'T': 0.601, 'W': 1.017, 'V': 0.812, 'Y': 0.875}, 'F': {'A': 0.663, 'C': 0.859, 'E': 0.743, 'D': 0.68, 'G': 0.813, 'F': 0.0, 'I': 0.414, 'H': 0.53, 'K': 0.775, 'M': 0.442, 'L': 0.439, 'N': 0.656, 'Q': 0.607, 'P': 0.732, 'S': 0.723, 'R': 0.871, 'T': 0.666, 'W': 0.379, 'V': 0.625, 'Y': 0.509}, 'I': {'A': 0.626, 'C': 0.856, 'E': 0.892, 'D': 0.835, 'G': 0.814, 'F': 0.414, 'I': 0.0, 'H': 0.673, 'K': 0.741, 'M': 0.602, 'L': 0.382, 'N': 0.717, 'Q': 0.611, 'P': 0.9, 'S': 0.602, 'R': 0.754, 'T': 0.469, 'W': 0.733, 'V': 0.239, 'Y': 0.578}, 'H': {'A': 0.711, 'C': 0.595, 'E': 0.545, 'D': 0.361, 'G': 0.82, 'F': 0.53, 'I': 0.673, 'H': 0.0, 'K': 0.669, 'M': 0.346, 'L': 0.758, 'N': 0.365, 'Q': 0.299, 'P': 0.883, 'S': 0.598, 'R': 0.684, 'T': 0.586, 'W': 0.602, 'V': 0.736, 'Y': 0.579}, 'K': {'A': 0.768, 'C': 1.153, 'E': 0.647, 'D': 0.648, 'G': 0.974, 'F': 0.775, 'I': 0.741, 'H': 0.669, 'K': 0.0, 'M': 0.844, 'L': 0.702, 'N': 0.604, 'Q': 0.412, 'P': 0.883, 'S': 0.656, 'R': 0.383, 'T': 0.605, 'W': 0.879, 'V': 0.777, 'Y': 0.71}, 'M': {'A': 0.616, 'C': 0.651, 'E': 0.578, 'D': 0.58, 'G': 0.902, 'F': 0.442, 'I': 0.602, 'H': 0.346, 'K': 0.844, 'M': 0.0, 'L': 0.639, 'N': 0.639, 'Q': 0.534, 'P': 1.024, 'S': 0.762, 'R': 0.903, 'T': 0.725, 'W': 0.637, 'V': 0.698, 'Y': 0.745}, 'L': {'A': 0.435, 'C': 1.093, 'E': 0.724, 'D': 0.803, 'G': 0.827, 'F': 0.439, 'I': 0.382, 'H': 0.758, 'K': 0.702, 'M': 0.639, 'L': 0.0, 'N': 0.8, 'Q': 0.682, 'P': 0.867, 'S': 0.729, 'R': 0.894, 'T': 0.665, 'W': 0.778, 'V': 0.53, 'Y': 0.786}, 'N': {'A': 0.708, 'C': 0.687, 'E': 0.647, 'D': 0.291, 'G': 0.576, 'F': 0.656, 'I': 0.717, 'H': 0.365, 'K': 0.604, 'M': 0.639, 'L': 0.8, 'N': 0.0, 'Q': 0.304, 'P': 0.675, 'S': 0.339, 'R': 0.635, 'T': 0.418, 'W': 0.744, 'V': 0.735, 'Y': 0.555}, 'Q': {'A': 0.641, 'C': 0.753, 'E': 0.532, 'D': 0.385, 'G': 0.78, 'F': 0.607, 'I': 0.611, 'H': 0.299, 'K': 0.412, 'M': 0.534, 'L': 0.682, 'N': 0.304, 'Q': 0.0, 'P': 0.849, 'S': 0.446, 'R': 0.447, 'T': 0.413, 'W': 0.737, 'V': 0.628, 'Y': 0.57}, 'P': {'A': 0.949, 'C': 1.184, 'E': 0.994, 'D': 0.747, 'G': 0.629, 'F': 0.732, 'I': 0.9, 'H': 0.883, 'K': 0.883, 'M': 1.024, 'L': 0.867, 'N': 0.675, 'Q': 0.849, 'P': 0.0, 'S': 0.734, 'R': 1.034, 'T': 0.805, 'W': 0.734, 'V': 1.021, 'Y': 0.676}, 'S': {'A': 0.634, 'C': 0.744, 'E': 0.811, 'D': 0.555, 'G': 0.452, 'F': 0.723, 'I': 0.602, 'H': 0.598, 'K': 0.656, 'M': 0.762, 'L': 0.729, 'N': 0.339, 'Q': 0.446, 'P': 0.734, 'S': 0.0, 'R': 0.662, 'T': 0.189, 'W': 0.924, 'V': 0.539, 'Y': 0.639}, 'R': {'A': 0.977, 'C': 1.028, 'E': 0.897, 'D': 0.793, 'G': 1.081, 'F': 0.871, 'I': 0.754, 'H': 0.684, 'K': 0.383, 'M': 0.903, 'L': 0.894, 'N': 0.635, 'Q': 0.447, 'P': 1.034, 'S': 0.662, 'R': 0.0, 'T': 0.555, 'W': 0.939, 'V': 0.735, 'Y': 0.626}, 'T': {'A': 0.657, 'C': 0.737, 'E': 0.844, 'D': 0.64, 'G': 0.601, 'F': 0.666, 'I': 0.469, 'H': 0.586, 'K': 0.605, 'M': 0.725, 'L': 0.665, 'N': 0.418, 'Q': 0.413, 'P': 0.805, 'S': 0.189, 'R': 0.555, 'T': 0.0, 'W': 0.883, 'V': 0.389, 'Y': 0.56}, 'W': {'A': 0.985, 'C': 0.97, 'E': 0.882, 'D': 0.76, 'G': 1.017, 'F': 0.379, 'I': 0.733, 'H': 0.602, 'K': 0.879, 'M': 0.637, 'L': 0.778, 'N': 0.744, 'Q': 0.737, 'P': 0.734, 'S': 0.924, 'R': 0.939, 'T': 0.883, 'W': 0.0, 'V': 0.932, 'Y': 0.474}, 'V': {'A': 0.646, 'C': 0.82, 'E': 0.96, 'D': 0.886, 'G': 0.812, 'F': 0.625, 'I': 0.239, 'H': 0.736, 'K': 0.777, 'M': 0.698, 'L': 0.53, 'N': 0.735, 'Q': 0.628, 'P': 1.021, 'S': 0.539, 'R': 0.735, 'T': 0.389, 'W': 0.932, 'V': 0.0, 'Y': 0.695}, 'Y': {'A': 0.973, 'C': 0.835, 'E': 0.973, 'D': 0.744, 'G': 0.875, 'F': 0.509, 'I': 0.578, 'H': 0.579, 'K': 0.71, 'M': 0.745, 'L': 0.786, 'N': 0.555, 'Q': 0.57, 'P': 0.676, 'S': 0.639, 'R': 0.626, 'T': 0.56, 'W': 0.474, 'V': 0.695, 'Y': 0.0}}

MobiDB_map = {'D_PA':'Polyampholite','D_WC':'Weak polyampholie','D_NPE':'D_NPE','D_PPE':'D_PPE'}

def binningSelect(keys,rows,table,db,cursor):
    #Use this for huge selects, the first entry of rows has to be the key by which the entries are selected
    key_name = rows[0]
    singletons,bins = binning(keys)

    if len(singletons) > 0:
        in_rows = {key_name:singletons}
        total_results = list(select(db,cursor,rows,table,in_rows=in_rows))
    else:
        total_results = []

    for ids,min_id,max_id in bins:
        between_rows = {key_name:(min_id,max_id)}

        results = select(db,cursor,rows,table,between_rows=between_rows)
        for row in results:
            if not row[0] in ids:
                continue
            total_results.append(row)

    return total_results

def select(db,cursor,rows,table,between_rows={},in_rows={},equals_rows={},null_columns=set(),n_trials=5):
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
            where_parts.append(in_row + ' IN (%s)' % ','.join(['%s']*len(in_rows[in_row]))) #There have to be as many %s placeholders in the statement as there are parameters for the IN clasue

        for bet_row in between_rows:
            (low,high) = between_rows[bet_row]
            params.append(low)
            params.append(high)
            where_parts.append(bet_row + ' BETWEEN %s AND %s')

        where_str = ' WHERE %s' % ' AND '.join(where_parts)

    statement = 'SELECT %s FROM %s%s' % (row_str,table,where_str)

    n = 0
    while n < n_trials: #Repeat the querry if fails for n_trials times
        try:
            cursor.execute(statement,params)
            results = cursor.fetchall()
            db.commit()
            break
        except:
            n+=1
    if n == n_trials:
        [e,f,g] = sys.exc_info()
        raise NameError('Invalid Select: %s,\nParam size:%s\n%s' % (statement[:500],str(len(params)),f))
    return results

def insert(db,cursor,table,columns,values,value_threshold=10000,n_trials=5):
    params = []

    columns_str = ','.join(columns)

    parts = []
    n=0
    value_strs = []

    if len(values) == 0:
        return

    for value_list in values:
        for value in value_list:
            params.append(value)
        value_str_part = '(%s)' % ','.join(['%s']*len(value_list))
        value_strs.append(value_str_part)
        n+=1
        if n == value_threshold:
            n = 0
            value_str = ','.join(value_strs)
            parts.append((value_str,params))
            value_strs = []
            params = []
    value_str = ','.join(value_strs)
    parts.append((value_str,params))

    for value_str,params in parts:
        statement = 'INSERT IGNORE INTO %s (%s) VALUES %s' % (table,columns_str,value_str)

        n = 0
        while n < n_trials: #Repeat the querry if fails for n_trials times
            try:
                cursor.execute(statement,params)
                db.commit()
                break
            except:
                [e,f,g] = sys.exc_info()
                n+=1
        if n == n_trials:
            raise NameError('Invalid Insert: %s,\n%s\n%s' % (statement[:500],e,f))
    return

def update(db,cursor,table,columns,values,value_threshold=10000):
    #Updating with an insert statement (see: https://stackoverflow.com/questions/20255138/sql-update-multiple-records-in-one-query)
    params = []

    column_str = ','.join(columns)

    update_str_parts = []
    for column in columns:
        update_str_parts.append('%s=VALUES(%s)' % (column,column))

    update_str = ','.join(update_str_parts)

    parts = []
    n=0
    value_strs = []

    if len(values) == 0:
        return

    for value_list in values:

        for value in value_list:
            params.append(value)
        value_str_part = '(%s)' % ','.join(['%s']*len(value_list))
        value_strs.append(value_str_part)
        n+=1
        if n == value_threshold:
            n = 0
            value_str = ','.join(value_strs)
            parts.append((value_str,params))
            value_strs = []
            params = []
            min_primary_value = None
            max_primary_value = None
    value_str = ','.join(value_strs)
    parts.append((value_str,params))

    for value_str,params in parts:
        statement = 'INSERT IGNORE INTO %s (%s) VALUES %s ON DUPLICATE KEY UPDATE %s' % (table,column_str,value_str,update_str)
        try:
            cursor.execute(statement,params)
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError('Invalid Update: %s,\n%s' % (statement[:500],f))

def getGeneScoreDict(gene_id_list,session_id,db,cursor,includeSequence=False):
    t0 = time.time()
    if len(gene_id_list) == 0:
        return {}
    max_g = max(gene_id_list)
    min_g = min(gene_id_list)
    #print max_g,min_g
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

    if not includeSequence:
        rows = ['Gene_Id','Uniprot_Ac','Genbank_Protein_Accession_Number','Uniprot_Id','Error_Code','Error','Species']
    else:
        rows = ['Gene_Id','Uniprot_Ac','Genbank_Protein_Accession_Number','Uniprot_Id','Error_Code','Error','Species','Sequence']
    table = 'Gene'
    results = binningSelect(gene_id_list,rows,table,db,cursor)
    
    gene_score_dict = {}
    for row in results:
        g_id = row[0]
        if g_id in gene_id_list and g_id in score_dict:
            if not includeSequence:
                gene_score_dict[g_id] = (row[1],row[2],row[3],row[4],row[5],row[6],score_dict[g_id])
            else:
                gene_score_dict[g_id] = (row[1],row[2],row[3],row[4],row[5],row[6],score_dict[g_id],row[7])
    t2 = time.time()

    #print "Time for getGeneScoreDict1: ",t1-t0
    #print "Time for getGeneScoreDict2: ",t2-t1

    #print gene_score_dict
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
        raise NameError("Error in geneScan: %s,\n%s" % (sql,f))
    
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
def geneCheck(genes_aac_list,species_map,session_id,db,cursor):
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
    #print genes_aac_list.keys()
    for row in results:
        #print row
        gene_id = row[0]
        u_id = row[1]
        if not u_id in genes_aac_list:
            continue
        stored_genes[u_id] = gene_id
        gene_id_list.add(gene_id)
        if gene_id > max_g_id:
            max_g_id = gene_id
    #print session_id
    t2 = time.time()

    t3 = time.time()

    #print session_id
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
                value_strs.append("('%s','%s','%s','%s','%s')" % (gene.replace("'","\\'"),','.join(genes_aac_list[gene][1]),genes_aac_list[gene][0],str(session_id),species_map[gene]))
            else:
                value_strs.append("('%s','%s','%s','%s',NULL)" % (gene.replace("'","\\'"),','.join(genes_aac_list[gene][1]),genes_aac_list[gene][0],str(session_id)))

        size = len(value_strs)
        #print size
        if size > 0:
            bound = 50000
            if size > bound:
                m = size//bound
                rest = size%bound
                if rest == 0:
                    part_size = size//m
                else:
                    m += 1
                    part_size = size//m
                parts = []
                for i in range(0,m-1):
                    parts.append(value_strs[(i*part_size):((i+1)*part_size)])
                parts.append(value_strs[(m-1)*part_size:])
            else:
                parts = [value_strs]
            for part in parts:

                sql = ("""INSERT INTO Gene(Uniprot_Ac,Genbank_Protein_Accession_Number,
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
    #print session_id
    for gene in new_gene_map:
        value_strs.append("('%s','%s')" % (str(new_gene_map[gene]),str(session_id)))

    for gene in stored_genes:
        value_strs.append("('%s','%s')" % (str(stored_genes[gene]),str(session_id)))

    #print(value_strs)

    sql = ("""INSERT IGNORE INTO RS_Gene_Session(
                Gene, Session)
                VALUES %s""") % ','.join(value_strs)
    try:
        cursor.execute(sql)
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Couldn't insert geneCheck: %s,%s" % (sql,f))

    #Return the mappings uniprot-id:gene_id(and 'more restrictive?')
    t8 = time.time()

    '''
    print "Time for geneCheck Part 1: %s" % str(t1-t0)
    print "Time for geneCheck Part 2: %s" % str(t2-t1)
    print "Time for geneCheck Part 3: %s" % str(t3-t2)
    print "Time for geneCheck Part 4: %s" % str(t4-t3)
    print "Time for geneCheck Part 5: %s" % str(t5-t4)
    print "Time for geneCheck Part 6: %s" % str(t6-t5)
    print "Time for geneCheck Part 7: %s" % str(t7-t6)
    print "Time for geneCheck Part 8: %s" % str(t8-t7)
    '''
    return stored_genes,new_gene_map,gene_id_list,new_gene_ids


def getLastAutoInc(db,cursor):
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


def getLigandDict(s_ids,db,cursor):
    t0 = time.time()
    if s_ids != None:
        if len(s_ids) == 0:
            return {}
        max_s = max(s_ids)
        min_s = min(s_ids)
        sql = "SELECT Ligand,Structure,Chain,Residue FROM RS_Ligand_Structure WHERE Structure BETWEEN %s AND %s" % (str(min_s),str(max_s))
             
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in getLigands: %s" % sql)
            db.rollback()

        if results == ():
            return {}

        structure_dict = {}
        ligand_id_list = set()
        for row in results:
            s_id = row[1]
            if not s_id in s_ids:
                continue
            if not s_id in structure_dict:
                structure_dict[s_id] = {row[0]:["Ligand","",row[3],row[2]]}
            elif not row[0] in structure_dict[s_id]:
                structure_dict[s_id][row[0]] = ["Ligand","",row[3],row[2]]
            ligand_id_list.add(row[0])
    else:
        sql = "SELECT Ligand,Structure,Chain,Residue FROM RS_Ligand_Structure"
             
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("Error in getLigands: %s" % sql)
            db.rollback()

        if results == ():
            return {}

        structure_dict = {}
        ligand_id_list = set()
        for row in results:
            s_id = row[1]
            if not s_id in structure_dict:
                structure_dict[s_id] = {row[0]:["Ligand","",row[3],row[2]]}
            elif not row[0] in structure_dict[s_id]:
                structure_dict[s_id][row[0]] = ["Ligand","",row[3],row[2]]
            ligand_id_list.add(row[0])

    t1 = time.time()


    max_l = max(ligand_id_list)
    min_l = min(ligand_id_list)
    #print max_l,min_l
    if len(ligand_id_list) > 0:
        sql = "SELECT Ligand_Id,Name FROM Ligand WHERE Ligand_Id BETWEEN %s AND %s" % (str(min_l),str(max_l))
             
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


        for s_id in structure_dict:
            for lig_id in structure_dict[s_id]:
                structure_dict[s_id][lig_id][1] = lig_map[lig_id]
    t2 = time.time()

    #print "Time for getLigandDict1: ",t1-t0
    #print "Time for getLigandDict2: ",t2-t1

    return structure_dict

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
def mutationCheck(gene_aaclist_map,stored_genes,stored_gene_ids,new_genes,new_gene_ids,database_session,tag_map,db,cursor,config,residue_id_backmap):
    #structure of stored_genes: {Uniprot-Id:gene_id}
    #structure of new_genes: {Uniprot-Id:gene_id}
    verbose = config.verbose

    gene_mut_map_new = {}
    stored_gene_new_pos = {}
    stored_gene_stored_pos = {}
    gene_mut_map_twins = {}
    t0 = time.time()

    stored_mutation_map = {}
    
    if len(stored_gene_ids) > 0:
        rows = ['Gene','Amino_Acid_Change','Mutation_Id']
        table = 'Mutation'

        results = binningSelect(stored_gene_ids,rows,table,db,cursor)

        for row in results:
            gene_id = row[0]
            if not gene_id in stored_gene_ids:
                continue
            aacs = row[1].split(",")
            aacbase = aacs[0]
            pos = aacbase[1:]
            mutation_id = row[2]

            if not gene_id in stored_mutation_map:
                stored_mutation_map[gene_id] = set()
                stored_mutation_map[gene_id].add(pos)
            else:
                stored_mutation_map[gene_id].add(pos)

    t1 = time.time()
    for gene in gene_aaclist_map:
        aaclist = gene_aaclist_map[gene][2]

        if gene in stored_genes:
            gene_id = stored_genes[gene]
            for aac in aaclist:
                aacbase = aac[:-1]
                pos = aac[1:-1]
                if gene_id in stored_mutation_map:
                    if pos in stored_mutation_map[gene_id]:
                        if gene not in stored_gene_stored_pos:
                            stored_gene_stored_pos[gene] = set()
                        stored_gene_stored_pos[gene].add(aacbase)
                    else:
                        if gene not in stored_gene_new_pos:
                            stored_gene_new_pos[gene] = set()
                        stored_gene_new_pos[gene].add(aacbase)
                else:
                    #print 'bingo'
                    if gene not in stored_gene_new_pos:
                        stored_gene_new_pos[gene] = set()
                    stored_gene_new_pos[gene].add(aacbase)
                    
        else:
            for aac in aaclist:
                aacbase = aac[:-1]
                pos = aac[1:-1]

                if gene not in gene_mut_map_new:
                    gene_mut_map_new[gene] = set()
                gene_mut_map_new[gene].add(aacbase)

    t2 = time.time()

    values = []
    for gene in gene_mut_map_new:
        gene_id = new_genes[gene]
        for aacbase in gene_mut_map_new[gene]:
            if gene in residue_id_backmap:
                values.append((gene_id,aacbase,residue_id_backmap[gene][int(aacbase[1:])]))
            else:
                values.append((gene_id,aacbase,None))

    for gene in stored_gene_new_pos:
        gene_id = stored_genes[gene]
        for aacbase in stored_gene_new_pos[gene]:
            if gene in residue_id_backmap:
                values.append((gene_id,aacbase,residue_id_backmap[gene][int(aacbase[1:])]))
            else:
                values.append((gene_id,aacbase,None))

    t3 = time.time()

    if len(values) > 0:
        columns = ['Gene','Amino_Acid_Change','Residue_Id']
        insert(db,cursor,'Mutation',columns,values,value_threshold=50000)

    t4 = time.time()

    t5 = time.time()

    t6 = time.time()

    update_map = {}

    rows = ['Gene','Amino_Acid_Change','Mutation_Id']
    table = 'Mutation'

    fused_gene_ids = new_gene_ids | stored_gene_ids

    results = binningSelect(fused_gene_ids,rows,table,db,cursor)

    t7 = time.time()

    for row in results:
        gene_id = row[0]
        if not (gene_id in fused_gene_ids):
            #print gene_id
            continue
        aacbase = row[1]
        pos = aacbase[1:]

        mutation_id = row[2]
        if not gene_id in update_map:
            update_map[gene_id] = {pos:mutation_id}
        else:
            update_map[gene_id][pos] = mutation_id

    t8 = time.time()

    mut_new_aa_map = {}
    gene_mut_map_new_id = {}
    stored_gene_new_pos_id = {}
    stored_gene_stored_pos_id = {}
    tag_update_map = {}

    for gene in gene_aaclist_map:
        aaclist = gene_aaclist_map[gene][2]
        if gene in stored_genes:
            gene_id = stored_genes[gene]
        else:
            gene_id = new_genes[gene]       

        for aac in aaclist:
            aac_base = aac[:-1]
            pos = aac_base[1:]
            aa2 = aac[-1]
            mut_id = update_map[gene_id][pos]
            if not mut_id in mut_new_aa_map:
                mut_new_aa_map[mut_id] = set()
                mut_new_aa_map[mut_id].add(aac[-1])
            else:
                mut_new_aa_map[mut_id].add(aac[-1])
            if gene in tag_map:
                if aac in tag_map[gene]:
                    tag_update_map[(mut_id,aa2)] = tag_map[gene][aac]
                else:
                    tag_update_map[(mut_id,aa2)] = None
            else:
                tag_update_map[(mut_id,aa2)] = None
            if gene in gene_mut_map_new and aac_base in gene_mut_map_new[gene]:
                if not gene in gene_mut_map_new_id:
                    gene_mut_map_new_id[gene] = (gene_id,{aac_base:mut_id})
                else:
                    gene_mut_map_new_id[gene][1][aac_base] = mut_id
            elif gene in stored_gene_new_pos:
                if not gene in stored_gene_new_pos_id:
                    stored_gene_new_pos_id[gene] = (gene_id,{aac_base:mut_id})
                else:
                    stored_gene_new_pos_id[gene][1][aac_base] = mut_id
            elif gene in stored_gene_stored_pos:
                if not gene in stored_gene_stored_pos_id:
                    stored_gene_stored_pos_id[gene] = (gene_id,{aac_base:mut_id})
                else:
                    stored_gene_stored_pos_id[gene][1][aac_base] = mut_id

    t9 = time.time()

    values = []

    for mut_id in mut_new_aa_map:
        for new_aa in mut_new_aa_map[mut_id]:
            if tag_update_map[(mut_id,new_aa)] == None:
                values.append((database_session,mut_id,new_aa,'NULL'))
            else:
                values.append((database_session,mut_id,new_aa,tag_update_map[(mut_id,new_aa)]))

    t10 = time.time()

    if verbose:
        print("Time for mutationCheck Part 1: %s" % (str(t1-t0)))
        print("Time for mutationCheck Part 2: %s" % (str(t2-t1)))
        print("Time for mutationCheck Part 3: %s" % (str(t3-t2)))
        print("Time for mutationCheck Part 4: %s" % (str(t4-t3)))
        print("Time for mutationCheck Part 5: %s" % (str(t5-t4)))
        print("Time for mutationCheck Part 6: %s" % (str(t6-t5)))
        print("Time for mutationCheck Part 7: %s" % (str(t7-t6)))
        print("Time for mutationCheck Part 8: %s" % (str(t8-t7)))
        print("Time for mutationCheck Part 9: %s" % (str(t9-t8)))
        print("Time for mutationCheck Part 10: %s" % (str(t10-t9)))

    if len(values) > 10000:
        process = multiprocessing.Process(target = backgroundInsertMS,args=(values,config))
        try:
            process.start()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in mutationCheck: %s" % (f))
        
        t11 = time.time()
        if verbose:
            print("Time for mutationCheck Part 11: %s" % (str(t11-t10)))
        return gene_mut_map_new_id,stored_gene_new_pos_id,process
    else:
        insert(db,cursor,'RS_Mutation_Session',['Session','Mutation','New_AA','Tag'],values,value_threshold=1000)
        t11 = time.time()
        if verbose:
            print("Time for mutationCheck Part 11: %s" % (str(t11-t10)))
        return gene_mut_map_new_id,stored_gene_new_pos_id,None

def backgroundInsertMS(values,config):
    db = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,config.db_name)
    cursor = db.cursor()
    insert(db,cursor,'RS_Mutation_Session',['Session','Mutation','New_AA','Tag'],values,value_threshold=1000)
    db.close()
    return

#called by serializedPipeline
def addIupred(iupred_map,gene_mut_map_new,stored_gene_new_pos,config):
    process = multiprocessing.Process(target = backgroundIU,args=(iupred_map,gene_mut_map_new,stored_gene_new_pos,config))
    process.start()
    return process

def backgroundIU(iupred_map,gene_mut_map_new,stored_gene_new_pos,config):
    #structure of gene_mut_map_new: {Uniprot_Ac:(gene_id,{AAC_Base:mutation_id})}
    #structure of iupred_map: {Uniprot_Ac:([regions],{Position:(aa1,iupred_score)})}

    binsize = 5000
    bins = set()
    max_mut_id = 0
    min_mut_id = None
    value_strs = {}
    value_strs_glob = {}

    def addToValues(g_m_map,max_mut_id,min_mut_id,bins):
        for u_ac in g_m_map:
            if not u_ac in iupred_map:
                print('Warning: Did not find %s in the iupred map' % u_ac)
                continue
            regions = iupred_map[u_ac][0]
            method = iupred_map[u_ac][2]
            for aac_base in g_m_map[u_ac][1]:
                m_id = g_m_map[u_ac][1][aac_base]
                aac_pos = int(aac_base[1:])
                if not aac_pos in iupred_map[u_ac][1]:
                    continue
                if method == 'MobiDB3.0'or method == 'mobidb-lite':
                    pos_region_type = 'globular'
                else:
                    pos_region_type = 'disorder'
                for [a,b,region_type] in regions:
                    if aac_pos > int(a) and aac_pos < int(b):
                        pos_region_type = region_type
                aa1,iupred_score = iupred_map[u_ac][1][aac_pos]

                bin_number = m_id//binsize
                bins.add(bin_number)
                if m_id > max_mut_id:
                    max_mut_id = m_id
                if min_mut_id == None or m_id < min_mut_id:
                    min_mut_id = m_id
                if not bin_number in value_strs:
                    value_strs[bin_number] = []
                if not bin_number in value_strs_glob:
                    value_strs_glob[bin_number] = []
                value_strs[bin_number].append("WHEN '%s' THEN '%s'" % (str(m_id),str(iupred_score)))

                value_strs_glob[bin_number].append("WHEN '%s' THEN '%s'" % (str(m_id),pos_region_type))
        return max_mut_id,min_mut_id,bins

    db = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,config.db_name)
    cursor = db.cursor()

    max_mut_id,min_mut_id,bins = addToValues(gene_mut_map_new,max_mut_id,min_mut_id,bins)
    max_mut_id,min_mut_id,bins = addToValues(stored_gene_new_pos,max_mut_id,min_mut_id,bins)

    if not len(value_strs) == 0:
        min_max_tuples = []

        for bin_number in bins:
            min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize))

        
        for (min_m,max_m) in min_max_tuples:
            bin_number = min_m//binsize
            part = value_strs[bin_number]
            glob_part = value_strs_glob[bin_number]

            sql = "UPDATE Mutation SET IUPRED = CASE Mutation_Id %s ELSE IUPRED END WHERE Mutation_Id BETWEEN %s AND %s" % (" ".join(part),min_m,max_m)
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in addIupred: %s,%s" % (sql,f))

            sql = "UPDATE Mutation SET IUPRED_Glob = CASE Mutation_Id %s ELSE IUPRED_Glob END WHERE Mutation_Id BETWEEN %s AND %s" % (" ".join(glob_part),min_m,max_m)
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in addIupred: %s,%s" % (sql,f))
    db.close()
    return

#called by serializedPipeline
def updateMutations(mutation_updates,db,cursor):
    binsize = 500
    bins = set()
    max_mut_id = 0
    min_mut_id = None

    value_strs = {}
    m_ids = []
    for (m_id,new_aac_base) in mutation_updates:
        
        m_ids.append(str(m_id))
        bin_number = m_id//binsize
        bins.add(bin_number)

        if not bin_number in value_strs:
            value_strs[bin_number] = []
        value_strs[bin_number].append("WHEN '%s' THEN '%s'" % (str(m_id),str(new_aac_base)))

        if m_id > max_mut_id:
            max_mut_id = m_id
        if min_mut_id == None or m_id < min_mut_id:
            min_mut_id = m_id

    if not len(m_ids) == 0:
        min_max_tuples = []

        for bin_number in bins:
            min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize))

        
        for (min_mut_id,max_mut_id) in min_max_tuples:
            bin_number = min_mut_id//binsize
            part = value_strs[bin_number]
        
            sql = "UPDATE Mutation SET Amino_Acid_Change = CASE Mutation_Id %s ELSE Amino_Acid_Change END WHERE Mutation_Id BETWEEN %s AND %s" % (" ".join(part),str(min_mut_id),str(max_mut_id))
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in updateMutations: %s,%s" % (sql[:1000],f))

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
            coverage = row[3]
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
            template = [pdb_id,seq_id,target_chain,coverage,resolution,[],r_value,score,sub_lig_dist,sub_ano,original_target_chain,original_chains]
            oligo_set = set([])
            for oligo in oligos:
                oligo_set.add(oligo)
            if not gene_id in gene_template_map:
                gene_template_map[gene_id] = ({template_id:template},{pdb_id:oligo_set})
            else:
                gene_template_map[gene_id][0][template_id] = template
                gene_template_map[gene_id][1][pdb_id] = oligo_set
    return gene_template_map

def getStructureDict(s_ids,db,cursor):
    t0 = time.time()
    structure_ligand_dict = getLigandDict(s_ids,db,cursor)
    if not s_ids == None:
        if len(s_ids) == 0:
            return {}

    if not s_ids == None:
        sql = "SELECT Structure_Id,PDB,Chain,Resolution,Homooligomer,Chains FROM Structure WHERE Structure_Id BETWEEN %s AND %s" % (str(min(s_ids)),str(max(s_ids)))
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("NameError in getTemplatesFromIdList: %s" % sql)
        structure_results = {}
        for row in results:
            s_id = row[0]
            if not s_id in s_ids:
                continue
            pdb_id = row[1]
            chain = row[2]
            resolution = row[3]
            oligos = row[4]
            chains = row[5]
            if s_id in structure_ligand_dict:
                ligands = list(structure_ligand_dict[s_id].values())
            else:
                ligands = []
            """
            iaps = chains.split(";")
            for iap in iaps:
                iap_info = iap.split(":")
                if len(iap_info) > 1:
                    ligands.append([iap_info[1],iap_info[0]])
            """

            structure = [pdb_id,chain,resolution,ligands,oligos,chains]
            structure_results[s_id] = structure

    else:
        sql = "SELECT Structure_Id,PDB,Chain,Resolution,Homooligomer,Chains FROM Structure"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            raise NameError("NameError in getTemplatesFromIdList: %s" % sql)
        structure_results = {}
        for row in results:
            s_id = row[0]
            pdb_id = row[1]
            chain = row[2]
            resolution = row[3]
            oligos = row[4]
            chains = row[5]
            if s_id in structure_ligand_dict:
                ligands = list(structure_ligand_dict[s_id].values())
            else:
                ligands = []
            """
            iaps = chains.split(";")
            for iap in iaps:
                iap_info = iap.split(":")
                if len(iap_info) > 1:
                    ligands.append([iap_info[1],iap_info[0]])
            """

            structure = [pdb_id,chain,resolution,ligands,oligos,chains]
            structure_results[s_id] = structure
    return structure_results

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
    if len(mutation_id_list) == 0:
        return {}

    rows = ['Mutation_Id','Amino_Acid_Change','Gene','IUPRED','IUPRED_Glob','Residue_Id']
    table = 'Mutation'
    results = binningSelect(mutation_id_list,rows,table,db,cursor)


    mutation_dict= {}

    for row in results:
        mut_id = row[0]
        if mut_id in mutation_id_list:
            mutation_dict[mut_id] = (row[1],row[2],row[3],row[4],row[5])
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
        print(("Couldn't update Session '%s'") % (str(session_id)))
        # Rollback in case there is any NameError
        db.rollback()

#called by serializedPipeline
def insertSession(time,ori_filename,db,cursor):
    sql = "INSERT IGNORE INTO Session(Input_File,Start) VALUES ('%s','%s')" % (ori_filename,str(time))
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
    #print gene_info_map
    ref_value_strs = []
    go_term_ids = set()
    reac_ids = set()
    seq_value_strs = []
    for gene_id in gene_info_map:
        (go_terms,pathways,sequence) = gene_info_map[gene_id]
        for go_id in go_terms:
            go_term_ids.add(go_id)
        for reac_id in pathways:
            reac_ids.add(reac_id)
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
        (go_terms,pathways,sequence) = gene_info_map[gene_id]
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
        (go_terms,pathways,sequence) = gene_info_map[gene_id]
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

    gene_id_strs = [str(x) for x in list(gene_info_map.keys())]


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
        stored_gene_ids.add(stored_genes[gene])
        gene_id_gene_map[stored_genes[gene]] = gene
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
        sequence_map[gene_id_gene_map[gene_id]] = gene_id_sequence_map[gene_id],'Stored','Stored'

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
        gene_id_string = ','.join([str(x) for x in list(gene_error_map.keys())])
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
            m = size//bound
            rest = size%bound
            if rest == 0:
                part_size = size//m
            else:
                m += 1
                part_size = size//m
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
                m = size//bound
                rest = size%bound
                if rest == 0:
                    part_size = size//m
                else:
                    m += 1
                    part_size = size//m
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
    return

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

    proc = multiprocessing.Process(target=backgroundInsertAS, args=(value_strings,max_m,min_m,AS_db,anno_anno_calculation))
    proc.start()

    return stored_annotations,proc

#called by serializedPipeline
def insertComplexes(complex_profiles,pdb_ids,db,cursor):
    value_strs = []
    for pdb in complex_profiles:
        ligand_profiles,metal_profiles,ion_profiles,chain_chain_profiles = complex_profiles[pdb]

        lig_str_parts = []
        for chain,res in ligand_profiles:
            deg,score = ligand_profiles[(chain,res)]
            lig_str_parts.append('%s_%s:%s_%s' % (chain,res,str(deg),str(score)))
        lig_str = ','.join(lig_str_parts)
        
        metal_str_parts = []
        for chain,res in metal_profiles:
            deg,score = metal_profiles[(chain,res)]
            metal_str_parts.append('%s_%s:%s_%s' % (chain,res,str(deg),str(score)))
        metal_str = ','.join(metal_str_parts)

        ion_str_parts = []
        for chain,res in ion_profiles:
            deg,score = ion_profiles[(chain,res)]
            ion_str_parts.append('%s_%s:%s_%s' % (chain,res,str(deg),str(score)))
        ion_str = ','.join(ion_str_parts)

        cc_str_parts = []
        for chain_a,chain_b in chain_chain_profiles:
            deg,score = chain_chain_profiles[(chain_a,chain_b)]
            cc_str_parts.append('%s_%s:%s_%s' % (chain_a,chain_b,str(deg),str(score)))
        cc_str = ','.join(cc_str_parts)

        value_strs.append("('%s','%s','%s','%s','%s')" % (pdb,lig_str,metal_str,ion_str,cc_str))

    if len(value_strs) > 0:
        value_strs_list = []
        i = 0
        part_size = 5000
        while i*part_size < len(value_strs):
            a = i*part_size
            b = min(((i+1)*part_size,len(value_strs)))
            i += 1
            sql = "INSERT IGNORE INTO Complex(PDB,Ligand_Profile,Metal_Profile,Ion_Profile,Chain_Chain_Profile) VALUES %s" % ','.join(value_strs[a:b])
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in insertComplexes: %s" % (f))

    rows = ['PDB','Ligand_Profile','Metal_Profile','Ion_Profile','Chain_Chain_Profile']
    table = 'Complex'
    if len(pdb_ids) > 0:
        results = select(db,cursor,rows,table,in_rows={'PDB':pdb_ids})
    else:
        results = []

    for row in results:
        pdb_id = row[0]
        if pdb_id in complex_profiles:
            continue
        lig_str = row[1]
        metal_str = row[2]
        ion_str = row[3]
        cc_str = row[4]
        lig_profile = {}
        if lig_str != '':
            for lig_str_parts in lig_str.split(','):
                part1,part2 = lig_str_parts.split(':')
                chain,res = part1.split('_')
                deg,score = part2.split('_')
                deg = int(deg)
                score = float(score)
                lig_profile[(chain,res)] = (deg,score)

        metal_profile = {}
        if metal_str != '':
            for metal_str_parts in metal_str.split(','):
                part1,part2 = metal_str_parts.split(':')
                chain,res = part1.split('_')
                deg,score = part2.split('_')
                deg = int(deg)
                score = float(score)
                metal_profile[(chain,res)] = (deg,score)

        ion_profile = {}
        if ion_str != '':
            for ion_str_parts in ion_str.split(','):
                part1,part2 = ion_str_parts.split(':')
                chain,res = part1.split('_')
                deg,score = part2.split('_')
                deg = int(deg)
                score = float(score)
                ion_profile[(chain,res)] = (deg,score)

        cc_profile = {}
        if cc_str != '':
            for cc_str_parts in cc_str.split(','):
                part1,part2 = cc_str_parts.split(':')
                chain,res = part1.split('_')
                deg,score = part2.split('_')
                deg = int(deg)
                score = float(score)
                cc_profile[(chain,res)] = (deg,score)

        complex_profiles[pdb_id] = lig_profile,metal_profile,ion_profile,cc_profile

    return complex_profiles

def getComplexMap(db,cursor,pdb_ids=None):
    table = 'Complex'
    rows = ['PDB','Ligand_Profile','Metal_Profile','Ion_Profile','Chain_Chain_Profile']

    results = select(db,cursor,rows,table)

    complex_map = {}
    for row in results:
        pdb_id = row[0]
        if pdb_ids != None:
            if not pdb_id in pdb_ids:
                continue
        lig_profile_str = row[1]
        metal_profile_str = row[2]
        ion_profile_str = row[3]
        cc_profile_str = row[4]

        lig_profile = {}
        if lig_profile_str != '':
            for lig_profile_part in lig_profile_str.split(','):
                part1,part2 = lig_profile_part.split(':')
                chain,res = part1.split('_')
                deg,score = part2.split('_')

                lig_profile[(chain,res)] = int(deg),float(score)
    
        metal_profile = {}

        if metal_profile_str != '':
            for metal_profile_part in metal_profile_str.split(','):
                part1,part2 = metal_profile_part.split(':')
                chain,res = part1.split('_')
                deg,score = part2.split('_')

                metal_profile[(chain,res)] = int(deg),float(score)

        ion_profile = {}

        if ion_profile_str != '':
            for ion_profile_part in ion_profile_str.split(','):
                part1,part2 = ion_profile_part.split(':')
                chain,res = part1.split('_')
                deg,score = part2.split('_')

                ion_profile[(chain,res)] = int(deg),float(score)

        cc_profile = {}

        if cc_profile_str != '':
            for cc_profile_part in cc_profile_str.split(','):
                part1,part2 = cc_profile_part.split(':')
                chain,chain_b = part1.split('_')
                deg,score = part2.split('_')

                cc_profile[(chain,chain_b)] = int(deg),float(score)

        complex_map[pdb_id] = (lig_profile,metal_profile,ion_profile,cc_profile)

    return complex_map


#called by serializedPipeline
def insertStructures(structurelist,db,cursor,smiles_path,inchi_path,pdb_path,structure_id_map):
    #structure_id_map contains only the structures of stored alignments, which means only for inputs, where the target protein and the template structure where already stored
    #   afterwards, it contains all structure_ids of all mapped structures
    #structurelist contains all other structures, which means there could be stored structures amoing them
    stored_structures = {}
    structure_dict = {}
    table = 'Structure'
    rows = ['Structure_Id','PDB','Chain','Resolution','Homooligomer','Chains']
    results = select(db,cursor,rows,table)
    
    for row in results:
        s_id = row[0]
        pdb_id = row[1]
        chain = row[2]
        resolution = row[3]
        oligos = row[4]
        chains = row[5]
        stored_structures[(pdb_id,chain)] = s_id #all structures, mapped or not go into this dictionary, this is important for not reinserting residues from interacting structures
        if (pdb_id,chain) in structurelist:
            structure_dict[s_id] = [pdb_id,chain,resolution,None,oligos,chains]
            structure_id_map[(pdb_id,chain)] = s_id
            del structurelist[(pdb_id,chain)]
        elif (pdb_id,chain) in structure_id_map:
            structure_dict[s_id] = [pdb_id,chain,resolution,None,oligos,chains]

    t0 = time.time()
    #structure of structurelist: {(pdb_id,chain): (resolution,oligo_counter,interaction_partner)}
    value_strs = []
    ligand_map = {}
    ligands = set()

    for (pdb_id,chain) in structurelist:
        (resolution,oligos,interaction_partner) = structurelist[(pdb_id,chain)]
        oligos = ''.join(oligos)
        
        chains = []
        ligand_map[(pdb_id,chain)] = []
        for iap in interaction_partner:
            ia_type = iap[0]
            if ia_type == "Ligand":
                ligand_map[(pdb_id,chain)].append(iap)
                ligands.add(iap[1])
            else:
                chains.append("%s:%s" % (iap[1],iap[0]))

        chain_str = ','.join(chains)
        value_strs.append("('%s','%s','%1.2f','%s','%s')" % (pdb_id,chain,float(resolution),oligos,chain_str))

    t1 = time.time()
    if len(value_strs) > 0:
        #For large amounts of templates, a single insert can destroy the connection to the database => Solution: divide large inserts into several smaller ones
        value_strs_list = []
        i = 0
        part_size = 5000
        while i*part_size < len(value_strs):
            a = i*part_size
            b = min(((i+1)*part_size,len(value_strs)))
            i += 1
            sql = "INSERT IGNORE INTO Structure(PDB,Chain,Resolution,Homooligomer,Chains) VALUES %s" % ','.join(value_strs[a:b])
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in insertTemplates: %s" % (f))

    t2 = time.time()

    if len(structurelist) > 0:
        table = 'Structure'
        rows = ['Structure_Id','PDB','Chain','Resolution','Homooligomer','Chains']
        results = select(db,cursor,rows,table)

        for row in results:
            s_id = row[0]
            pdb_id = row[1]
            chain = row[2]
            resolution = row[3]
            oligos = row[4]
            chains = row[5]

            if not (pdb_id,chain) in structurelist:
                continue

            structure_id_map[(pdb_id,chain)] = s_id
            structure_dict[s_id] = [pdb_id,chain,resolution,None,oligos,chains]

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

    t5 = time.time()

    update_ligands = []

    for (pdb_id,chain)in ligand_map:
        iaps = ligand_map[(pdb_id,chain)]
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
                        (smiles,inchi) = pdb.getSI(pdb_id,name,res,chain,pdb_path)
                        update_ligands.append((name,smiles,inchi))
                    value_strs.append("('%s','%s','%s')" % (name,smiles,inchi))

    if len(update_ligands) > 0:
        pdb.updateLigandDB(update_ligands,smiles_path,inchi_path)

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

    for (pdb_id,chain) in structurelist:
        s_id = structure_id_map[(pdb_id,chain)]
        (resolution,oligos,interaction_partner) = structurelist[(pdb_id,chain)]
        ligand_ids = set()
        for iap in interaction_partner:
            ia_type = iap[0]
            if ia_type == "Ligand":
                name = iap[1]
                if name in stored_ligands:
                    ligand_id = stored_ligands[name]
                elif name in new_ligand_map:
                    ligand_id = new_ligand_map[name]
                else:
                    continue
                if ligand_id in ligand_ids:
                    continue
                ligand_ids.add(ligand_id)
                res = iap[2]
                chain = iap[3]
                value_strs.append("('%s','%s','%s','%s')" % (str(ligand_id),str(s_id),str(chain),str(res)))

    t9 = time.time()

    if len(value_strs) > 0:
        sql = "INSERT IGNORE INTO RS_Ligand_Structure(Ligand,Structure,Chain,Residue) VALUES %s" % ','.join(value_strs)
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
    return structure_id_map,stored_structures,structure_dict

#called by serializedPipeline
def insertInteractingChains(interacting_structure_dict,db,cursor,smiles_path,inchi_path,pdb_path):
    if len(interacting_structure_dict) == 0:
       return {}

    interacting_structure_ids = {}

    sql = "SELECT Structure_Id,PDB,Chain FROM Structure"
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in insertTemplates: %s\n%s" % (sql,f))
    for row in results:
        pdb_id = row[1]
        chain = row[2]
        if not (pdb_id,chain) in interacting_structure_dict:
            continue
        s_id = row[0]
        interacting_structure_ids[(pdb_id,chain)] = s_id
        del interacting_structure_dict[(pdb_id,chain)]
    

    t0 = time.time()

    value_strs = []
    ligand_map = {}
    ligands = set()

    for (pdb_id,chain) in interacting_structure_dict:
        (resolution,interaction_partner,oligos) = interacting_structure_dict[(pdb_id,chain)]
        oligos = ''.join(oligos)
        
        chains = []
        ligand_map[(pdb_id,chain)] = []
        for iap in interaction_partner:
            ia_type = iap[0]
            if ia_type == "Ligand":
                ligand_map[(pdb_id,chain)].append(iap)
                ligands.add(iap[1])
            else:
                chains.append("%s:%s" % (iap[1],iap[0]))

        chain_str = ','.join(chains)
        value_strs.append("('%s','%s','%1.2f','%s','%s')" % (pdb_id,chain,float(resolution),oligos,chain_str))

    t1 = time.time()
    if len(value_strs) > 0:
        #For large amounts of templates, a single insert can destroy the connection to the database => Solution: divide large inserts into several smaller ones
        value_strs_list = []
        i = 0
        part_size = 5000
        while i*part_size < len(value_strs):
            a = i*part_size
            b = min(((i+1)*part_size,len(value_strs)))
            i += 1
            sql = "INSERT IGNORE INTO Structure(PDB,Chain,Resolution,Homooligomer,Chains) VALUES %s" % ','.join(value_strs[a:b])
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in insertTemplates: %s" % (f))

    t2 = time.time()

    if len(interacting_structure_dict) > 0:
        sql = "SELECT Structure_Id,PDB,Chain FROM Structure"
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in insertTemplates: %s\n%s" % (sql,f))

        for row in results:
            pdb_id = row[1]
            chain = row[2]
            if not (pdb_id,chain) in interacting_structure_dict:
                continue
            s_id = row[0]

            interacting_structure_ids[(pdb_id,chain)] = s_id

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
    for (pdb_id,chain)in ligand_map:
        iaps = ligand_map[(pdb_id,chain)]
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
                        (smiles,inchi) = pdb.getSI(pdb_id,name,res,chain,pdb_path)
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

    for (pdb_id,chain) in interacting_structure_dict:
        s_id = interacting_structure_ids[(pdb_id,chain)]
        (resolution,interaction_partner,oligos) = interacting_structure_dict[(pdb_id,chain)]
        for iap in interaction_partner:
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
                value_strs.append("('%s','%s','%s','%s')" % (str(ligand_id),str(s_id),str(chain),str(res)))

    t9 = time.time()

    if len(value_strs) > 0:
        sql = "INSERT IGNORE INTO RS_Ligand_Structure(Ligand,Structure,Chain,Residue) VALUES %s" % ','.join(value_strs)
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

    return interacting_structure_ids

#called by serializedPipeline
def insertAlignments(alignment_list,structure_id_map,stored_structures,db,cursor):
    #structure of alignment_list: [(gene_id,pdb_id,chain,seq_id,coverage,alignment_pir)]
    #structure of structure_id_map: {(pdb_id,chain):s_id}
    value_strs = []
    for (gene_id,pdb_id,chain,seq_id,coverage,alignment_pir) in alignment_list:
        if (pdb_id,chain) in structure_id_map:
            s_id = structure_id_map[(pdb_id,chain)]
        else:
            s_id = stored_structures[(pdb_id,chain)]
        value_strs.append("('%s','%s','%s','%s','%s')" % (str(gene_id),str(s_id),str(seq_id),str(coverage),alignment_pir))


    if len(value_strs) > 0:
        #For large amounts of templates, a single insert can destroy the connection to the database => Solution: divide large inserts into several smaller ones
        value_strs_list = []
        i = 0
        part_size = 500
        while i*part_size < len(value_strs):
            a = i*part_size
            b = min(((i+1)*part_size,len(value_strs)))
            i += 1
            sql = "INSERT IGNORE INTO Alignment(Gene,Structure,Sequence_Identity,Coverage,Alignment) VALUES %s" % ','.join(value_strs[a:b])
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in insertAlignments: %s" % (f))

#called by serializedPipeline
def insertResidues(annotation_map,db,cursor,structure_id_map,interacting_structure_ids,stored_structures):
    #stored_structures: contains the s-ids from all structures, which are already stored in the database
    #   [structure]: {(pdb-id,chain):s_id}
    #structure_id_map: contains the s-ids from all mapped structures

    values = []

    if len(annotation_map) == 0:
        return {},{}
    
    structure_ids = set()
    target_structure_ids = set()

    for (pdb,chain) in structure_id_map:
        structure_ids.add(structure_id_map[(pdb,chain)])
        if (pdb,chain) in stored_structures: #all residues belonging to stored structures must be filtered from the residue_dict; the dict gets updated later by 'stored_residue_dict' in serializedPipeline
            continue
        target_structure_ids.add(structure_id_map[(pdb,chain)])

    for (pdb,chain) in interacting_structure_ids:
        structure_ids.add(interacting_structure_ids[(pdb,chain)])

    for (pdb,chain) in annotation_map:

        if (pdb,chain) in stored_structures: #all residues belonging to stored structures must not inserted twice
            continue

        anno_map,residue_residue_map = annotation_map[(pdb,chain)]
        if (pdb,chain) in structure_id_map:
            s_id = structure_id_map[(pdb,chain)]
        else:
            s_id = interacting_structure_ids[(pdb,chain)]
        
        for res_id in anno_map:
            (one_letter,lig_dists,chain_dists,rsa,ssa,homomer_map,profile_str,centrality_scores,avg_b_factor,modres) = anno_map[res_id]

            if centrality_scores[0] == None:
                centrality_score_str = None
            else:
                centrality_score_str = ';'.join([str(x)[:40] for x in centrality_scores])

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

            if len(lig_dists) == 0:
                lig_dist_string = "-"
            if len(chain_dists) == 0:
                chain_dist_string = "-"
            res_id = res_id.replace("'","\\'")
            
            rsa = float(rsa)

            homo_strs = []
            for homo_chain in homomer_map:
                min_d = homomer_map[homo_chain]
                homo_strs.append('%s:%1.2f' % (homo_chain,min_d))
            homo_str = ','.join(homo_strs)

            values.append([s_id,res_id,one_letter,lig_dist_string,chain_dist_string,rsa,ssa,homo_str,profile_str,centrality_score_str,avg_b_factor,modres])

    columns = ['Structure','Number','Amino_Acid','Sub_Lig_Dist','Sub_Chain_Distances','Relative_Surface_Access','Secondary_Structure_Assignment',
                'Homomer_Distances','Interaction_Profile','Centralities','B_Factor','Modres']

    insert(db,cursor,'Residue',columns,values,value_threshold=500)

    residue_dict = {}
    structure_residue_map = {}
    if len(annotation_map) > 0:
        rows = ['Structure','Residue_Id','Number','Amino_Acid','Sub_Lig_Dist','Sub_Chain_Distances','Relative_Surface_Access','Secondary_Structure_Assignment','Homomer_Distances','Interaction_Profile','Centralities','Modres']
        table = 'Residue'

        results = binningSelect(structure_ids,rows,table,db,cursor)
        
        for row in results:
            s_id = row[0]
            r_id = row[1]
            
            if not s_id in structure_ids:
                continue
            res_nr = row[2]
            if not s_id in structure_residue_map:
                structure_residue_map[s_id] = {}
            structure_residue_map[s_id][res_nr] = r_id
            if not s_id in target_structure_ids:
                continue
            residue_dict[r_id] = (s_id,) + row[2:]

    residue_residue_values = []
    for (pdb,chain) in annotation_map:
        anno_map,residue_residue_map = annotation_map[(pdb,chain)]
        if (pdb,chain) in structure_id_map:
            s_id = structure_id_map[(pdb,chain)]
        else:
            s_id = interacting_structure_ids[(pdb,chain)]

        for res_nr_1 in residue_residue_map:
            r_id_1 = structure_residue_map[s_id][res_nr_1]
            for chain_2 in residue_residue_map[res_nr_1]:
                if (pdb,chain_2) in structure_id_map:
                    s_id_2 = structure_id_map[(pdb,chain_2)]
                else:
                    s_id_2 = interacting_structure_ids[(pdb,chain_2)]
                for res_nr_2 in residue_residue_map[res_nr_1][chain_2]:
                    distance = residue_residue_map[res_nr_1][chain_2][res_nr_2]
                    #print s_id_2,res_nr_2
                    r_id_2 = structure_residue_map[s_id_2][res_nr_2]
                    residue_residue_values.append("('%s','%s','%s','%s','%s')" % (str(r_id_1),str(s_id),str(r_id_2),str(s_id_2),str(distance)))

    size = len(residue_residue_values)
    if size > 0:
        bound = 500000
        if size > bound:
            m = size//bound
            rest = size%bound
            if rest == 0:
                part_size = size//m
            else:
                m += 1
                part_size = size//m
            parts = []
            for i in range(0,m-1):
                parts.append(residue_residue_values[(i*part_size):((i+1)*part_size)])
            parts.append(residue_residue_values[(m-1)*part_size:])
        else:
            parts = [residue_residue_values]
        for part in parts:
            sql = "INSERT IGNORE INTO RS_Residue_Residue(Residue_1,Structure_1,Residue_2,Structure_2,Distance) VALUES %s" % ','.join(part)
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Error in insertAnnotations: %s" % (f))

    return structure_residue_map,residue_dict


def processAlignmentData(alignment):
    lines = alignment.split("\n")
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
    return target_seq,template_seq

#called by serializedPipeline
def getAlignments(stored_genes,db,cursor,verbose=False):
    if verbose:
        t0 = time.time()

    rows = ['Gene','Structure','Sequence_Identity','Coverage,Alignment']
    table = 'Alignment'
    results = binningSelect(stored_genes,rows,table,db,cursor)

    if verbose:
        t1 = time.time()
        print("Time for part 1 in getAlignments: %s" % (str(t1-t0)))

    gene_structure_alignment_map = {}
    structure_ids = set()

    for row in results:
        #print row
        gene_id = row[0]
        if not gene_id in stored_genes:
            continue
        structure_id =row[1]
        seq_id = row[2]
        coverage = row[3]
        alignment = row[4]

        structure_ids.add(structure_id)

        target_seq,template_seq = processAlignmentData(alignment)

        if not gene_id in gene_structure_alignment_map:
            gene_structure_alignment_map[gene_id] = {}

        gene_structure_alignment_map[gene_id][structure_id] = (target_seq,template_seq,coverage,seq_id)

    if verbose:
        t2 = time.time()
        print("Time for part 2 in getAlignments: %s" % (str(t2-t1)))

    structure_map,id_structure_map = getStructure_map(structure_ids,db,cursor)

    if verbose:
        t3 = time.time()
        print("Time for part 3 in getAlignments: %s" % (str(t3-t2)))

    return gene_structure_alignment_map,structure_map,id_structure_map

def getStructure_map(structure_ids,db,cursor):
    structure_map = {}
    id_structure_map = {}
    if len(structure_ids) > 0:
        min_s = min(structure_ids)
        max_s = max(structure_ids)
        sql = "SELECT Structure_Id,PDB,Chain FROM Structure WHERE Structure_Id BETWEEN %s AND %s" % (min_s,max_s)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in getAlignments: %s,\n%s" % (sql,f))

        for row in results:
            s_id = row[0]
            if not s_id in structure_ids:
                continue
            pdb_id = row[1]
            chain = row[2]
            structure_map[(pdb_id,chain)] = s_id
            id_structure_map[s_id] = (pdb_id,chain)
    return structure_map,id_structure_map


def getSeqIds(s_ids,g_ids,db,cursor):
    seq_id_map = {}
    if len(s_ids) > 0 and len(g_ids) > 0:
        min_s = min(s_ids)
        max_s = max(s_ids)
        min_g = min(g_ids)
        max_g = max(g_ids)

        sql = "SELECT Structure,Gene,Sequence_Identity,Coverage FROM Alignment WHERE Structure BETWEEN %s AND %s AND Gene BETWEEN %s AND %s" % (min_s,max_s,min_g,max_g)
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in getAlignments: %s,\n%s" % (sql,f))

        for row in results:
            s_id = row[0]
            if not s_id in s_ids:
                continue
            g_id = row[1]
            if not g_id in g_ids:
                continue

            seq_id = row[2]
            cov = row[3]

            if not g_id in seq_id_map:
                seq_id_map[g_id] = {}
            seq_id_map[g_id][s_id] = (seq_id,cov)

    return seq_id_map

#called by serializedPipeline
def insertClassifications(db,cursor,residue_dict,structure_dict,gene_structure_map,m_r_map,config,complexMap):
    if config.verbose:
        t0 = time.time()
    mutation_dict = getMutationDict(m_r_map,db,cursor)

    if config.verbose:
        t1 = time.time()
        print('insertClassifications part 1: %s' % (t1-t0))

    print(len(residue_dict),len(structure_dict),len(gene_structure_map),len(m_r_map),len(mutation_dict))

    min_l_dict,min_m_dict,min_i_dict,min_c_dict,min_r_dict,min_d_dict,min_homo_dict,mutation_surface_dict,mutation_sec_dict,mutation_inter_dict,structure_classification_map,mutation_oligo_chain_map,centrality_dict = createStructureDicts(residue_dict,structure_dict,gene_structure_map,m_r_map,mutation_dict,ligand_filter=config.ligand_filter,n_of_processes=config.proc_n)

    if config.verbose:
        t2 = time.time()
        print('insertClassifications part 2: %s' % (t2-t1))

    class_dict,structure_classification_map,interaction_dict = createClassDict(min_l_dict,min_m_dict,min_i_dict,min_c_dict,min_r_dict,min_d_dict,min_homo_dict,mutation_surface_dict,mutation_sec_dict,centrality_dict,mutation_dict,structure_classification_map,complexMap,config.verbose)

    if config.verbose:
        t3 = time.time()
        print('insertClassifications part 3: %s' % (t3-t2))

    table = 'Mutation'
    rows = ['Mutation_Id','Amino_Acid_Change','Gene','Location','Class','RIN_Class','Simple_Class','RIN_Simple_Class','Interactions','Confidence','Secondary_Structure','Recommended_Structure','Max_Seq_Structure','Mapped_Structures','RIN_Profile','Modres_Score','Modres_Propensity']
    values = createClassValues(class_dict,mutation_sec_dict,mutation_dict,interaction_dict,mutation_inter_dict,n_of_processes=config.proc_n)

    if config.verbose:
        t4 = time.time()
        print('insertClassifications part 4: %s' % (t4-t3))

    update(db,cursor,table,rows,values)

    if config.verbose:
        t5 = time.time()
        print('insertClassifications part 5: %s' % (t5-t4))

def paraClassValues(inqueue,outqueue,lock,class_dict,mutation_sec_dict,mutation_dict,interaction_dict,mutation_inter_dict):
    with lock:
        inqueue.put(None)

    outs = []
    outs_package_size = 100000
    n = 0

    while True:
        with lock:
            m_ids = inqueue.get()
        if m_ids == None:
            break
        for m in m_ids:
            n += 1
            if m in mutation_sec_dict:
                mv_sec_ass = majority_vote(mutation_sec_dict[m])
            else:
                mv_sec_ass = None

            (Class,conf,weighted_sc,conf_sc,best_res,max_seq_res,amount_of_structures,
            weighted_c,conf_c,
            weighted_d,conf_d,
            weighted_r,conf_r,
            weighted_l,conf_l,
            weighted_m,conf_m,
            weighted_i,conf_i,
            weighted_h,conf_h,
            max_seq_id,
            weighted_raw,weighted_cent,weighted_norm,
            weighted_lig_degree,weighted_lig_score,
            weighted_metal_degree,weighted_metal_score,
            weighted_ion_degree,weighted_ion_score,
            weighted_prot_degree,weighted_prot_score,
            weighted_rna_degree,weighted_rna_score,
            weighted_dna_degree,weighted_dna_score,
            weighted_modres,modres_prop) = class_dict[m]
            simple_class = simplifyClass(Class,weighted_sc)

            if best_res != None:
                [r_id,qual,res_aa,res_nr,pdb_id,chain,resolution,cov,seq_id,rsa,min_lig,min_metal,min_ion,iacs] = best_res

                recommended_structure = '%s:%s %s:%s;%s;%s;%s' % (pdb_id,chain,res_nr,res_aa,str(seq_id),str(cov),str(resolution))
            else:
                resolution = '-'
                cov = '-'
                seq_id = '-'
                recommended_structure = '-'

            if max_seq_res != None:
                [max_seq_r_id,max_seq_qual,max_seq_res_aa,max_seq_res_nr,max_seq_pdb_id,max_seq_chain,max_seq_resolution,max_seq_cov,max_seq_seq_id,max_seq_rsa,max_min_lig,max_min_metal,max_min_ion,max_iacs] = max_seq_res

                max_seq_structure = '%s:%s %s:%s;%s;%s;%s' % (max_seq_pdb_id,max_seq_chain,max_seq_res_nr,max_seq_res_aa,str(max_seq_seq_id),str(max_seq_cov),str(max_seq_resolution))
            else:
                max_seq_resolution = '-'
                max_seq_cov = '-'
                max_seq_seq_id = '-'
                max_seq_structure = '-'

            (aachange,gene_id) = mutation_dict[m][:2]

            if m in interaction_dict:
                interaction_str = str(interaction_dict[m])
            else:
                interaction_str = None

            if m in mutation_inter_dict:
                profiles = []
                for profile_str,qual in mutation_inter_dict[m]:
                    profile = rin.Interaction_profile(profile_str=profile_str)
                    profiles.append((qual,profile))
                average_profile = rin.calculateAverageProfile(profiles)
                avg_profile_str = average_profile.encode()
                raw_rin_class,raw_rin_simple_class = average_profile.getClass()

                if raw_rin_class == 'No interaction':
                    rin_class = weighted_sc
                else:
                    rin_class = raw_rin_class

                if raw_rin_simple_class == 'No interaction':
                    rin_simple_class = weighted_sc
                else:
                    rin_simple_class = raw_rin_simple_class
            else:
                avg_profile_str = None
                rin_class = None
                rin_simple_class = None

            out = [m,aachange,gene_id,weighted_sc,Class,rin_class,simple_class,rin_simple_class,interaction_str,conf,mv_sec_ass,recommended_structure,max_seq_structure,amount_of_structures,avg_profile_str,weighted_modres,modres_prop]
            outs.append(out)
            if n == outs_package_size:
                with lock:
                    outqueue.put(outs)
                outs = []
                n = 0
    if outs != []:
        with lock:
            outqueue.put(outs)

def createClassValues(class_dict,mutation_sec_dict,mutation_dict,interaction_dict,mutation_inter_dict,n_of_processes=1):
    manager = multiprocessing.Manager()
    lock = manager.Lock()

    inqueue = manager.Queue()
    outqueue = manager.Queue()

    package_size = 100
    package = []

    i = 0
    n = 0

    for m in class_dict:
        i += 1
        package.append(m)

        if i == package_size:
            inqueue.put(package)
            package = []
            i = 0
            n += 1

    if package != []:
        inqueue.put(package)
        n += 1

    processes = {}
    for i in range(1,min(n+1,n_of_processes + 1)):
        p = multiprocessing.Process(target=paraClassValues, args=(inqueue,outqueue,lock,class_dict,mutation_sec_dict,mutation_dict,interaction_dict,mutation_inter_dict))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()

    with lock:
        outqueue.put(None)

    values = []
    while True:
        outs = outqueue.get()
        if outs == None:
            break
        for out in outs:
            values.append(out)
    return values

#called by serializedPipeline
def getStoredResidues(mappings,new_gene_stored_structure_mappings,stored_structures,gene_template_alignment_map,gene_mut_map_new,db,cursor,verbose = False):
    #gene_mut_map_new: {Uniprot_Id:(gene_id,{AAC_Base:mutation_id})}
    t0 = time.time()
    structure_res_map = {}

    for (m_id,s_id,res_nr,g_id) in mappings:
        if not s_id in structure_res_map:
            structure_res_map[s_id] = {}
        structure_res_map[s_id][res_nr] = None
        

    for (m_id,(pdb_id,chain),res_nr,gene_id) in new_gene_stored_structure_mappings:
        if not (pdb_id,chain) in stored_structures:
            continue
        s_id = stored_structures[(pdb_id,chain)]
        if not s_id in structure_res_map:
            structure_res_map[s_id] = {}
        structure_res_map[s_id][res_nr] = None

    for u_ac in gene_template_alignment_map:
        for (pdb_id,chain) in gene_template_alignment_map[u_ac]:
            if not (pdb_id,chain) in stored_structures:
                continue
            s_id = stored_structures[(pdb_id,chain)]
            if not s_id in structure_res_map:
                structure_res_map[s_id] = {}
            sub_infos = gene_template_alignment_map[u_ac][(pdb_id,chain)]
            for aacbase in sub_infos:
                m_id = gene_mut_map_new[u_ac][1][aacbase]
                res_nr = sub_infos[aacbase][0]
                if res_nr == None:
                    continue
                structure_res_map[s_id][res_nr] = None

    t1 = time.time()
    residue_dict = {}
    if len(structure_res_map) > 0:

        rows = ['Structure','Residue_Id','Number','Amino_Acid','Sub_Lig_Dist','Sub_Chain_Distances','Relative_Surface_Access','Secondary_Structure_Assignment','Homomer_Distances','Interaction_Profile','Centralities','Modres']
        table = 'Residue'
        results = binningSelect(structure_res_map.keys(),rows,table,db,cursor)
        
        for row in results:
            s_id = row[0]
            if not s_id in structure_res_map:
                continue
            res_nr = row[2]
            if not res_nr in structure_res_map[s_id]:
                continue
            r_id = row[1]
            structure_res_map[s_id][res_nr] = r_id
            residue_dict[r_id] = (s_id,) + row[2:]
    if verbose:
        t2 = time.time()
        print("Time for insertMapping Part 2: %s" % str(t2-t1))

    
    m_r_map = {}
    value_strs = []
    for (m_id,s_id,res_nr,g_id) in mappings:
        if res_nr == None:
            continue
        r_id = structure_res_map[s_id][res_nr]
        if r_id == None:
            continue
        if not m_id in m_r_map:
            m_r_map[m_id] = set()
        m_r_map[m_id].add(r_id)

    for (m_id,(pdb_id,chain),res_nr,g_id) in new_gene_stored_structure_mappings:
        if not (pdb_id,chain) in stored_structures:
            continue
        s_id = stored_structures[(pdb_id,chain)]
        if res_nr == None: #This could happen in the past, when the position is mapped to a gap in the alignment
            continue
        r_id = structure_res_map[s_id][res_nr]
        if r_id == None:
            continue
        if not m_id in m_r_map:
            m_r_map[m_id] = set()
        m_r_map[m_id].add(r_id)

    return m_r_map,residue_dict

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
        print("Error in class definition")
        print(lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc)
        print(sys.exc_info())
    if mut_class == "":
        print(lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc)
        
    return mut_class,t_id

def getClass(lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc):
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
                mut_class = "Surface"
            else:
                mut_class = "Core"
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
        print("Error in class definition")
        print(lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc)
        print(sys.exc_info())
    if mut_class == "":
        print(lig_sub_dist,chain_sub_dist,chain_type,rel_sur_acc)
        
    return mut_class

def majority_vote(secs):
    class_dict = {}
    for (sec,qual) in secs:
        if not sec in class_dict:
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
    #print aac
    try:
        if aac.count(',') < 1:
            aa1 = aac[0]
            aa2 = aac[-1]
            
            chemical_distance = chem_dist_matrix[aa1][aa2]
        else:
            aa1 = aac[0]
            aa2s = aac.split(",")
            aa2s[0] = aa2s[0][-1]
            chem_dists = []
            for aa2 in aa2s:
                chem_dists.append(chem_dist_matrix[aa1][aa2])
            chemical_distance = float(sum(chem_dists))/float(len(chem_dists))
    except:
        #print aac
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

    print("TESTPRINT: ",len(template_mutation_map))

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

def split(id_set,min_id,max_id):
    #print 'split ',min_id,max_id
    if len(id_set) == 0:
        return []
    if len(id_set) == 1:
        return [(id_set,min_id,max_id)]
    if len(id_set) == 2:
        #print id_set
        return [(set([min_id]),min_id,min_id),(set([max_id]),max_id,max_id)]
    split_id = (min_id+max_id)//2
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
    #print 'after split ',min_l,max_l,min_r,max_r
    return [(l,min_l,max_l),(r,min_r,max_r)]

def binning(id_set,tablesize=None,binsize=500000,density=0.5):
    if tablesize == 0:
        return set(),[]

    if tablesize != None:
        if tablesize > 0:
            factor = 1
            density = float(len(id_set))/(float(tablesize)*factor)
        else:
            density = 0.
    else:
        factor = 1

    bins = {}
    for i in id_set:
        bin_number = i//binsize
        #print bin_number
        if not bin_number in bins:
            bins[bin_number] = set()
        bins[bin_number].add(i)
    sorted_bins = []
    for bin_number in sorted(bins.keys()):
        id_set = bins[bin_number]
        sorted_bins.append((id_set,min(id_set),max(id_set)))

    '''
    size_distribution = {}
    for (id_set,min_id,max_id) in sorted_bins:
        s = len(id_set)
        if not s in size_distribution:
            size_distribution[s] = 1
        else:
            size_distribution[s] += 1
    for s in sorted(size_distribution.keys()):
        print 'For size ',s,' got number of bins: ',size_distribution[s]
    '''

    dense_enough = False
    while not dense_enough:
        dense_enough = True
        split_bins = []
        for (id_set,min_id,max_id) in sorted_bins:
            if len(id_set) < ((max_id - min_id)*density):
                for (iset,mi,ma) in split(id_set,min_id,max_id):
                    split_bins.append((iset,mi,ma))
                    dense_enough = False
            else:
                split_bins.append((id_set,min_id,max_id))

        sorted_bins = split_bins

    '''
    size_distribution = {}
    for (id_set,min_id,max_id) in sorted_bins:
        s = len(id_set)
        if not s in size_distribution:
            size_distribution[s] = 1
        else:
            size_distribution[s] += 1
    for s in sorted(size_distribution.keys()):
        print 'For size ',s,' got number of bins: ',size_distribution[s]
    '''

    fusion_done = False
    while not fusion_done:
        fusion_done = True
        fused_bins = []
        last_bin_fused = False
        for pos,(id_set,min_id,max_id) in enumerate(sorted_bins):
            if pos == 0 or last_bin_fused:
                last_bin_fused = False
                if pos == (len(sorted_bins) -1):
                    fused_bins.append((id_set,min_id,max_id))
                continue
            pre_id_set,pre_min_id,pre_max_id = sorted_bins[pos-1]
            if ((len(id_set)+len(pre_id_set))) > ((max_id - pre_min_id)*density*2*factor):

                fused_bins.append(((pre_id_set | id_set),pre_min_id,max_id))
                last_bin_fused = True
                fusion_done = False
            else:
                fused_bins.append((pre_id_set,pre_min_id,pre_max_id))
                last_bin_fused = False
                if pos == (len(sorted_bins) -1):
                    fused_bins.append((id_set,min_id,max_id))

        sorted_bins = fused_bins

    '''
    size_distribution = {}
    for (id_set,min_id,max_id) in sorted_bins:
        s = len(id_set)
        if not s in size_distribution:
            size_distribution[s] = 1
        else:
            size_distribution[s] += 1
    for s in sorted(size_distribution.keys()):
        print 'For size ',s,' got number of bins: ',size_distribution[s]
    '''

    print('sorted bins size: ',len(sorted_bins))
    singletons = set()
    non_singleton_bins = []
    for (id_set,min_id,max_id) in sorted_bins:
        if min_id == max_id:
            singletons.add(min_id)
        else:
            non_singleton_bins.append((id_set,min_id,max_id))
    return singletons,non_singleton_bins
    

#called by output
def minDistOut(outfolder,session_name,session_id,db,cursor,ligand_filter=None,intertable=False,templ_qual_separation=False,overwrite=False,processes=1,verbose=False):
    outfile = '%s/%s' % (outfolder,session_name)
    if verbose:
        t0 = time.time()

    if os.path.isfile("%s.classification.tsv" % outfile):
        print("A classification file already exists:\n%s.classification.tsv\n" % outfile)
        if overwrite:
            print("Overwrite flag is active, classification gets recalculated and old files will be overwritten.")
        else:
            print("The classification gets skipped, if you wish to reclassify, use the --overwrite flag.")
            files = os.listdir(outfolder)
            class_files = []
            for fname in files:
                if fname.count('.classification') > 0:
                    class_files.append('%s/%s' % (outfolder,fname))
            return class_files,[]

    if ligand_filter != None:
        f = open(ligand_filter,'r')
        lines = f.readlines()
        f.close()
        ligand_filter = set()
        for line in lines:
            ligand_filter.add(line.replace('\n','').replace(' ',''))
    global chem_dist_matrix

    if verbose:
        t1 = time.time()
        print("Time for mindistout part 1: ",t1-t0)

    table = 'RS_Mutation_Session'
    rows = ['Mutation','New_AA','Tag']
    eq_rows = {'Session':session_id}
    results = select(db,cursor,rows,table,equals_rows=eq_rows)

    if verbose:
        t2 = time.time()
        print("Time for mindistout part 2: ",t2-t1)

    new_aa_map = {}
    tag_map = {}
    for row in results:
        new_aa_map[row[0]] = row[1]
        tag_map[row[0]] = row[2]

    mutation_dict = getMutationDict(set(tag_map.keys()),db,cursor)
    gene_id_list = set()
    for m in mutation_dict:
        gene_id_list.add(mutation_dict[m][1])

    gene_score_dict = getGeneScoreDict(gene_id_list,session_id,db,cursor)
    #print gene_score_dict

    table = 'RS_Gene_Session'
    rows = ['Gene']
    eq_rows = {'Session':session_id}
    results = select(db,cursor,rows,table,equals_rows=eq_rows)

    gene_trunksize = 10000

    gene_ids = []
    gene_id_trunks = []
    n = 0
    for row in results:
        gene_ids.append(row[0])
        n+=1
        if n == gene_trunksize:
            gene_id_trunks.append(gene_ids[:])
            gene_ids = []
            n = 0
    if n != 0:
        gene_id_trunks.append(gene_ids[:])

    class_files = []
    interfiles = []

    option = ''
    class_file = "%s.classification%s.tsv" % (outfile,option)
    class_files.append(class_file)
    if overwrite and os.path.isfile(class_file):
        os.remove(class_file)

    if intertable:
        interfile = "%s.interaction_profiles%s.tsv" % (outfile,option)
        interfiles.append(interfile)
        if overwrite and os.path.isfile(interfile):
            os.remove(interfile)

    stat_file = "%s.statistics%s.tsv" % (outfile,option)
    if overwrite and os.path.isfile(stat_file):
        os.remove(stat_file)

    total_class_dict = {}

    trunk_number = 0
    for gene_ids in gene_id_trunks:

        binsize = 500000
        bin_threshs = []

        m_r_map = {}
        res_bins = {}
        res_ids = set()

        if verbose:
            t3 = time.time()
            print("Time for mindistout part 3: ",t3-t2)

        if not len(tag_map) == 0:
            #part 3.1
            if len(gene_ids) > 0:
                rows = ['Mutation','Residue']
                table = 'RS_Mutation_Residue'
                in_rows = {'Gene':gene_ids}
                results = select(db,cursor,rows,table,in_rows=in_rows)

                #between_rows = {'Gene':(min(gene_ids),max(gene_ids))}
                #results = select(db,cursor,rows,table,between_rows=between_rows)
            else:
                results = ()

            if verbose:
                t31 = time.time()
                print('Time for mindistout part 3.1: ',t31-t3)

            if not results == ():
                #minimal_bin_distance = len(results)/5000
                #print minimal_bin_distance
                for row in results:
                    m_id = row[0]
                    if not m_id in tag_map:
                        continue
                    r_id = row[1]
                    bin_number = r_id//binsize
                    #print bin_number
                    if not bin_number in res_bins:
                        res_bins[bin_number] = set()
                    res_bins[bin_number].add(r_id)

                    '''
                    if len(bin_threshs) == 0:
                        bin_threshs.append([r_id,r_id])
                    else:
                        for pos,[t0,t1] in enumerate(bin_threshs):
                            if r_id < t0:
                                if t0 - r_id > minimal_bin_distance:
                                    bin_threshs.insert(pos,[r_id,r_id])
                                    break
                                else:
                                    bin_threshs[pos][0] = r_id
                                    break
                            elif r_id < t1:
                                break
                            elif r_id - t1 > minimal_bin_distance:
                                if pos+1 == len(bin_threshs):
                                    bin_threshs.append([r_id,r_id])
                                    break
                                else:
                                    continue
                            elif pos+1 == len(bin_threshs):
                                bin_threshs[pos][1] = r_id
                                break
                            elif r_id + minimal_bin_distance > bin_threshs[pos+1][0]:
                                bin_threshs[pos][1] = bin_threshs[pos+1][1]
                                del bin_threshs[pos+1]
                                break
                            else:
                                bin_threshs[pos][1] = r_id
                                break
                    '''

                    res_ids.add(r_id)
                    if not m_id in m_r_map:
                        m_r_map[m_id] = set()
                    m_r_map[m_id].add(r_id)

            #print bin_threshs,len(bin_threshs)
        else:
            if verbose:
                t31 = time.time()

        sql = 'SELECT COUNT(*) FROM Residue'
        try:
            cursor.execute(sql)
            results = cursor.fetchall()
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in mindistout: %s,\n%s" % (sql,f))
            db.rollback()

        tablesize = results[0][0]

        singletons,res_bins = binning(res_ids,tablesize=tablesize)

        if verbose:
            t4 = time.time()
            print("Time for mindistout part 4: ",t4-t31)

        residue_dict = {}
        s_ids = set()

        #print(res_ids)

        #'''
        #part 5
        #print res_bins
        if not len(res_ids) == 0:
            '''
            max_r = max(res_ids)
            min_r = min(res_ids)

            min_max_tuples = []
            if max_r - min_r > binsize:
                for bin_number in res_bins:
                    min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize,res_bins[bin_number]))
            else:
                min_max_tuples = [(min_r,max_r,res_ids)]

            
            for (min_r,max_r) in min_max_tuples:
            #for [min_r,max_r] in bin_threshs:
            '''
            #for bin_number in res_bins:
            #    res_ids = res_bins[bin_number]
            #    max_r = max(res_ids)
            #    min_r = min(res_ids)
            if verbose:
                t50 = time.time()
            rows = ['Residue_Id','Structure','Number','Amino_Acid','Sub_Lig_Dist','Sub_Chain_Distances','Relative_Surface_Access','Secondary_Structure_Assignment','Homomer_Distances',
                                 'Interaction_Profile','Centralities']
            table = 'Residue'
            if len(singletons) > 0:
                in_rows = {'Residue_Id':singletons}
                results = select(db,cursor,rows,table,in_rows=in_rows)
                #if verbose:
                #    t51 = time.time()
                #    print 'Select residue singletons in total: ',len(singletons),'; Time: ', t51-t50

                if not results == ():
                    for row in results:
                        r_id = row[0]
                        if not r_id in res_ids:
                            continue
                        s_id = row[1]
                        s_ids.add(s_id)
                        residue_dict[r_id] = row[1:]

            for res_ids,min_r,max_r in res_bins:

                if verbose:
                    t51 = time.time()
                between_rows = {'Residue_Id':(min_r,max_r)}

                results = select(db,cursor,rows,table,between_rows=between_rows)

                #if verbose:
                #    t52 = time.time()
                #    print 'Select residue from ',min_r,' to ',max_r,'; in total: ',len(res_ids),'; Time: ', t52-t51

                if not results == ():
                    for row in results:
                        r_id = row[0]
                        if not r_id in res_ids:
                            continue
                        s_id = row[1]
                        s_ids.add(s_id)
                        residue_dict[r_id] = row[1:]

                #if verbose:
                #    t53 = time.time()
                #    print 'After Select residue Time: ', t53-t52
        #'''
        if verbose:
            t5 = time.time()
            print("Time for mindistout part 5: ",t5-t4)

        if verbose:
            t6 = time.time()
            print("Time for mindistout part 6: ",t6-t5)

        if verbose:
            t7 = time.time()
            print("Time for mindistout part 7: ",t7-t6)

        structure_dict = getStructureDict(s_ids,db,cursor)

        pdb_ids = set()
        for s_id in structure_dict:
            pdb_id = structure_dict[s_id][0]
            pdb_ids.add(pdb_id)

        complexMap = getComplexMap(db,cursor,pdb_ids=pdb_ids)
        if verbose:
            t8 = time.time()
            print("Time for mindistout part 8: ",t8-t7)
            t9 = time.time()

        gene_structure_map = {}

        if len(gene_id_list) > 0:
            rows = ['Gene','Structure','Sequence_Identity','Coverage']
            table = 'Alignment'
            between_rows = {'Gene':(min(gene_id_list),max(gene_id_list))}
            results = select(db,cursor,rows,table,between_rows=between_rows)

            t9 = time.time()
            #print "Time for mindistout part 9: ",t9-t8

            for row in results:
                g_id = row[0]
                if not g_id in gene_id_list:
                    continue
                s_id = row[1]
                seq_id = row[2]
                cov = row[3]
                if not g_id in gene_structure_map:
                    gene_structure_map[g_id] = {}
                gene_structure_map[g_id][s_id] = (seq_id,cov)

        if verbose:
            t10 = time.time()
            print("Time for mindistout part 10: ",t10-t9)

        #if templ_qual_separation: #TODO
        #    options_map = {"":total_results,".bioassembly":residue_dict_ba,".seqid90":residue_dict_seq,".seqid90.bioassembly":residue_dict_seq_ba}
        #else:
        #    options_map = {"":residue_dict}

        

        #limit amount of processes in order to limit memory usage
        processes = min((processes,10))

        print('Amount of processes going into createStructureDicts: ',processes)

        #for option in options_map:
            #residue_dict = options_map[option]

        
        min_l_dict,min_m_dict,min_i_dict,min_c_dict,min_r_dict,min_d_dict,min_homo_dict,mutation_surface_dict,mutation_sec_dict,mutation_inter_dict,structure_classification_map,mutation_oligo_chain_map,centrality_dict = createStructureDicts(residue_dict,structure_dict,gene_structure_map,m_r_map,mutation_dict,ligand_filter=ligand_filter,n_of_processes=processes)
        if verbose:
            t11 = time.time()
            print("Time for mindistout part 11: ",t11-t10)

        class_dict,structure_classification_map,interaction_dict = createClassDict(min_l_dict,min_m_dict,min_i_dict,min_c_dict,min_r_dict,min_d_dict,min_homo_dict,mutation_surface_dict,mutation_sec_dict,centrality_dict,mutation_dict,structure_classification_map,complexMap,verbose)
        total_class_dict.update(class_dict)
        if verbose:
            t12 = time.time()
            print("Time for mindistout part 12: ",t12-t11)

        if intertable:
            inter_dict = createInterDict(mutation_inter_dict)

        if trunk_number > 0:
            header=False
        else:
            header=True
        trunk_number+=1

        writeClassFile(class_file,mutation_surface_dict,mutation_sec_dict,mutation_dict,gene_score_dict,class_dict,tag_map,header=header)
        try:
            if intertable:
                writeInterFile(interfile,inter_dict,mutation_dict,gene_score_dict,new_aa_map,tag_map,class_dict,header=header)
        except:
            print('Error in writing the intertable')
        
    writeStatFile(stat_file,mutation_dict,total_class_dict,tag_map)

    return class_files,interfiles

#called by output.py
def classificationOutput(outfolder,session_name,session_id,db,cursor,ligand_filter=None,overwrite=False,verbose=False):
    outfile = '%s/%s' % (outfolder,session_name)
    if verbose:
        t0 = time.time()

    if os.path.isfile("%s.classification.tsv" % outfile):
        print("A classification file already exists:\n%s.classification.tsv\n" % outfile)
        if overwrite:
            print("Overwrite flag is active, classification gets recalculated and old files will be overwritten.")
        else:
            print("The classification gets skipped, if you wish to reclassify, use the --overwrite flag.")
            files = os.listdir(outfolder)
            class_files = []
            for fname in files:
                if fname.count('.classification') > 0:
                    class_files.append('%s/%s' % (outfolder,fname))
            return class_files,[]

    if ligand_filter != None:
        f = open(ligand_filter,'r')
        lines = f.readlines()
        f.close()
        ligand_filter = set()
        for line in lines:
            ligand_filter.add(line.replace('\n','').replace(' ',''))

    if verbose:
        t1 = time.time()
        print("Time for classificationOutput part 1: ",t1-t0)

    table = 'RS_Mutation_Session'
    rows = ['Mutation','New_AA','Tag']
    eq_rows = {'Session':session_id}
    results = select(db,cursor,rows,table,equals_rows=eq_rows)

    if verbose:
        t2 = time.time()
        print("Time for classificationOutput part 2: ",t2-t1)

    new_aa_map = {}
    tag_map = {}
    for row in results:
        new_aa_map[row[0]] = row[1]
        tag_map[row[0]] = row[2]

    mutation_dict = getMutationDict(set(tag_map.keys()),db,cursor)
    gene_id_list = set()
    for m in mutation_dict:
        gene_id_list.add(mutation_dict[m][1])

    gene_score_dict = getGeneScoreDict(gene_id_list,session_id,db,cursor)

    class_files = []

    class_file = "%s.classification.tsv" % (outfile)
    class_files.append(class_file)
    if overwrite and os.path.isfile(class_file):
        os.remove(class_file)

    stat_file = "%s.statistics.tsv" % (outfile)
    if overwrite and os.path.isfile(stat_file):
        os.remove(stat_file)

    if verbose:
        t3 = time.time()
        print("Time for classificationOutput part 3: ",t3-t2)

    rows = ['Mutation_Id','Amino_Acid_Change','Gene','Location','Class','Simple_Class','RIN_Class','RIN_Simple_Class','Interactions','Confidence','Secondary_Structure','Recommended_Structure','Max_Seq_Structure','Mapped_Structures','IUPRED_Glob']
    table = 'Mutation'
    results = binningSelect(mutation_dict.keys(),rows,table,db,cursor)

    if verbose:
        t4 = time.time()
        print("Time for classificationOutput part 4: ",t4-t3)

    startline = "Uniprot-Ac\tUniprot Id\tRefseq\tPDB-ID (Input)\tResidue-Id\tAmino Acid\tPosition\tSpecies\tTag\tWeighted Surface/Core\tClass\tSimple Class\tRIN Class\tRIN Simple Class\tIndividual Interactions\tConfidence Value\tSecondary Structure\tRecommended Structure\tSequence-ID\tCoverage\tResolution\tMax Seq Id Structure\tMax Sequence-ID\tMax Seq Id Coverage\tMax Seq Id Resolution\tAmount of mapped structures"

    lines = [startline]

    stat_dict = {}
    for row in results:
        m = row[0]
        aac = row[1]
        gene_id = row[2]
        weighted_sc = row[3]
        Class =row[4]
        simple_class = row[5]
        rin_class = row[6]
        rin_simple_class = row[7]
        interaction_str = row[8]
        conf = row[9]
        mv_sec_ass = row[10]
        recommended_structure_str = row[11]
        max_seq_structure_str = row[12]
        amount_of_structures = row[13]
        disorder_state = row[14]

        input_res_id = mutation_dict[m][4]
        if input_res_id == None:
            input_res_id = ''

        if recommended_structure_str != None and recommended_structure_str != '-':
            recommended_structure,seq_id,cov,resolution = recommended_structure_str.split(';')
        else:
            resolution = '-'
            cov = '-'
            seq_id = '-'
            recommended_structure = '-'

        if disorder_state != None:
            if Class == None and (disorder_state == 'disorder' or disorder_state[0] == 'D'):
                Class = 'Disorder'
                simple_class = 'Disorder'
                rin_class = 'Disorder'
                rin_simple_class = 'Disorder'

        if max_seq_structure_str != None and max_seq_structure_str != '-':
            max_seq_structure,max_seq_seq_id,max_seq_cov,max_seq_resolution = max_seq_structure_str.split(';')
            stat_dict[m] = (Class,float(max_seq_seq_id))
        else:
            max_seq_resolution = '-'
            max_seq_cov = '-'
            max_seq_seq_id = '-'
            max_seq_structure = '-'

            stat_dict[m] = (Class,0.)

        (u_ac,gpan,u_id,error_code,error,species,gene_score) = gene_score_dict[gene_id]

        input_pdb_id = ''
        if len(u_ac) == 6 and u_ac[4] == ':':
            input_pdb_id = u_ac
            u_ac = ''

        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (u_ac,u_id,gpan,input_pdb_id,input_res_id,aac[0],aac[1:],species,tag_map[m],weighted_sc,Class,simple_class,rin_class,rin_simple_class,interaction_str,str(conf),mv_sec_ass,recommended_structure,seq_id,cov,resolution,max_seq_structure,max_seq_seq_id,max_seq_cov,max_seq_resolution,str(amount_of_structures)))
    f = open(class_file,'a')
    f.write("\n".join(lines))
    f.close()

    if verbose:
        t5 = time.time()
        print("Time for classificationOutput part 5: ",t5-t4)

    writeStatFile(stat_file,mutation_dict,{},tag_map,stat_dict=stat_dict)

    if verbose:
        t6 = time.time()
        print("Time for classificationOutput part 6: ",t6-t5)

    return class_files

#method for comparing/developing the new classifiction
def diffSurfs(mutation_surface_dict,g_u_dict,mutation_dict,outfile):
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

def processClassData(sld,scd,chains,ligand_filter=None):
    min_ld = None
    min_md = None
    min_id = None
    min_cd = None
    min_dd = None
    min_rd = None

    slds = sld.split(",")
    ldists = {}
    mdists = {}
    idists = {}
    for ld in slds:
        ldi = ld.split(":")
        if len(ldi) > 1:
            lig_name,res,chain = ldi[0].split('_')
            if ligand_filter == None:
                if lig_name in metal_atoms:
                    mdists[float(ldi[1])] = lig_name,res,chain
                elif lig_name in ion_atoms:
                    idists[float(ldi[1])] = lig_name,res,chain
                else:
                    ldists[float(ldi[1])] = lig_name,res,chain
            elif not lig_name in ligand_filter:
                if lig_name in metal_atoms:
                    mdists[float(ldi[1])] = lig_name,res,chain
                elif lig_name in ion_atoms:
                    idists[float(ldi[1])] = lig_name,res,chain
                else:
                    ldists[float(ldi[1])] = lig_name,res,chain
    min_lig = None
    min_metal = None
    min_ion = None

    if len(ldists) > 0:
        min_ld = min(ldists.keys())
        min_lig = ldists[min_ld]

    if len(mdists) > 0:
        min_md = min(mdists.keys())
        min_metal = mdists[min_md]

    if len(idists) > 0:
        min_id = min(idists.keys())
        min_ion = idists[min_id]

    chain_dict = {}
    for chain_tuple in chains.split(','):
        [ichain,chaintype] = chain_tuple.split(':')
        chain_dict[ichain] = chaintype

    scds = scd.split(",")
    cdists = {}
    ddists = {}
    rdists = {}

    iacs = {}

    for sd in scds:
        sdi = sd.split(":")
        if len(sdi) > 1:
            ichain = sdi[0][0]
            mc_d = float(sdi[1])
            
            if not ichain in chain_dict:
                chaintype = 'Protein'
            else:
                chaintype = chain_dict[ichain]

            if chaintype == "Protein" or chaintype == 'Peptide':
                cdists[mc_d] = ichain
            elif chaintype == "RNA":
                rdists[mc_d] = ichain
            elif chaintype == "DNA":
                ddists[mc_d] = ichain

    if len(cdists) > 0:
        min_cd = min(cdists.keys())
        iacs['Protein'] = cdists[min_cd]
    if len(rdists) > 0:
        min_rd = min(rdists.keys())
        iacs['RNA'] = rdists[min_rd]
    if len(ddists) > 0:
        min_dd = min(ddists.keys())
        iacs['DNA'] = ddists[min_dd]

    return min_ld,min_md,min_id,min_cd,min_rd,min_dd,min_lig,min_metal,min_ion,iacs

def paraStructureDicts(inqueue,outqueue,lock,residue_dict,structure_dict,gene_structure_map,m_r_map,mutation_dict,ligand_filter):
    residue_calc = {}

    with lock:
        inqueue.put(None)

    outs = []
    outs_package_size = 100000
    n = 0

    while True:
        with lock:
            m_ids = inqueue.get()
        if m_ids == None:
            break
        for m_id in m_ids:
            if not m_id in mutation_dict:
                continue
            (aachange,gene_id,iupred_score,glob,input_res_id) = mutation_dict[m_id]

            for r_id in m_r_map[m_id]:
                n += 1

                row = residue_dict[r_id]
                s_id = row[0]
                residue = row[1]
                res_aa = row[2]
                sld = row[3]
                scd = row[4]
                homomer_dists = row[7]
                centralities = row[9]
                modres= row[10]
                if centralities != None and centralities != 'NULL':
                    Raw_Centrality,Centrality,Norm_Centrality,Cleaned_Raw_Centrality,Cleaned_Centrality,Cleaned_Norm_Centrality = [float(x) for x in centralities.split(';')]
                else:
                    Raw_Centrality = None
                    Centrality = None
                    Norm_Centrality = None
                    Cleaned_Raw_Centrality = None
                    Cleaned_Centrality = None
                    Cleaned_Norm_Centrality = None
                structure = structure_dict[s_id]
                
                pdb_id = structure[0]
                chain = structure[1]
                resolution = structure[2]
                ligands = structure[3]
                oligos = structure[4]
                chains = structure[5]


                if not s_id in gene_structure_map[gene_id]:
                    print('Warning in paraStructureDicts: ',gene_id,s_id,r_id,m_id)
                    #This seems to happen, when there is a mutation-residue mapping in the db, but not the corresponding alignment between structure and protein
                    continue
                (seq_id,cov) = gene_structure_map[gene_id][s_id]

                qual = templateFiltering.qualityScore(resolution,cov,seq_id)

                if not r_id in residue_calc:
                    try:
                        rel_sur_acc = float(row[5])
                        if rel_sur_acc > 1.0:
                            rel_sur_acc = 0.01*rel_sur_acc
                    except:
                        rel_sur_acc = None
                    try:
                        sec_str_ass = str(row[6])
                    except:
                        sec_str_ass = None

                    profile_str = row[8]
                    
                    min_ld,min_md,min_id,min_cd,min_rd,min_dd,min_lig,min_metal,min_ion,iacs = processClassData(sld,scd,chains,ligand_filter=ligand_filter)

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
                    if min_id != None:
                        minimal_distances.append(min_id)

                    if len(minimal_distances) == 0:
                        min_minimal_distances = 2.0
                    else:
                        min_minimal_distances = min(minimal_distances)

                    if homomer_dists == None or homomer_dists == '':
                        homomer_dist = None
                    else:
                        h_dists = []
                        for h_dist_str in homomer_dists.split(','):
                            h_dist = float(h_dist_str.split(':')[1].split('(')[0])
                            h_dists.append(h_dist)
                        homomer_dist = min(h_dists)

                    residue_calc[r_id] = (min_minimal_distances,rel_sur_acc,sec_str_ass,homomer_dist,min_ld,min_md,min_id,min_cd,min_rd,min_dd,
                                        profile_str,Raw_Centrality,Centrality,Norm_Centrality,min_lig,min_metal,min_ion,iacs)
                    

                else:
                    (min_minimal_distances,rel_sur_acc,sec_str_ass,homomer_dist,min_ld,min_md,min_id,min_cd,min_rd,min_dd,
                                        profile_str,Raw_Centrality,Centrality,Norm_Centrality,min_lig,min_metal,min_ion,iacs) = residue_calc[r_id]

                if rel_sur_acc > surface_threshold:
                    sc = "Surface"
                else:
                    sc = "Core"
                Class,conf = getWeightedClass(sc,1.0,min_cd,1.0,min_dd,1.0,min_rd,1.0,min_ld,1.0,min_md,1.0,min_id,1.0)
                simpleClass =  simplifyClass(Class,sc)

                out = (m_id,r_id,Class,simpleClass,qual,res_aa,residue,pdb_id,chain,resolution,cov,seq_id,
                                profile_str,
                                Raw_Centrality,Centrality,Norm_Centrality,
                                oligos,chains,
                                rel_sur_acc,
                                sec_str_ass,
                                homomer_dist,
                                min_ld,
                                min_md,
                                min_id,
                                min_cd,
                                min_rd,
                                min_dd,
                                min_minimal_distances,
                                min_lig,min_metal,min_ion,iacs,modres
                        )
                outs.append(out)
                if n == outs_package_size:
                    with lock:
                        outqueue.put(outs)
                    outs = []
                    n = 0
    if outs != []:
        with lock:
            outqueue.put(outs)
    

def createStructureDicts(residue_dict,structure_dict,gene_structure_map,m_r_map,mutation_dict,ligand_filter=None,n_of_processes=1):
    #print('MP start Method: ',multiprocessing.get_start_method())
    manager = multiprocessing.Manager()
    lock = manager.Lock()

    min_l_dict = {}
    min_m_dict = {}
    min_i_dict = {}
    min_c_dict = {}
    min_r_dict = {}
    min_d_dict = {}
    min_homo_dict = {}

    mutation_surface_dict = {}

    mutation_sec_dict = {}

    mutation_inter_dict = {}

    mutation_oligo_chain_map = {}

    centrality_dict = {}

    residue_calc = {}

    #needed in order to pick the recommended structure
    structure_classification_map = {}

    inqueue = manager.Queue()
    outqueue = manager.Queue()

    #print 'Amount of positions going into createStructureDicts: ',len(m_r_map)

    package_size = 100
    package = []

    i = 0
    n = 0

    for m_id in m_r_map:
        i += 1
        
        structure_classification_map[m_id] = {}
        mutation_inter_dict[m_id] = []
        mutation_oligo_chain_map[m_id] = [0,0]
        mutation_surface_dict[m_id] = []
        mutation_sec_dict[m_id] = []
        centrality_dict[m_id] = []
        min_l_dict[m_id] = []
        min_m_dict[m_id] = []
        min_i_dict[m_id] = []
        min_c_dict[m_id] = []
        min_r_dict[m_id] = []
        min_d_dict[m_id] = []
        min_homo_dict[m_id] = []
        package.append(m_id)
        if i == package_size:
            inqueue.put(package)
            package = []
            i = 0
            n += 1

    if package != []:
        inqueue.put(package)
        n += 1

    print('Amount of packages: ',n)

    processes = {}
    for i in range(1,min(n+1,n_of_processes + 1)):
        p = multiprocessing.Process(target=paraStructureDicts, args=(inqueue,outqueue,lock,residue_dict,structure_dict,gene_structure_map,m_r_map,mutation_dict,ligand_filter))
        processes[i] = p
        p.start()
    for i in processes:
        processes[i].join()
    with lock:
        outqueue.put(None)

    while True:
        outs = outqueue.get()
        if outs == None:
            break
        for out in outs:
            (m_id,r_id,Class,simpleClass,qual,res_aa,residue,pdb_id,chain,resolution,cov,seq_id,
                                profile_str,
                                Raw_Centrality,Centrality,Norm_Centrality,
                                oligos,chains,
                                rel_sur_acc,
                                sec_str_ass,
                                homomer_dist,
                                min_ld,
                                min_md,
                                min_id,
                                min_cd,
                                min_rd,
                                min_dd,
                                min_minimal_distances,
                                min_lig,min_metal,min_ion,iacs,modres
                        ) = out

            structure_classification_map[m_id][r_id] = (Class,simpleClass,qual,res_aa,residue,pdb_id,chain,resolution,cov,seq_id,rel_sur_acc,min_lig,min_metal,min_ion,iacs,modres)
            if profile_str != 'NULL':
                mutation_inter_dict[m_id].append((profile_str,qual))
            if Centrality != None and Centrality != 'NULL':
                centrality_dict[m_id].append((Raw_Centrality,Centrality,Norm_Centrality,qual)) 

            if len(oligos) > mutation_oligo_chain_map[m_id][0]:
                mutation_oligo_chain_map[m_id][0] = len(oligos)
            if len(chains.split(',')) > mutation_oligo_chain_map[m_id][1]:
                mutation_oligo_chain_map[m_id][1] = len(chains.split(','))

            if not min_minimal_distances < 1.2:
                #print rel_sur_acc
                if rel_sur_acc != None and cov > 0.0:
                    mutation_surface_dict[m_id].append((rel_sur_acc,qual,cov))

                mutation_sec_dict[m_id].append((sec_str_ass,qual))

                if not min_ld == None:
                    min_l_dict[m_id].append((min_ld,qual))
                if not min_md == None:
                    min_m_dict[m_id].append((min_md,qual))
                if not min_id == None:
                    min_i_dict[m_id].append((min_id,qual))
                if not min_cd == None:
                    min_c_dict[m_id].append((min_cd,qual))
                if not min_rd == None:
                    min_r_dict[m_id].append((min_rd,qual))
                if not min_dd == None:
                    min_d_dict[m_id].append((min_dd,qual))

                if not homomer_dist == None:
                    min_homo_dict[m_id].append((homomer_dist,qual))


    #print mutation_surface_dict

    return min_l_dict,min_m_dict,min_i_dict,min_c_dict,min_r_dict,min_d_dict,min_homo_dict,mutation_surface_dict,mutation_sec_dict,mutation_inter_dict,structure_classification_map,mutation_oligo_chain_map,centrality_dict

def createInterDict(mutation_inter_dict,chain_type='sc'):

    inter_dict = {}

    for m_id in mutation_inter_dict:
        profiles = mutation_inter_dict[m_id]
        total_qual = 0.0
        ion_qual = 0.0
        metal_qual = 0.0
        lig_qual = 0.0
        chain_qual = 0.0
        average_profile = [0.0]*14
        for profile_str,qual in profiles:
            if None in (profile_str,qual):
                continue
            profile = rin.Interaction_profile(profile_str=profile_str)
            total_qual += qual
            
            Ion_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type,'ion')
            Ion_Interaction_Score = profile.getChainSpecificCombiScore(chain_type,'ion')
            Metal_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type,'metal')
            Metal_Interaction_Score = profile.getChainSpecificCombiScore(chain_type,'metal')
            Ligand_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type,'ligand')
            Ligand_Interaction_Score = profile.getChainSpecificCombiScore(chain_type,'ligand')
            Chain_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type,'interchain')
            Chain_Interaction_Score = profile.getChainSpecificCombiScore(chain_type,'interchain')
            Short_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type,'neighbor')
            Short_Interaction_Score = profile.getChainSpecificCombiScore(chain_type,'neighbor')
            Medium_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type,'short')
            Medium_Interaction_Score = profile.getChainSpecificCombiScore(chain_type,'short')
            Long_Interaction_Degree = profile.getChainSpecificCombiDegree(chain_type,'long')
            Long_Interaction_Score = profile.getChainSpecificCombiScore(chain_type,'long')

            if Ligand_Interaction_Degree > 0:
                lig_qual+= qual
            if Metal_Interaction_Degree > 0:
                metal_qual += qual
            if Ion_Interaction_Degree > 0:
                ion_qual += qual
            if Chain_Interaction_Degree > 0:
                chain_qual += qual

            average_profile[0] += Ion_Interaction_Degree*qual
            average_profile[1] += Ion_Interaction_Score*qual
            average_profile[2] += Metal_Interaction_Degree*qual
            average_profile[3] += Metal_Interaction_Score*qual
            average_profile[4] += Ligand_Interaction_Degree*qual
            average_profile[5] += Ligand_Interaction_Score*qual
            average_profile[6] += Chain_Interaction_Degree*qual
            average_profile[7] += Chain_Interaction_Score*qual
            average_profile[8] += Short_Interaction_Degree*qual
            average_profile[9] += Short_Interaction_Score*qual
            average_profile[10] += Medium_Interaction_Degree*qual
            average_profile[11] += Medium_Interaction_Score*qual
            average_profile[12] += Long_Interaction_Degree*qual
            average_profile[13] += Long_Interaction_Score*qual

        if total_qual > 0.0:
            average_profile[8] = average_profile[8]/total_qual
            average_profile[9] = average_profile[9]/total_qual
            average_profile[10] = average_profile[10]/total_qual
            average_profile[11] = average_profile[11]/total_qual
            average_profile[12] = average_profile[12]/total_qual
            average_profile[13] = average_profile[13]/total_qual

        if ion_qual > 0.0:
            average_profile[0] = average_profile[0]/ion_qual
            average_profile[1] = average_profile[1]/ion_qual
        if metal_qual > 0.0:
            average_profile[2] = average_profile[2]/metal_qual
            average_profile[3] = average_profile[3]/metal_qual
        if lig_qual > 0.0:
            average_profile[4] = average_profile[4]/lig_qual
            average_profile[5] = average_profile[5]/lig_qual
        if chain_qual > 0.0:
            average_profile[6] = average_profile[6]/chain_qual
            average_profile[7] = average_profile[7]/chain_qual
        inter_dict[m_id] = average_profile

    return inter_dict

def median(l):
    #print l
    n = len(l)
    l = sorted(l)
    if n % 2 == 0:
        med = (l[(n//2)-1]+l[n//2])/2.0
    else:
        med = l[(n-1)//2]
    return med

def getWeightedClass(weighted_sc,conf_sc,weighted_c,conf_c,weighted_d,conf_d,weighted_r,conf_r,weighted_l,conf_l,weighted_m,conf_m,weighted_i,conf_i):
    dt = short_distance_threshold
    edt = long_distance_threshold

    contacts = []
    confs = []

    # Determine the near macromolecule with the highest conf.
    # If there is no near macromolecule, also consider far ones.
    # In case of multiple maximum confs use the first one added to macros.
    macros = []
    macro_confs = []
    macros_far = []
    macro_confs_far = []

    if weighted_c and weighted_c < edt:
        if weighted_c < dt:
            macros.append('Protein')
            macro_confs.append(conf_c)
        macros_far.append('Protein far')
        macro_confs_far.append(conf_c)

    if weighted_d and weighted_d < edt:
        if weighted_d < dt:
            macros.append('DNA')
            macro_confs.append(conf_d)
        macros_far.append('DNA far')
        macro_confs_far.append(conf_d)

    if weighted_r and weighted_r < edt:
        if weighted_r < dt:
            macros.append('RNA')
            macro_confs.append(conf_r)
        macros_far.append('RNA far')
        macro_confs_far.append(conf_r)

    contact_macro = None
    conf_macro = None

    if len(macros) > 0:
        best_macro_idx = np.min(np.argmax(macro_confs))
        contacts.append(macros[best_macro_idx])
        confs.append(macro_confs[best_macro_idx])
    elif len(macros_far) > 0:
        best_macro_idx = np.min(np.argmax(macro_confs_far))
        contacts.append(macros_far[best_macro_idx])
        confs.append(macro_confs_far[best_macro_idx])

    # Additional interactions
    if weighted_l and weighted_l < edt:
        if weighted_l < dt:
            contacts.append('Ligand')
        else:
            contacts.append('Ligand far')
        confs.append(conf_l)

    if weighted_m and weighted_m < edt:
        if weighted_m < dt:
            contacts.append('Metal')
        else:
            contacts.append('Metal far')
        confs.append(conf_m)

    if weighted_i and weighted_i < edt:
        if weighted_i < dt:
            contacts.append('Ion')
        else:
            contacts.append('Ion far')
        confs.append(conf_i)

    if len(contacts) == 0:
        return weighted_sc,conf_sc

    clas = " and ".join(contacts)
    if len(contacts) == 4:
        clas = "Quadruple Interaction: " + clas
    elif len(contacts) == 3:
        clas = "Triple Interaction: " + clas
    elif len(contacts) == 2:
        clas = "Double Interaction: " + clas
    elif len(contacts) == 1:
        clas = clas + " Interaction"

    conf = sum(confs) / len(confs)

    return clas,conf

def createClassDict(min_l_dict,min_m_dict,min_i_dict,min_c_dict,min_r_dict,min_d_dict,min_homo_dict,mutation_surface_dict,mutation_sec_dict,centrality_dict,mutation_dict,structure_classification_map,complexMap,verbose):
    distance_threshold = long_distance_threshold
    class_dict = {}

    t0 = 0.0
    t1= 0.0
    t2 = 0.0
    t3 = 0.0
    t4 = 0.0
    interaction_dict = {}
    for m in mutation_surface_dict:
        if len(mutation_surface_dict[m]) == 0:
            continue
        t0 += time.time()
        aac,g,iupred,glob,input_res_id = mutation_dict[m]

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

        t1 += time.time()

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
                weighted_c = None
                conf_c = 0.0
        else:
            weighted_c = None
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
                weighted_d = None
                conf_d = 0.0
        else:
            weighted_d = None
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
                weighted_r = None
                conf_r = 0.0
        else:
            weighted_r = None
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
                weighted_l = None
                conf_l = 0.0
        else:
            weighted_l = None
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
                weighted_m = None
                conf_m = 0.0
        else:
            weighted_m = None
            conf_m = 0.0

        if m in min_i_dict:
            nom = 0.0
            denom = 0.0
            n = 0.0
            qs = []
            for (min_i,qual) in min_i_dict[m]:
                if min_i < distance_threshold:
                    nom += min_i*qual
                    denom += qual
                    n += 1.0
                    qs.append(qual)
            if denom > 0.0:
                weighted_i = nom/denom
                conf_i = (1.0-1.0/(n+1.0))*(max(qs)+median(qs))/2
            else:
                weighted_i = None
                conf_i = 0.0
        else:
            weighted_i = None
            conf_i = 0.0
        
        if m in min_homo_dict:
            nom = 0.0
            denom = 0.0
            n = 0.0
            qs = []
            for (min_h,qual) in min_homo_dict[m]:
                nom += min_h*qual
                denom += qual
                n += 1.0
                qs.append(qual)
            if denom > 0.0:
                weighted_h = nom/denom
                conf_h = (1.0-1.0/(n+1.0))*(max(qs)+median(qs))/2
            else:
                weighted_h = None
                conf_h = 0.0
        else:
            weighted_h = None
            conf_h = 0.0

        if m in centrality_dict:
            raw_nom = 0.0
            nom = 0.0
            norm_nom = 0.0
            denom = 0.0
            n = 0.0
            for (raw,cent,norm,qual) in centrality_dict[m]:
                raw_nom += raw*qual
                nom += cent*qual
                norm_nom += norm*qual
                denom += qual
                n += 1.0
            if denom > 0.0:
                weighted_raw = raw_nom/denom
                weighted_cent = nom/denom
                weighted_norm = norm_nom/denom
            else:
                weighted_raw = 0.0
                weighted_cent = 0.0
                weighted_norm = 0.0
        else:
            weighted_raw = None
            weighted_cent = None
            weighted_norm = None

        t2 += time.time()


        mut_class,conf = getWeightedClass(weighted_sc,conf_sc,weighted_c,conf_c,weighted_d,conf_d,weighted_r,conf_r,weighted_l,conf_l,weighted_m,conf_m,weighted_i,conf_i)
        
        simple_mut_class = simplifyClass(mut_class,weighted_sc)

        fe_mut_class = excludeFarClasses(mut_class,weighted_sc) 

        t3 += time.time()

        

        #if m in structure_classification_map:
        best_res_iaa = None
        best_res = None
        best_res_iaa_isc = None
        best_res_isc = None
        best_res_icc = None
        max_seq_res = None
        amount_of_structures = len(structure_classification_map[m])
        max_seq_id = 0.

        total_lig_qual = 0.
        total_lig_degree = 0
        total_lig_score = 0.

        total_metal_qual = 0.
        total_metal_degree = 0
        total_metal_score = 0.

        total_ion_qual = 0.
        total_ion_degree = 0
        total_ion_score = 0.

        total_prot_qual = 0.
        total_prot_degree = 0
        total_prot_score = 0.

        total_rna_qual = 0.
        total_rna_degree = 0
        total_rna_score = 0.

        total_dna_qual = 0.
        total_dna_degree = 0
        total_dna_score = 0.

        interaction_dict[m] = {}
        weighted_modres = 0.
        modres_prop = 0.
        for r_id in structure_classification_map[m]:
            residue_class,residue_simple_class,qual,res_aa,res_nr,pdb_id,chain,resolution,cov,seq_id,rsa,min_lig,min_metal,min_ion,iacs,modres = structure_classification_map[m][r_id]

            if modres:
                weighted_modres += qual
                modres_prop += 1.
            else:
                weighted_modres -= qual

            ligand_profile,metal_profile,ion_profile,cc_profile = complexMap[pdb_id]

            fe_res_class = excludeFarClasses(residue_class,residue_simple_class) 

            if fe_res_class.count('Interaction') > 0:
                if not fe_res_class in interaction_dict[m]:
                    interaction_dict[m][fe_res_class] = (qual,pdb_id,chain,res_nr)
                elif qual > interaction_dict[m][fe_res_class][0]:
                    interaction_dict[m][fe_res_class] = (qual,pdb_id,chain,res_nr)

            res_infos = [r_id,qual,res_aa,res_nr,pdb_id,chain,resolution,cov,seq_id,rsa,min_lig,min_metal,min_ion,iacs]
            if seq_id > max_seq_id:
                max_seq_id = seq_id

            if max_seq_res == None:
                max_seq_res = res_infos
            elif seq_id > max_seq_res[8]:
                max_seq_res = res_infos

            #if (mut_class == residue_class):# and simple_mut_class == residue_simple_class):
            if fe_res_class == fe_mut_class:
                if best_res == None:
                    best_res = res_infos
                else:
                    if qual > best_res[1]:
                        best_res = res_infos
                if res_aa == aac[0]:
                    if best_res_iaa == None:
                        best_res_iaa = res_infos
                    else:
                        if qual > best_res_iaa[1]:
                            best_res_iaa = res_infos
                elif best_res_icc == None:
                    best_res_icc = res_infos
                else:
                    if qual > best_res_icc[1]:
                        best_res_icc = res_infos
                    
            elif simple_mut_class == residue_simple_class:
                if best_res_isc == None:
                    best_res_isc = res_infos
                else:
                    if qual > best_res_isc[1]:
                        best_res_isc = res_infos
                if res_aa == aac[0]:
                    if best_res_iaa_isc == None:
                        best_res_iaa_isc = res_infos
                    else:
                        if qual > best_res_iaa_isc[1]:
                            best_res_iaa_isc = res_infos

            if min_lig != None:
                lig_name,lig_res,lig_chain = min_lig

                if (lig_chain,lig_res) in ligand_profile:
                    lig_degree,lig_score = ligand_profile[(lig_chain,lig_res)]
                    total_lig_qual += qual
                    total_lig_degree += lig_degree*qual
                    total_lig_score += lig_score*qual

            if min_metal != None:
                metal_name,metal_res,metal_chain = min_metal
                if (metal_chain,metal_res) in metal_profile:
                    metal_degree,metal_score = metal_profile[(metal_chain,metal_res)]
                    total_metal_qual += qual
                    total_metal_degree += metal_degree*qual
                    total_metal_score += metal_score*qual

            if min_ion != None:
                ion_name,ion_res,ion_chain = min_ion
                if (ion_chain,ion_res) in ion_profile:
                    ion_degree,ion_score = ion_profile[(ion_chain,ion_res)]
                    total_ion_qual += qual
                    total_ion_degree += ion_degree*qual
                    total_ion_score += ion_score*qual

            if 'Protein' in iacs:
                ichain = iacs['Protein']
                if (chain,ichain) in cc_profile:
                    prot_degree,prot_score = cc_profile[(chain,ichain)]
                    total_prot_qual += qual
                    total_prot_degree += prot_degree*qual
                    total_prot_score += prot_score*qual


            if 'RNA' in iacs:
                ichain = iacs['RNA']
                if (chain,ichain) in cc_profile:
                    rna_degree,rna_score = cc_profile[(chain,ichain)]
                    total_rna_qual += qual
                    total_rna_degree += rna_degree*qual
                    total_rna_score += rna_score*qual

            if 'DNA' in iacs:
                ichain = iacs['DNA']
                if (chain,ichain) in cc_profile:
                    dna_degree,dna_score = cc_profile[(chain,ichain)]
                    total_dna_qual += qual
                    total_dna_degree += dna_degree*qual
                    total_dna_score += dna_score*qual

        modres_prop = modres_prop/float(len(structure_classification_map[m]))
            
        if total_lig_qual  > 0.:
            weighted_lig_degree = total_lig_degree/total_lig_qual
            weighted_lig_score = total_lig_score/total_lig_qual
        else:
            weighted_lig_degree = 0
            weighted_lig_score = 0.

        if total_metal_qual  > 0.:
            weighted_metal_degree = total_metal_degree/total_metal_qual
            weighted_metal_score = total_metal_score/total_metal_qual
        else:
            weighted_metal_degree = 0
            weighted_metal_score = 0.

        if total_ion_qual  > 0.:
            weighted_ion_degree = total_ion_degree/total_ion_qual
            weighted_ion_score = total_ion_score/total_ion_qual
        else:
            weighted_ion_degree = 0
            weighted_ion_score = 0.

        if total_prot_qual  > 0.:
            weighted_prot_degree = total_prot_degree/total_prot_qual
            weighted_prot_score = total_prot_score/total_prot_qual
        else:
            weighted_prot_degree = 0
            weighted_prot_score = 0.

        if total_rna_qual  > 0.:
            weighted_rna_degree = total_rna_degree/total_rna_qual
            weighted_rna_score = total_rna_score/total_rna_qual
        else:
            weighted_rna_degree = 0
            weighted_rna_score = 0.

        if total_dna_qual  > 0.:
            weighted_dna_degree = total_dna_degree/total_dna_qual
            weighted_dna_score = total_dna_score/total_dna_qual
        else:
            weighted_dna_degree = 0
            weighted_dna_score = 0.

        if best_res_iaa != None:
            best_res = best_res_iaa
        if best_res == None:
            best_res = best_res_icc
        if best_res == None:
            best_res = best_res_iaa_isc
        if best_res == None:
            best_res = best_res_isc


        class_dict[m] = (mut_class,conf,weighted_sc,conf_sc,best_res,max_seq_res,amount_of_structures,
                        weighted_c,conf_c,
                        weighted_d,conf_d,
                        weighted_r,conf_r,
                        weighted_l,conf_l,
                        weighted_m,conf_m,
                        weighted_i,conf_i,
                        weighted_h,conf_h,
                        max_seq_id,
                        weighted_raw,weighted_cent,weighted_norm,
                        weighted_lig_degree,weighted_lig_score,
                        weighted_metal_degree,weighted_metal_score,
                        weighted_ion_degree,weighted_ion_score,
                        weighted_prot_degree,weighted_prot_score,
                        weighted_rna_degree,weighted_rna_score,
                        weighted_dna_degree,weighted_dna_score,
                        weighted_modres,modres_prop
                        )
        t4 += time.time()

    #print 'createClassDict part 1: ',t1-t0
    #print 'createClassDict part 2: ',t2-t1
    #print 'createClassDict part 3: ',t3-t2
    #print 'createClassDict part 4: ',t4-t3

    t41 = time.time()
    for m in mutation_dict:
        if m in class_dict:
            continue
        aac,g,iupred,region,input_res_id = mutation_dict[m]
        #print glob
        if region == None:
            continue
        glob = True
        if region == 'disorder':
            glob = False
        elif region == 'globular':
            glob = True
        elif region in MobiDB_map:
            glob = False
        elif verbose:
            print('Unknown disorder region: ',region,' ,Mutation-ID: ',m)
        if not glob:
            mut_class = 'Disorder'
            conf = iupred
            weighted_sc = 'Surface'
            conf_sc = iupred
            best_res = None
            class_dict[m] = (mut_class,conf,weighted_sc,conf_sc,best_res,None,0,
                None,0,
                None,0,
                None,0,
                None,0,
                None,0,
                None,0,
                None,0,
                0.,
                None,None,None,
                None,None,
                None,None,
                None,None,
                None,None,
                None,None,
                None,None,
                None,None)
    t5 = time.time()

    #print 'createClassDict part 5: ',t5-t41

    return class_dict,structure_classification_map,interaction_dict

def getNeighborhood(class_dict,db,cursor,num_of_neighbors=10):
    neighborhood = {}
    if len(class_dict) == 0:
        return {}
    for m in class_dict:
        (mut_class,conf,weighted_sc,conf_sc,best_res,max_seq_res,amount_of_structures,
        weighted_c,conf_c,
        weighted_d,conf_d,
        weighted_r,conf_r,
        weighted_l,conf_l,
        weighted_m,conf_m,
        weighted_i,conf_i,
        weighted_h,conf_h,
        max_seq_id,
        weighted_raw,weighted_cent,weighted_norm,
        weighted_lig_degree,weighted_lig_score,
        weighted_metal_degree,weighted_metal_score,
        weighted_ion_degree,weighted_ion_score,
        weighted_prot_degree,weighted_prot_score,
        weighted_rna_degree,weighted_rna_score,
        weighted_dna_degree,weighted_dna_score,
        weighted_modres,modres_prop) = class_dict[m]
        if best_res == None:
            continue
        [r_id,qual,res_aa,res_nr,pdb_id,chain,resolution,cov,seq_id,rsa,min_lig,min_metal,min_ion,iacs] = best_res
        neighborhood[r_id] = []

    rows = ['Residue_1','Residue_2','Distance']
    table = 'RS_Residue_Residue'
    in_rows = {'Residue_1':list(neighborhood.keys())}
    results = select(db,cursor,rows,table,in_rows=in_rows)

    res_ids = set()
    res_bins = set()
    binsize = 500000

    for row in results:
        r_id_1 = row[0]
        r_id_2 = row[1]

        res_ids.add(r_id_2)
        bin_number = r_id_2//binsize
        res_bins.add(bin_number)

        d = row[2]
        #sorted inserting
        if neighborhood[r_id_1] == []:
            neighborhood[r_id_1].append([r_id_2,d])
        else:
            ins = False
            for i in range(len(neighborhood[r_id_1])):
                if d < neighborhood[r_id_1][i][1]:
                    neighborhood[r_id_1].insert(i,[r_id_2,d])
                    ins = True
                    break
            if not ins:
                neighborhood[r_id_1].append([r_id_2,d])

    residue_dict = {}
    if not len(res_ids) == 0:
        max_r = max(res_ids)
        min_r = min(res_ids)

        min_max_tuples = []
        if max_r - min_r > binsize:
            for bin_number in res_bins:
                min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize))
        else:
            min_max_tuples = [(min_r,max_r)]

        for (min_r,max_r) in min_max_tuples:
        #for [min_r,max_r] in bin_threshs:
            rows = ['Residue_Id','Structure','Number','Amino_Acid','Sub_Lig_Dist','Sub_Chain_Distances','Relative_Surface_Access']
            table = 'Residue'
            between_rows = {'Residue_Id':(min_r,max_r)}
            results = select(db,cursor,rows,table,between_rows=between_rows)

            if not results == ():
                for row in results:
                    r_id = row[0]
                    if not r_id in res_ids:
                        continue
                    residue_dict[r_id] = row[1:]

    for r1 in neighborhood:
        for neighbor in neighborhood[r1]:
            r2 = neighbor[0]
            neighbor += residue_dict[r2]

    return neighborhood

def excludeFarClasses(c,sc):
    if c == "Surface" or c == "Core" or c == 'Disorder':
        return c

    interactions = re.sub(r' Interaction$','',re.sub(r'^[^:]*: ','',c)).split(' and ')
    non_far_interactions = [ x for x in interactions if ' far' not in x ]

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

def simplifyClass(c,sc):
    if c == "Surface" or c == "Core" or c == 'Disorder':
        return c

    interactions = re.sub(r' Interaction$','',re.sub(r'^[^:]*: ','',c)).split(' and ')
    non_far_interactions = set([ x for x in interactions if ' far' not in x ])

    priorities = ['Metal', 'Ligand', 'DNA', 'RNA', 'Protein', 'Ion']

    if len(non_far_interactions) == 0:
        return sc

    for interaction in priorities:
        if interaction in non_far_interactions:
            return interaction + ' Interaction'

    print('Unknown Class: %s' % c)
    return sc

def writeInterFile(outfile,inter_dict,mutation_dict,gene_score_dict,new_aa_map,tag_map,class_dict,header=True):
    startline = "Uniprot\tAAC\tSpecie\tTag\tLigand_Interaction_Degree\tLigand_Interaction_Score\tChain_Interaction_Degree\tChain_Interaction_Score\tShort_Interaction_Degree\tShort_Interaction_Score\tMedium_Interaction_Degree\tMedium_Interaction_Score\tLong_Interaction_Degree\tLong_Interaction_Score\tClass\tComplex class"
    if header:
        lines = [startline]
    else:
        lines = []
    for m in inter_dict:
        aac,gene_id = mutation_dict[m][0:2]


        new_aa = new_aa_map[m]
        aac = "%s%s" % (aac.split(',')[0],new_aa)
        #(u_ac,u_id,species) = g_u_dict[mutation_dict[m][1]]
        (u_ac,gpan,u_id,error_code,error,species,gene_score) = gene_score_dict[gene_id]

        (Class,conf,weighted_sc,conf_sc,best_res,max_seq_res,amount_of_structures,
        weighted_c,conf_c,
        weighted_d,conf_d,
        weighted_r,conf_r,
        weighted_l,conf_l,
        weighted_m,conf_m,
        weighted_i,conf_i,
        weighted_h,conf_h,
        max_seq_id,
        weighted_raw,weighted_cent,weighted_norm,
        weighted_lig_degree,weighted_lig_score,
        weighted_metal_degree,weighted_metal_score,
        weighted_ion_degree,weighted_ion_score,
        weighted_prot_degree,weighted_prot_score,
        weighted_rna_degree,weighted_rna_score,
        weighted_dna_degree,weighted_dna_score,
        weighted_modres,modres_prop) = class_dict[m]
        simple_class = simplifyClass(Class,weighted_sc)

        if m in inter_dict:
            interstr = '\t'.join([str(x) for x in inter_dict[m]])
        else:
            interstr = '\t'.join((['None']*14))

        line = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (u_ac,aac,species,tag_map[m],interstr,simple_class,Class)
        lines.append(line)

    f = open(outfile,'a')
    f.write("\n".join(lines))
    f.close()


def writeClassFile(outfile,mutation_surface_dict,mutation_sec_dict,mutation_dict,gene_score_dict,class_dict,tag_map,header=True):
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

        (Class,conf,weighted_sc,conf_sc,best_res,max_seq_res,amount_of_structures,
        weighted_c,conf_c,
        weighted_d,conf_d,
        weighted_r,conf_r,
        weighted_l,conf_l,
        weighted_m,conf_m,
        weighted_i,conf_i,
        weighted_h,conf_h,
        max_seq_id,
        weighted_raw,weighted_cent,weighted_norm,
        weighted_lig_degree,weighted_lig_score,
        weighted_metal_degree,weighted_metal_score,
        weighted_ion_degree,weighted_ion_score,
        weighted_prot_degree,weighted_prot_score,
        weighted_rna_degree,weighted_rna_score,
        weighted_dna_degree,weighted_dna_score,
        weighted_modres,modres_prop) = class_dict[m]
        simple_class = simplifyClass(Class,weighted_sc)

        if best_res != None:
            [r_id,qual,res_aa,res_nr,pdb_id,chain,resolution,cov,seq_id,rsa,min_lig,min_metal,min_ion,iacs] = best_res

            recommended_structure = '%s:%s %s:%s' % (pdb_id,chain,res_nr,res_aa)
        else:
            resolution = '-'
            cov = '-'
            seq_id = '-'
            recommended_structure = '-'

        if max_seq_res != None:
            [max_seq_r_id,max_seq_qual,max_seq_res_aa,max_seq_res_nr,max_seq_pdb_id,max_seq_chain,max_seq_resolution,max_seq_cov,max_seq_seq_id,max_seq_rsa,max_min_lig,max_min_metal,max_min_ion,max_iacs] = max_seq_res

            max_seq_structure = '%s:%s %s:%s' % (max_seq_pdb_id,max_seq_chain,max_seq_res_nr,max_seq_res_aa)
        else:
            max_seq_resolution = '-'
            max_seq_cov = '-'
            max_seq_seq_id = '-'
            max_seq_structure = '-'

        aac,gene_id = mutation_dict[m][0:2]
        input_res_id = mutation_dict[m][4]
        if input_res_id == None:
            input_res_id = ''

        (u_ac,gpan,u_id,error_code,error,species,gene_score) = gene_score_dict[gene_id]

        input_pdb_id = ''
        if len(u_ac) == 6 and u_ac[4] == ':':
            input_pdb_id = u_ac
            u_ac = ''
        interaction_str = '-' # TODO

        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (u_ac,u_id,gpan,input_pdb_id,input_res_id,aac[0],aac[1:],species,tag_map[m],weighted_sc,Class,simple_class,interaction_str,str(conf),mv_sec_ass,recommended_structure,str(seq_id),str(cov),str(resolution),max_seq_structure,str(max_seq_seq_id),str(max_seq_cov),str(max_seq_resolution),str(amount_of_structures)))
    f = open(outfile,'a')
    f.write("\n".join(lines))
    f.close()

def writeStatFile(out_file,mutation_dict,class_dict,tag_map,stat_dict=None):
    seq_id_threshold = 0.95
    startline = 'Tag\tTotal proteins\tTotal positions\tMapped positions\tMapped into at least one corresponding structure (seq-id > %s%%)\tMapped only in homologs (seq-id <= %s%%)\tUnmapped, Disorder\tUnmapped, Globular' % (str(seq_id_threshold),str(seq_id_threshold))
    outmap = {'All':[set(),0,0,0,0,0,0]}
    if stat_dict != None:
        class_dict = stat_dict
        max_seq_key = 1
    else:
        max_seq_key = 21
    for m in tag_map:
        if tag_map[m] != None:
            raw_tags = tag_map[m].split(',')
        else:
            raw_tags = []
        tags = set(['All'])
        for tag in raw_tags:
            if tag[0] == '#':
                tag = tag[1:].split(':')[0]
            tags.add(tag)

        for tag in tags:
            if not tag in outmap:
                outmap[tag] = [set(),0,0,0,0,0,0]
            g = mutation_dict[m][1]
            outmap[tag][0].add(g)
            outmap[tag][1] += 1
            
            if m in class_dict:
                clas = class_dict[m][0]
                max_seq_id = class_dict[m][max_seq_key]
                if clas != 'Disorder' and clas != None:
                    outmap[tag][2] += 1
                    if max_seq_id > seq_id_threshold:
                        outmap[tag][5] += 1
                    else:
                        outmap[tag][6] += 1
                elif clas == 'Disorder':
                    outmap[tag][3] += 1
                else:
                    outmap[tag][4] += 1
            else:
                outmap[tag][4] += 1
    if None in outmap:
        del outmap[None]

    lines = [startline]
    for tag in outmap:
        tot_prot = len(outmap[tag][0])
        tot_pos = outmap[tag][1]
        mapped = outmap[tag][2]
        dis = outmap[tag][3]
        unmapped = outmap[tag][4]
        mapped_to_corr = outmap[tag][5]
        mapped_to_homolog = outmap[tag][6]

        if float(tot_pos) == 0.0:
            continue
        line = '%s\t%s\t%s\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)\t%s (%s%%)' % (tag,str(tot_prot),str(tot_pos),str(mapped),str(100.*float(mapped)/float(tot_pos)),str(mapped_to_corr),str(100.*float(mapped_to_corr)/float(tot_pos)),str(mapped_to_homolog),str(100.*float(mapped_to_homolog)/float(tot_pos)),str(dis),str(100.*float(dis)/float(tot_pos)),str(unmapped),str(100.*float(unmapped)/float(tot_pos)))
        lines.append(line)
    f = open(out_file,'w')
    f.write("\n".join(lines))
    f.close()

#Is this still used???
#called by output
def prodAnoOut(output_file,session_id,db_name,db_adress,db_user_name,db_password,top=None,ligand_filter=None,proteome=False):
    if ligand_filter != None:
        f = open(ligand_filter,'r')
        lines = f.readlines()
        f.close()
        ligand_filter = set()
        for line in lines:
            ligand_filter.add(line.replace('\n','').replace(' ',''))    

    db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    cursor = db.cursor()

    mutation_id_list = set()
    template_id_list = set()

    print(db_name,session_id)

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

    sql = "SELECT Gene FROM RS_Gene_Session WHERE Session = '%s'" % str(session_id)
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in prodAnoOut: %s,\n%s" % (sql,f))

    gene_ids = []
    for row in results:
        gene_ids.append(str(row[0]))

    binsize = 500000
    """
    mutation_id_list = set()
    max_mut_id = 0
    min_mut_id = None

    
    bins = set()

    for row in results:
        mut_id = row[0]

        bin_number = mut_id//binsize
        bins.add(bin_number)

        if mut_id > max_mut_id:
            max_mut_id = mut_id
        if min_mut_id == None or mut_id < min_mut_id:
            min_mut_id = mut_id
        mutation_id_list.add(mut_id)
    """
    m_r_map = {}
    res_bins = set()
    res_ids = set()

    t3 = time.time()
    #print "Time for prodAnoOut part 3: ",t3-t2

    if not len(tag_map) == 0:
        
        sql = "SELECT Mutation,Residue FROM RS_Mutation_Residue WHERE Gene in (%s)" % (','.join(gene_ids))
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
                m_id = row[0]
                if not m_id in tag_map:
                    continue
                r_id = row[1]

                bin_number = r_id//binsize
                res_bins.add(bin_number)

                res_ids.add(r_id)
                if not m_id in m_r_map:
                    m_r_map[m_id] = set()
                m_r_map[m_id].add(r_id)

    #print m_r_map

    residue_dict = {}
    s_ids = set()
    if not len(res_ids) == 0:
        max_r = max(res_ids)
        min_r = min(res_ids)

        min_max_tuples = []
        if max_r - min_r > binsize:
            for bin_number in res_bins:
                min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize))
        else:
            min_max_tuples = [(min_r,max_r)]

        #print min_max_tuples
        
        for (min_r,max_r) in min_max_tuples:
            sql = """SELECT Residue_Id,
                            Structure,
                            Number,
                            Sub_Lig_Dist,
                            Sub_Chain_Distances,
                            Relative_Surface_Access,
                            Secondary_Structure_Assignment,
                            Homomer_Distances,
                            Interaction_Profile,
                            FROM Residue WHERE Residue_Id BETWEEN %s AND %s""" % (str(min_r),str(max_r))
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
                    r_id = row[0]
                    if not r_id in res_ids:
                        continue
                    s_id = row[1]
                    s_ids.add(s_id)
                    residue_dict[r_id] = row[1:]

    #print residue_dict

    mutation_dict = getMutationDict(set(tag_map.keys()),db,cursor)
    gene_id_list = set()
    for m in mutation_dict:
        gene_id_list.add(mutation_dict[m][1])

    gene_score_dict = getGeneScoreDict(gene_id_list,session_id,db,cursor)
    structure_dict = getStructureDict(s_ids,db,cursor)

    seq_id_map = getSeqIds(s_ids,gene_id_list,db,cursor)

    cursor.close()
    db.close()

    
    lines = produceAnotationOutput(session_id,residue_dict,gene_score_dict,structure_dict,mutation_dict,ligand_filter,new_aa_map,tag_map,m_r_map,seq_id_map)

    #(gene_name,species,sp_id,gpan,aachange,pdb_id,residue,chain,ligtext,sub_lig_dist,lig_sub_dist,sub_chain_dists,chain_type,chain_sub_dist,str(chains),str(rel_sur_acc),sec_str_ass,mut_class,str(chemical_distance),str(blosum_value),tag)

    startline = "\t".join(("Gene",'Species',"SP ID","REFSEQS","Position","PDB","Residue_id","Chain","Sequence Identity","Coverage","Ligands(Name,Residue Nr)","Sub-Lig-Dists","Min Sub-Lig-Dist","Sub Chain Distances","Min Sub-Chain-Dist","Chain Type","Relative Surface Access","Secondary Structure Assignment","Class","Chemical Distance","Blosum62 Value","Tag"))
    #print len(lines)
    lines2 = [startline]
    for line in lines:
        try:
            line2 = "\t".join(line)
            lines2.append(line2)
        except:
            print("Error in prodAnoOut: ",line)
    page = "\n".join(lines2)
    f = open(output_file, "wb")
    f.write(page)
    f.close()

    print("%s done" % output_file)
    return output_file

def produceAnotationOutput(session_id,residue_dict,gene_score_dict,structure_dict,mutation_dict,ligand_filter,new_aa_map,tag_map,m_r_map,seq_id_map):
    global chem_dist_matrix
    outlines = []
    #print mutation_dict
    for m_id in mutation_dict:
        if not m_id in m_r_map:
            continue #unmapped mutations
        (aachange,gene_id,iupred_score,glob,input_res_id) = mutation_dict[m_id]
        tag = tag_map[m_id]
        new_aa = new_aa_map[m_id]
        if ord(aachange[-1]) >47 and ord(aachange[-1]) < 59:
            aachange = "%s%s" % (aachange.split(',')[0],new_aa)
        else:
            aachange = "%s%s" % (aachange.split(',')[0][:-1],new_aa)
        (gene_name,gpan,sp_id,error_code,error,species,gene_score) = gene_score_dict[gene_id]
        for r_id in m_r_map[m_id]:
            row = residue_dict[r_id]
            s_id = row[0]
            residue = row[1]
            sub_lig_dist = row[2]
            sub_chain_dists = row[3]
            rel_sur_acc = row[4]
            sec_str_ass = row[5]
            homomer_dists = row[6]
            structure = structure_dict[s_id]
            
            pdb_id = structure[0]
            chain = structure[1]
            resolution = structure[2]
            ligands = structure[3]
            oligos = structure[4]
            chains = structure[5]

            (seq_id,coverage) = seq_id_map[gene_id][s_id]

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
                interacting_chain_id = ds[0][1][0]
            else:
                chain_sub_dist = "NONE"
                interacting_chain_id = "NONE"

        
            ligand_texts = []        
            for ligand in ligands:
                if ligand[0] == "Ligand":
                    t = "%s,%s" %(ligand[1],ligand[2])
                    ligand_texts.append(t)
            ligtext = ";".join(ligand_texts)

            chain_type_dict = {}
            for c_pair in chains.split(","):
                chain_type_dict[c_pair.split(":")[0]] = c_pair.split(":")[1]
            if interacting_chain_id == "NONE":
                class_chain_type = "NONE"
            else:
                class_chain_type = chain_type_dict[interacting_chain_id]

            mut_class = getClass(lig_sub_dist,chain_sub_dist,class_chain_type,rel_sur_acc)

            chemical_distance = getChemicalDistance(aachange)
            blosum_value = getBlosumValue(aachange)        

            if gpan == None:
                gpan = 'None'

            line = (gene_name,species,sp_id,gpan,aachange,pdb_id,residue,chain,str(seq_id),str(coverage),ligtext,str(sub_lig_dist),str(lig_sub_dist),sub_chain_dists,str(chain_sub_dist),str(chains),str(rel_sur_acc),str(sec_str_ass),mut_class,str(chemical_distance),str(blosum_value),tag)
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

    #print "Time for calculating Gene Scores: ",t1-t0
    #print "Time for updating Gene Scores into Database: ",t2-t1

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

    go_list = sorted(list(go_dict.values()),key = lambda x: x[4],reverse = True)
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

    path_list = sorted(list(path_dict.values()),key = lambda x: x[4],reverse = True)
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

def getStructureClassification(pdb_id,chain,db,cursor):
    sql = "SELECT Structure_Id,Chains FROM Structure WHERE PDB = '%s' AND Chain = '%s'" % (pdb_id,chain)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getStructureClassification: %s" % sql)
    if results == ():
        print("Did not found the Structure")
        return {}
    s_id = results[0][0]
    chains = results[0][1]

    #print s_id,chains

    sql = "SELECT Number,Sub_Lig_Dist,Sub_Chain_Distances,Relative_Surface_Access FROM Residue WHERE Structure = '%s'" % str(s_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in getStructureClassification: %s" % sql)

    chain_type_dict = {}
    for c_pair in chains.split(","):
        chain_type_dict[c_pair.split(":")[0]] = c_pair.split(":")[1]
    res_class_map = {}

    for row in results:

        residue = row[0]
        sub_lig_dist = row[1]
        sub_chain_dists = row[2]
        rel_sur_acc = row[3]
        

        ds = []
        if sub_lig_dist != None:
            lig_dists = sub_lig_dist.split(",")
            for lig_dist in lig_dists:
                ldp = lig_dist.split(":")
                if len(ldp) > 1:
                    lig_name = ldp[0].split('_')[0]

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

        if chain_type == "NONE":
            class_chain_type = chain_type
        else:
            class_chain_type = chain_type_dict[chain_type[0]]


        mut_class = getClass(lig_sub_dist,chain_sub_dist,class_chain_type,rel_sur_acc)
        res_class_map[residue] = mut_class
    return res_class_map

def structuralCoverageDistribution(session_id,db,cursor):

    sql = "SELECT Mutation FROM RS_Mutation_Session WHERE Session = '%s'" % str(session_id)
    try:
        cursor.execute(sql)
        results = cursor.fetchall()
        db.commit()
    except:
        raise NameError("Error in structuralCoverageDistribution: %s" % sql)


    mutation_id_list = set()
    max_mut_id = 0
    min_mut_id = None

    binsize = 200000
    bins = set()

    for row in results:
        mut_id = row[0]

        bin_number = mut_id//binsize
        bins.add(bin_number)

        if mut_id > max_mut_id:
            max_mut_id = mut_id
        if min_mut_id == None or mut_id < min_mut_id:
            min_mut_id = mut_id
        mutation_id_list.add(mut_id)

    m_r_map = {}
    


    min_max_tuples = []
    if max_mut_id - min_mut_id > binsize:
        for bin_number in bins:
            min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize))
    else:
        min_max_tuples = [(min_mut_id,max_mut_id)]

    
    for (min_mut_id,max_mut_id) in min_max_tuples:
        results = ()
        sql = "SELECT Mutation,Residue FROM RS_Mutation_Residue WHERE Mutation BETWEEN %s AND %s" % (str(min_mut_id),str(max_mut_id))
        try:
            # Execute the SQL command
            cursor.execute(sql)
            results = cursor.fetchall()
            # Commit your changes in the database
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in structuralCoverageDistribution: %s,\n%s" % (sql,f))

        if not results == ():
            for row in results:
                m_id = row[0]
                if not m_id in mutation_id_list:
                    continue
                r_id = row[1]

                if not m_id in m_r_map:
                    m_r_map[m_id] = set()
                m_r_map[m_id].add(r_id)

    histo = [0]
    for m in mutation_id_list:
        if not m in m_r_map:
            histo[0] += 1
        else:
            size = len(m_r_map[m])
            if size >= len(histo):
                histo += [0]*(1+size-len(histo))

            histo[size] += 1

    return histo

def proteinStructureCoverage(session_id,db,cursor):
    sql = "SELECT Gene FROM RS_Gene_Session WHERE Session = '%s'" % str(session_id)
    try:
        # Execute the SQL command
        cursor.execute(sql)
        results = cursor.fetchall()
        # Commit your changes in the database
        db.commit()
    except:
        [e,f,g] = sys.exc_info()
        raise NameError("Error in proteinStructureCoverage: %s,\n%s" % (sql,f))

    gene_ids = set()

    max_g_id = 0
    min_g_id = None

    binsize = 5000
    bins = set()

    for row in results:
        g_id = row[0]

        bin_number = g_id//binsize
        bins.add(bin_number)

        if g_id > max_g_id:
            max_g_id = g_id
        if min_g_id == None or g_id < min_g_id:
            min_g_id = g_id
        gene_ids.add(g_id)


    min_max_tuples = []
    if max_g_id - min_g_id > binsize:
        for bin_number in bins:
            min_max_tuples.append((bin_number*binsize,(bin_number+1)*binsize))
    else:
        min_max_tuples = [(min_g_id,max_g_id)]

    gr_35 = set()
    gr_95 = set()
    for (min_g_id,max_g_id) in min_max_tuples:

        sql = "SELECT Gene,Sequence_Identity FROM Alignment WHERE Gene BETWEEN %s AND %s" % (min_g_id,max_g_id)
        try:
            # Execute the SQL command
            cursor.execute(sql)
            results = cursor.fetchall()
            # Commit your changes in the database
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in proteinStructureCoverage: %s,\n%s" % (sql,f))

        for row in results:
            g_id = row[0]
            if not g_id in gene_ids:
                continue
            seq_id = float(row[1])

            #print seq_id

            if seq_id > 0.98:
                gr_95.add(g_id)
            gr_35.add(g_id)

    #print len(gene_ids)

    no_str = set()
    for g_id in gene_ids:
        if not g_id in gr_35:
            no_str.add(g_id)

    disorder_map = {}

    for (min_g_id,max_g_id) in min_max_tuples:
        sql = "SELECT Gene,IUPRED_Glob FROM Mutation WHERE Gene BETWEEN %s AND %s"  % (min_g_id,max_g_id)
        try:
            # Execute the SQL command
            cursor.execute(sql)
            results = cursor.fetchall()
            # Commit your changes in the database
            db.commit()
        except:
            [e,f,g] = sys.exc_info()
            raise NameError("Error in proteinStructureCoverage: %s,\n%s" % (sql,f))

        for row in results:
            g_id = row[0]
            if not g_id in no_str:
                continue
            if not g_id in disorder_map:
                disorder_map[g_id] = [0,0]
            try:
                glob = int(row[1])
                if glob == 0:
                    disorder_map[g_id][0] += 1
                elif glob == 1:
                    disorder_map[g_id][1] += 1
            except:
                pass

    majorly_disordered = 0
    for g_id in disorder_map:
        disorder,order = disorder_map[g_id]
        if disorder == 0:
            continue
        if float(disorder)/float(order+disorder) >= 0.8:
            majorly_disordered += 1

    print("Structure: ",len(gr_95))
    print("Homolog: ",len(gr_35)-len(gr_95))
    print("Disordered: ",majorly_disordered)
    print("No Structure: ",len(no_str) - majorly_disordered)




