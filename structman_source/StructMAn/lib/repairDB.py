#!/usr/bin/python3
import sys
import getopt
import os
import traceback
import database
import gzip

import subprocess

main_file_path = os.path.abspath(sys.argv[0])
sys.path.append(main_file_path.rsplit('/',2)[0])
import structman

#Tries to reclassify all positions with no classification in the database
#Usefull when,
#   - There was an error in insertClassifications
def reclass_null(config):
    db = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,config.db_name)
    cursor = db.cursor()
    #get all positions without classifications (together with their protein id)
    columns = ['Mutation_Id','Gene']
    table = 'Mutations'
    null_columns = set(['Class'])
    results = database.select(db,cursor,columns,table,null_columns=null_columns)

    prot_ids = set()
    for row in results:
        prot_ids.add(row[1])

    #for all proteins search for alignments in the database

    rows = ['Gene','Structure','Alignment']
    table = 'Alignment'
    results = database.binningSelect(prot_ids,rows,table,db,cursor)

    #try to map the positions into the alignments

    #for all positions, which could be mapped, reclassify

    return

def reset(cursor,keep_structures = False):

    sql_commands = ['SET FOREIGN_KEY_CHECKS=0;', 
    'TRUNCATE Gene;',
    'TRUNCATE Mutation;',
    'TRUNCATE GO_Term;',
    'TRUNCATE Pathway;',
    'TRUNCATE Session;',
    'TRUNCATE Alignment;',
    'TRUNCATE RS_Gene_Session;',
    'TRUNCATE RS_Gene_GO_Term;',
    'TRUNCATE RS_Mutation_Session;',
    'TRUNCATE RS_Gene_Pathway;',
]

    if not keep_structures:
        sql_commands += [
        'TRUNCATE Ligand;',
        'TRUNCATE Structure;',
        'TRUNCATE Residue;',
        'TRUNCATE Complex;',
        'TRUNCATE RS_Ligand_Structure;',
        'TRUNCATE RS_Residue_Residue;'
        ]

    sql_commands.append('SET FOREIGN_KEY_CHECKS=1;')

    for sql in sql_commands:
        try:
            # Execute the SQL command
            cursor.execute(sql)
            # Commit your changes in the database
            #db.commit()
        except:
            # Rollback in case there is any error
            #db.rollback()
            [e,f,g] = sys.exc_info()
            g = traceback.format_exc(g)
            print("Error: ",e,f,g)
    return

def reduceToStructures(config,infile):
    f = open(infile,'r')
    lines = f.readlines()
    f.close()

    complexes = set()
    for line in lines:
        words = line.replace('\t',' ').split()
        pdb_tuple = words[0]
        if pdb_tuple.count(':') != 1:
            continue
        pdb_id = pdb_tuple.split(':')[0]
        complexes.add(pdb_id)

    table = 'Structure'
    rows = ['Structure_Id','PDB']
    results = database.select(config,rows,table)

    list_of_doom = []
    for row in results:
        if row[1] in complexes:
            continue
        list_of_doom.append(row[0])

    results = database.select(config,['Complex_Id','PDB'],'Complex')

    list_of_doom_c = []
    for row in results:
        if row[1] in complexes:
            continue
        list_of_doom_c.append(row[0])

    statement = 'DELETE FROM Residue WHERE Structure IN (%s)' % (','.join(['%s']*len(list_of_doom)))
    db,cursor = config.getDB()

    cursor.execute(statement,list_of_doom)
    results = cursor.fetchall()
    db.commit()

    db.close()

    statement = 'DELETE FROM Structure WHERE Structure_Id IN (%s)' % (','.join(['%s']*len(list_of_doom)))
    db,cursor = config.getDB()

    cursor.execute(statement,list_of_doom)
    results = cursor.fetchall()
    db.commit()

    db.close()

    statement = 'DELETE FROM Complex WHERE Complex_Id IN (%s)' % (','.join(['%s']*len(list_of_doom_c)))
    db,cursor = config.getDB()

    cursor.execute(statement,list_of_doom_c)
    results = cursor.fetchall()
    db.commit()

    db.close()
    return

def export(config,target_path):
    target_path = '%s/%s.sql.gz' % (os.path.abspath(target_path),config.db_name)
    cmds = ' '.join(['mysqldump','-n','-u',config.db_user_name,'-h',config.db_adress,'--password=%s' % config.db_password,'--single-transaction','--databases',config.db_name,
            '|','pigz','--best','-p','8','>','"%s"' % target_path])

    p = subprocess.Popen(cmds,shell = True)
    p.wait()

    return

def empty(config):
    db,cursor = config.getDB()
    reset(cursor)
    db.close()
    return

def clear(config):
    db,cursor = config.getDB()
    reset(cursor,keep_structures = True)
    db.close()
    return

def destroy(config):
    db,cursor = config.getDB()
    sql = 'DROP DATABASE %s' % config.db_name
    cursor.execute(sql)
    db.close()

def load(config):
    db,cursor = config.getDB(server_connection = True)

    sql = 'CREATE DATABASE %s' % config.db_name

    try:
        cursor.execute(sql)
        db.close()
    except:
        [e,f,g] = sys.exc_info()
        g = traceback.format_exc()
        print('\n'.join([str(e),str(f),str(g)]))
        db.close()
        #empty(config)
        return

    new_lines = []
    if config.database_source_path[-3:] == '.gz':
        f = gzip.open(config.database_source_path,'rb')
        binary = True
    else:
        f = open(config.database_source_path,'r')
        binary = False
    lines = f.readlines()
    f.close()

    for line in lines:
        if binary:
            line = line.decode('ascii')
        if line[:4] == 'USE ':
            new_lines.append('USE `%s`;\n' % config.db_name)
        else:
            new_lines.append(line)

    db_file = 'tmp_database_file.sql'
    f = open(db_file,'w')
    f.write(''.join(new_lines))
    f.close()

    cmds = ' '.join(['mysql','-u',config.db_user_name,'-h',config.db_adress,'--password=%s' % config.db_password,'<',db_file])

    p = subprocess.Popen(cmds,shell = True)
    p.wait()

    os.remove(db_file)

    #Old version, remove if new version is working
    '''
    if config.database_source_path[-3:] == '.gz':
        f = gzip.open(config.database_source_path,'r')
        text = f.read().decode('ascii')
    else:
        f = open(config.database_source_path,'r')
        text = f.read()
    f.close()

    db,cursor = config.getDB()

    for cmd in text.split(';'):
        if cmd == '':
            continue
        try:
            cursor.execute('%s;' % cmd)
        except:
            pass

    db.close()
    '''
    return

#destroys and reinitializes the database
def reinit(config):
    empty(config)
    destroy(config)
    load(config)

def main(config,reclassify_null=False):
    #called with --rcn
    if reclassify_null:
        reclass_null(config)

    return

if __name__ == "__main__":
    argv = sys.argv[1:]
    try:
        opts,args = getopt.getopt(argv,"c:",['rcn'])
    except getopt.GetoptError:
        print("Illegal Input")
        sys.exit(2)

    reclassify_null = False
    config_path = None

    for opt,arg in opts:
        if opt == '-c':
            config_path = arg
        elif opt == '--rcn':
            reclassify_null = True

    config = structman.Config(config_path,None,None,None,None,False)

    main(config,reclassify_null=reclassify_null)
