#!/usr/bin/python3
import sys
import getopt
import os

import database
import gzip

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
    'TRUNCATE Complex;',
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
        db.close()
        empty(config)
        return


    if config.database_source_path[-3:] == '.gz':
        f = gzip.open(config.database_source_path,'r')
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
