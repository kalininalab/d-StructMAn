#!/usr/bin/python3
import sys
import getopt
import os

import database

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
