#!/usr/bin/python3
#curl --metalink ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz > uniprot_trembl.fasta.gz
#curl --metalink ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot_varsplic.fasta.gz > uniprot_sprot_varsplic.fasta.gz
#curl --metalink ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz > uniprot_sprot.fasta.gz

import pymysql as MySQLdb
import os
import sys
import getopt
import gzip

sys.path.append('/wibicom/SHARED_DATA/agress/structman')
sys.path.append('/wibicom/SHARED_DATA/agress/structman/lib')
import structman
import database

def parseSeqFile(seq_file_path,config):
    f = gzip.open(seq_file_path,'r')
    lines = f.read().decode('ascii').split('\n')
    f.close()

    seq_map = {}
    for line in lines:
        if line == '':
            continue
        if line[0] == '>':
            u_ac = line.split('|')[1]
            seq_map[u_ac] = ''
        else:
            seq_map[u_ac] += (line)

    values = []
    for u_ac in seq_map:
        values.append((u_ac,seq_map[u_ac]))

    db = MySQLdb.connect(config.db_adress,config.db_user_name,config.db_password,config.mapping_db)
    cursor = db.cursor()

    database.insert(db,cursor,'Sequences',['Uniprot_Ac','Sequence'],values)

    db.close()

if __name__ == "__main__":

    argv = sys.argv[1:]
    try:
        opts,args = getopt.getopt(argv,"c:i:f:",[])
    except getopt.GetoptError:
        print("Illegal Input")
        sys.exit(2)

    infolder = None
    infile = None
    for opt,arg in opts:
        if opt == '-c':
            config = arg
        if opt == '-f':
            infolder = arg
        if opt == '-i':
            infile = arg

    if infolder == None:
        #infolder = os.path.abspath(sys.argv[0]).rsplit('/',1)[0]
        pass #for now
    else:
        infolder = os.path.abspath(infolder)

    config = os.path.abspath(config)

    config = structman.Config(config)

    if infile != None:
        parseSeqFile(infile,config)

    if infolder != None:
        seq_files = ['uniprot_sprot_varsplic.fasta.gz','uniprot_sprot.fasta.gz','uniprot_trembl.fasta.gz',]

        for seq_file in seq_files:
            seq_file_path = '%s/%s' % (infolder,seq_file)
            parseSeqFile(seq_file_path,config)
    



