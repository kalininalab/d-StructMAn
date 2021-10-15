#!/usr/bin/python2
import MySQLdb
import os
import sys
import getopt
import gzip

def parseConfig(config):

    global db_adress
    global db_user_name
    global db_password

    
    f = open(config, "r")
    lines = f.read().split('\n')
    f.close()
    for line in lines:
        if len(line) == 0:
            continue
        if line[0] == '#':
            continue
        words = line.split("=")
        if len(words) > 1:
            opt = words[0]
            arg = words[1].replace("\n","")
            
            if opt == 'db_adress':
                db_adress = arg
            elif opt == 'db_user_name':
                db_user_name = arg
            elif opt == 'db_password':
                db_password = arg


if __name__ == "__main__":

    argv = sys.argv[1:]
    try:
        opts,args = getopt.getopt(argv,"c:i:",[])
    except getopt.GetoptError:
        print("Illegal Input")
        sys.exit(2)

    for opt,arg in opts:
        if opt == '-c':
            config = arg
        if opt == '-i':
            infile = arg

    infile = os.path.abspath(infile)
    config = os.path.abspath(config)

    parseConfig(config)

    db_name = 'struct_man_db_uniprot'

    db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    cursor = db.cursor()

    if infile[-3:] == '.gz':
        f = gzip.open(infile,'r')
        lines = f.readlines()
        f.close()
    else:
        f = open(infile,'r')
        lines = f.readlines()
        f.close()

    gene_sequence_db = {}
    for line in lines:
        if line[0] == '>':
            u_ac = line.split('|')[1]
            gene_sequence_db[u_ac] = ''
        else:
            gene_sequence_db[u_ac] += line[:-1]

    value_strs = []
    for gene in gene_sequence_db:
        value_strs.append("('%s','%s')" % (gene.replace("'","\\'"),gene_sequence_db[gene]))
        

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

            sql = ("""INSERT INTO Sequences(Uniprot_Ac,
                    Sequence)
                    VALUES %s""") % ','.join(part)
            try:
                cursor.execute(sql)
                db.commit()
            except:
                [e,f,g] = sys.exc_info()
                raise NameError("Couldn't insert Sequence: %s,%s" % (e,f))

    db.close()


