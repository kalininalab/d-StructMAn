#!/usr/bin/python2
#This is the wrapper script, including and automatizing all functions StructMAn has to offer
import sys
import getopt
import os
import time
import MySQLdb

cwd = os.getcwd()
db_adress = ""
db_user_name = ""
db_password = ""
db_name = ""
config = ''
infile = ''
overwrite = False
annotationtable = False
num_of_cores = None


def parseConfig(config):
    global db_adress
    global db_user_name
    global db_password
    global db_name
    global annotationtable

    f = open(config,'r')
    lines = f.read().split('\n')
    f.close()

    for line in lines:
        words = line.split('=')
        if words[0] == 'db_adress':
            db_adress = words[1]
        if words[0] == 'db_user_name':
            db_user_name = words[1]
        if words[0] == 'db_password':
            db_password = words[1]
        if words[0] == 'db_name':
            db_name = words[1]
        if words[0] == 'do_anno':
            if words[1] == 'True':
                annotationtable=True

def main(db_name,db_adress,db_password,db_user_name,infile,main_file_path,config,overwrite=False,anno=False,intertable=False):

    #check if the infile got already processed
    db = MySQLdb.connect(db_adress,db_user_name,db_password,db_name)
    cursor = db.cursor()
    session_id = database.getSessionId(infile,db,cursor)
    db.close()

    #create the output folder
    [trunk,infilename] = infile.rsplit('/',1)
    out_folder = "%s/Output" % trunk
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    #run the main pipeline
    if session_id == None:
        session_id = serializedPipeline.main(infile,config,out_folder,main_file_path,num_of_cores)

    #run the output scripts
    session_name = infilename.rsplit('.',1)[0]
    months = {1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}
    time_struct = time.gmtime()
    year = str(time_struct[0])
    month = months[time_struct[1]]
    day = str(time_struct[2])
    date = "%s%s%s" % (day,month,year)
    outpath = "%s/%s_%s" % (out_folder,session_name,date)
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    output.main(db_name,db_adress,db_password,db_user_name,session_id,outpath,config_path=config,overwrite=overwrite,anno=anno,intertable=intertable)


if __name__ == "__main__":

    main_file_path = os.path.abspath(sys.argv[0])

    sys.path.append("%s/lib" % main_file_path.rsplit('/',1)[0])

    import output
    import database
    import serializedPipeline

    argv = sys.argv[1:]
    try:
        opts,args = getopt.getopt(argv,"c:i:n:",['overwrite','annotationtable'])
    except getopt.GetoptError:
        print("Illegal Input")
        sys.exit(2)

    for opt,arg in opts:
        if opt == '-c':
            config = arg
        if opt == '-i':
            infile = arg
        if opt == '-n':
            num_of_cores = int(arg)
        if opt == '--overwrite':
            overwrite = True
        if opt == '--annotationtable':
            annotationtable = True

    #make infile to global path
    infile = os.path.abspath(infile)

    if config == '':
        if os.path.exists('%s/config.txt' % infile.rsplit('/',1)[0]):
            config = '%s/config.txt' % infile.rsplit('/',1)[0]
        elif os.path.exists('%s/config.txt' % cwd):
            config = '%s/config.txt' % cwd
        elif os.path.exists('%s/config.txt' % main_file_path.rsplit('/',1)[0]):
            config = '%s/config.txt' % main_file_path.rsplit('/',1)[0]
        else:
            print("No config file found, please use -c [Path to config]")
            sys.exit(2)

    print("Using following config file: %s" % config)
    config = os.path.abspath(config)

    parseConfig(config)

    main(db_name,db_adress,db_password,db_user_name,infile,main_file_path,config,overwrite=overwrite,anno=annotationtable)
