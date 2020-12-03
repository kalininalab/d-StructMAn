#!/bin/bash
# Description: Docker entrypoint script to configure the StructMAn container. All ENV variables are set in the "docker-compose.yml" file
# Author: Sanjay kumar Srikakulam (sanjaysrikakulam@gmail.com)

set -euo pipefail

echo "***** Container configuration starting for the host $HOSTNAME *****"
cat "/usr/structman_library/sources/StructMAn_db/db_split*" > "/usr/structman_library/sources/StructMAn_db/struct_man_db.sql.gz"
echo "Cat done"
# Get the default StructMAn database name from the source file stored under /usr/structman_library/sources/StructMAn_db/
if [[ $(ls -L /usr/structman_library/sources/StructMAn_db/*.sql.gz | wc -l) == 1 ]]; then
    for sql_file in /usr/structman_library/sources/StructMAn_db/*.sql.gz; do
        STRUCTMAN_DB=$(basename $sql_file | cut -d "." -f 1)
        STRUCTMAN_DB_FOLDER_NAME=$(cat $sql_file | grep "Database" | awk '{print $NF}' | tr -d '`')
    done
elif [[ $(ls -L /usr/structman_library/sources/StructMAn_db/*.sql.gz | wc -l) == 0 ]]; then
    echo "===>    There does not seem to be any default StructMAn database file in the source directory (/usr/structman_library/sources/StructMAn_db/). Aborting the container setup!    <==="
    exit 1
fi

# Adding "structman.py" to "/usr/local/bin" as a symlink to make it a command line utility
configure_structman() {
    if [[ -f /usr/structman_library/sources/StructMAn/structman.py ]]; then
        if [[ ! -e /usr/local/bin/structman.py ]]; then
            ln -s /usr/structman_library/sources/StructMAn/* /usr/local/bin/
        fi
    else
        echo "===>    **structman.py** script could not be found. Aborting the container setup!    <==="
        exit 1
    fi
}

# Create necessary directories for the MySQL server setup with proper permissions and ownership
create_mysql_dirs(){
    # Setup the MySQL data directory
    mkdir -p /var/lib/mysql
    chmod -R 777 /var/lib/mysql
    chown -R mysql:mysql /var/lib/mysql

    # Setup the MySQL run directory
    mkdir -p /var/run/mysqld
    chmod -R 777 /var/run/mysqld
    rm -rf var/run/mysqld/mysqld.sock.lock
    chown -R mysql:mysql /var/run/mysqld

    # Setup the MySQL log directory
    mkdir -p /var/log/mysql
    chmod -R 777 /var/log/mysql
}

# Remove old configuration files and make the new ones
configure_mysql() {
    rm -rf /etc/mysql/my.cnf &>/dev/null
    rm -rf /etc/mysql/mysql.conf.d/mysqld.cnf &>/dev/null
    rm -rf /etc/mysql/mysql.cnf &>/dev/null

    # Remove debian system's maintenance password
    sed 's/password = .*/password = /g' -i /etc/mysql/debian.cnf

    # Configure MySQL server
    if [[ ! -f /etc/mysql/my.cnf ]]; then
        echo -e "# Include the configuration files from these directories
!includedir /etc/mysql/conf.d/
!includedir /etc/mysql/mysql.conf.d/" > /etc/mysql/my.cnf
    fi

    if [[ ! -f /etc/mysql/mysql.conf.d/mysqld.cnf ]]; then
        echo -e "# MySQL Server configuration
[mysqld_safe]
socket = /var/run/mysqld/mysqld.sock
nice = 0

[mysqld]
skip-host-cache
skip-external-locking
skip-name-resolve
datadir = /var/lib/mysql
tmpdir = /tmp
port = 3306
bind-address = 0.0.0.0
user = mysql
socket = /var/run/mysqld/mysqld.sock
pid-file = /var/run/mysqld/mysqld.pid
log-error = /var/log/mysql/error.log
basedir = /usr
lc-messages-dir = /usr/share/mysql
init_connect = 'SET NAMES utf8'
character-set-server = utf8
collation-server = utf8_unicode_ci
key_buffer_size = 16M
max_allowed_packet = 16M
thread_stack = 192K
thread_cache_size = 8
myisam-recover-options = BACKUP
query_cache_limit = 1M
query_cache_size = 16M
expire_logs_days = 10
max_binlog_size = 100M

[client]
default-character-set = utf8
socket = /var/run/mysqld/mysqld.sock

[mysql]
default-character-set = utf8
socket = /var/run/mysqld/mysqld.sock" > /etc/mysql/mysql.conf.d/mysqld.cnf
    fi

    echo "===>    MySQL server configuration is done    <==="
}

# Initialize MySQL and ceate a new MySQL user for the application
initialize_mysql_and_create_users() {
    echo "===>    Initializing MySQL database server!       <==="
    mysqld --initialize-insecure --user=mysql > /dev/null 2>&1 &
    sleep 60

    # Start mysql server
    /usr/bin/mysqld_safe >/dev/null 2>&1 &
    sleep 60

    # Create a debian-sys-maint user and a StructMAn db user
    echo "===>    Creating a debian system maintenance user and a MySQL StructMAn user  <==="

    mysql --user=root << EOF
CREATE USER IF NOT EXISTS 'debian-sys-maint'@'localhost' IDENTIFIED BY '';
GRANT ALL PRIVILEGES on *.* TO 'debian-sys-maint'@'localhost' IDENTIFIED BY '' WITH GRANT OPTION;

CREATE USER IF NOT EXISTS "$MYSQL_STRUCTMAN_USER_NAME"@'%' IDENTIFIED BY "$MYSQL_STRUCTMAN_USER_PASSWORD";
GRANT ALL PRIVILEGES ON *.* TO "$MYSQL_STRUCTMAN_USER_NAME"@'%' IDENTIFIED BY "$MYSQL_STRUCTMAN_USER_PASSWORD" WITH GRANT OPTION;
FLUSH PRIVILEGES;
EOF

    echo "===>    MySQL user $MYSQL_STRUCTMAN_USER_NAME has been created with the password $MYSQL_STRUCTMAN_USER_PASSWORD    <==="

    # MySQL server shutdown
    /usr/bin/mysqladmin --defaults-file=/etc/mysql/debian.cnf shutdown

    if [[ $?==0 ]]; then
        echo "===>    MySQL server shutdown successful after the user creation step!    <==="
    else
        echo "===>    MySQL server could not be shutdown after user creation step. Aborting!    <==="
        exit 1
    fi
}

# Create a default database and import the default StructMAn SQL file if not exists
configure_database() {
    # Start MySQL server
    /usr/bin/mysqld_safe >/dev/null 2>&1 &
    sleep 60
    
    echo "===>    creating and importing the default database!    <==="
    gunzip < /usr/structman_library/sources/StructMAn_db/$STRUCTMAN_DB.sql.gz | mysql --user=root

    if [[ $? == 0 ]]; then
        echo "===>    MySQL StructMAn default database has been created and imported successfully!    <==="
    else
        echo "===>    MySQL StructMAn default database could not be created or imported. Aborting the container setup!   <==="
        exit 1
    fi

    sleep 10

    # MySQL server shutdown
    /usr/bin/mysqladmin --defaults-file=/etc/mysql/debian.cnf shutdown
    
    if [[ $?==0 ]]; then
        echo "===>    MySQL server shutdown successful after importing the default StructMAn database!    <==="
    else
        echo "===>    MySQL server could not be shutdown after importing the default StructMAn database. Aborting!    <==="
        exit 1
    fi
}

build_mmseqs_index() {
    (cd /usr/structman_library/sources/StructMAn/lib/base/blast_db/; mmseqs createindex /usr/structman_library/sources/StructMAn/lib/base/blast_db/pdbba_search_db_mmseqs2 /usr/structman_library/sources/StructMAn/lib/base/blast_db/tmp/)
    
    if [[ $?==0 ]]; then
        echo "===>    MMseqs2 index has been created for the database!    <==="
    else
        echo "===>    MMseqs2 index could not be created for the database!    <==="
        exit 1
    fi
}

# Initializing MySQL server configuration
configure_structman
create_mysql_dirs
configure_mysql
build_mmseqs_index

# Initialize MySQL server and create users
if [[ ! -d /var/lib/mysql/mysql ]]; then
    initialize_mysql_and_create_users
fi

# Configure default StructMAn database
if [[ ! -d /var/lib/mysql/$STRUCTMAN_DB_FOLDER_NAME ]]; then
    configure_database
fi

echo "***** Container setup and configuration is done, starting < $@ > on the container < $HOSTNAME > *****"

exec "$@"
