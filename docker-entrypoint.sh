#!/bin/bash
# Description: Docker entrypoint script to configure the StructMAn container. All ENV variables are set in the "docker-compose.yml" file
# Author: Sanjay kumar Srikakulam (sanjaysrikakulam@gmail.com)

set -euo pipefail

# Get the default StructMAn database name from the source file stored under /usr/structman_library/sources/StructMAn_db/
if [[ $(ls -L /usr/structman_library/sources/StructMAn_db/*.sql | wc -l) == 1 ]]; then
    for sql_file in /usr/structman_library/sources/StructMAn_db/*.sql; do
        STRUCTMAN_DB=$(basename $sql_file | cut -d "." -f 1)
    done
elif [[ $(ls -L /usr/structman_library/sources/StructMAn_db/*.sql | wc -l) == 0 ]]; then
    echo "===>    There does not seem to be any default StructMAn database file in the source directory (/usr/structman_library/sources/StructMAn_db/). Aborting the container setup!    <==="
    exit 1
fi

echo "***** Container configuration starting for the host $HOSTNAME *****"

# Remove old configuration files and make the new ones and change the file/folder permissions/ownership
configure_mysql() {
    rm -rf /etc/mysql/my.cnf &>/dev/null
    rm -rf /etc/mysql/mysql.conf.d/mysqld.cnf &>/dev/null
    
    # Setup the MySQL data directory
    mkdir -p /var/lib/mysql
    chmod -R 0777 /var/lib/mysql
    chown -R mysql:mysql /var/lib/mysql

    # Setup the MySQL run directory
    mkdir -p /var/run/mysqld
    chmod -R 0777 /var/run/mysqld
    chown -R mysql:mysql /var/run/mysqld
    rm -rf var/run/mysqld/mysqld.sock.lock

    # Setup the MySQL log directory
    mkdir -p /var/log/mysql
    chmod -R 0777 /var/log/mysql

    # Remove debian system's maintenance password
    sed 's/password = .*/password = /g' -i /etc/mysql/debian.cnf

    # Configure MySQL server
    if [[ ! -f /etc/mysql/my_structman.cnf ]]; then
        echo -e "# Include the configuration files from these directories
!includedir /etc/mysql/conf.d/
!includedir /etc/mysql/mysql_custom_conf.d/
!includedir /etc/mysql/mysql.conf.d/" > /etc/mysql/my_structman.cnf
    fi

    if [[ ! -f /etc/mysql/mysql.conf.d/structman_mysqld.cnf ]]; then
        echo -e "# MySQL Server configuration
[mysqld_safe]
socket=/var/run/mysqld/mysqld.sock
nice=0

[mysqld]
skip-host-cache
skip-name-resolve
datadir=/var/lib/mysql
tmpdir=/tmp
port=3306
bind-address=0.0.0.0
user=mysql
socket=/var/run/mysqld/mysqld.sock
pid-file=/var/run/mysqld/mysqld.pid
general_log_file=/var/log/mysql/query.log
slow_query_log_file=/var/log/mysql/slow.log
log-error=/var/log/mysql/error.log
basedir=/usr
init_connect='SET NAMES utf8'
character-set-server=utf8
collation-server=utf8_unicode_ci

[client]
default-character-set=utf8
socket=/var/run/mysqld/mysqld.sock

[mysql]
default-character-set=utf8
socket=/var/run/mysqld/mysqld.sock" > /etc/mysql/mysql.conf.d/structman_mysqld.cnf

    fi
    echo "===>	MySQL server configuration is done	<==="
}

# Initialize MySQL and ceate a new MySQL user for the application
initialize_mysql_and_create_users() {
    # Initialize MySQL data directory
    if [[ ! -d /var/lib/mysql/mysql ]]; then
        echo "===>    Initializing MySQL database server!	<==="
        mysqld --initialize-insecure --user=mysql > /dev/null 2>&1 &
        sleep 5

        # Start mysql server
        service mysql start > /dev/null 2>&1 &
        sleep 5

        # Create a debian-sys-maint user and a StructMAn db user
        echo "===>    Creating a debian system maintenance user and a MySQL StructMAn user  <==="
        mysql --user=root << EOF
CREATE USER IF NOT EXISTS 'debian-sys-maint'@'localhost' IDENTIFIED BY '';
GRANT ALL PRIVILEGES on *.* TO 'debian-sys-maint'@'localhost' IDENTIFIED BY '' WITH GRANT OPTION;

CREATE USER IF NOT EXISTS "$MYSQL_STRUCTMAN_USER_NAME"@'localhost' IDENTIFIED BY "$MYSQL_STRUCTMAN_USER_PASSWORD";
GRANT ALL PRIVILEGES ON *.* TO "$MYSQL_STRUCTMAN_USER_NAME"@'localhost' IDENTIFIED BY "$MYSQL_STRUCTMAN_USER_PASSWORD" WITH GRANT OPTION;
FLUSH PRIVILEGES;
EOF

        echo "===>    MySQL user $MYSQL_STRUCTMAN_USER_NAME has been created with the password $MYSQL_STRUCTMAN_USER_PASSWORD    <==="
    fi
}


# Create a default database and import the default StructMAn SQL file if not exists
configure_database() {
    if [[ ! -d /var/lib/mysql/$STRUCTMAN_DB ]]; then
        mysql --user=root < /usr/structman_library/sources/StructMAn_db/$STRUCTMAN_DB.sql

        if [[ $? == 0 ]]; then
            echo "===>    MySQL StructMAn default database has been created and imported successfully!    <==="
        else
            echo "===>    MySQL StructMAn default database could not be created or imported. Aborting the container setup!   <==="
            exit 1
        fi

        #Shutdown MySQL server
         service mysql stop &>/dev/null

         echo "===>    Container cleanup in progress!    <==="
         killall mysql &>/dev/null
    fi
}

# Adding "structman.py" to "/usr/local/bin" as a symlink to make it a command line utility
configure_structman() {
    if [[ -f /usr/structman_library/sources/StructMAn/structman.py ]]; then
        if [[ ! -e /usr/local/bin/structman.py ]]; then
            ln -s /usr/structman_library/sources/StructMAn/structman.py /usr/local/bin/
        fi
    else
            echo "===>    structman.py script could not be found!    <==="
    fi
}

# Initializing MySQL server configuration
configure_mysql
initialize_mysql_and_create_users
configure_database
configure_structman

echo "***** Container setup and configuration done, starting  $@ on $HOSTNAME *****"

exec "$@"
