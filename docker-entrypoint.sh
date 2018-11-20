#!/bin/bash
# Description: Docker entrypoint script to configure the StructMAn container. All ENV variables are set in the "docker-compose.yml" file
# Author: Sanjay kumar Srikakulam (sanjaysrikakulam@gmail.com)

set -euo pipefail

echo "*** Container configuration starting for $HOSTNAME ***"

# Remove old configuration files and make the new ones and change the file/folder permissions/ownership
configure_mysql() {
    rm -rf /etc/mysql/my.cnf &>/dev/null
    rm -rf /etc/mysql/mysql.conf.d/mysqld.cnf &>/dev/null
    mkdir -p /var/lib/mysql /var/run/mysqld &>/dev/null
    chown -R mysql:mysql /var/lib/mysql /var/run/mysqld /var/log/mysql &>/dev/null
    chmod 777 /var/run/mysqld &>/dev/null

    echo -e "# Include the configuration files from these directories
!includedir /etc/mysql/conf.d/
!includedir /etc/mysql/mysql_custom_conf.d/
!includedir /etc/mysql/mysql.conf.d/" > /etc/mysql/my.cnf

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
socket=/var/run/mysqld/mysqld.sock" > /etc/mysql/mysql.conf.d/mysqld.cnf

    echo "MySQL server configuration is done"
}

# Create and configure new MySQL user and StructMAn database
configure_mysql_user_and_database() {
    # Change mysql root user's password
    mysqld_safe &
    sleep 5
    mysql --user=root << EOF
CREATE DATABASE IF NOT EXISTS $MYSQL_STRUCTMAN_DATABASE CHARACTER SET utf8 COLLATE utf8_general_ci;
USE $MYSQL_STRUCTMAN_DATABASE;
SOURCE /usr/structman_library/sources/StructMAn_db/struct_man_db.sql;
CREATE USER IF NOT EXISTS "$MYSQL_STRUCTMAN_USER_NAME"@localhost IDENTIFIED BY "$MYSQL_STRUCTMAN_USER_PASSWORD"; GRANT ALL ON *.* TO "$MYSQL_STRUCTMAN_USER_NAME"@localhost;
FLUSH PRIVILEGES;
EOF

    # Shutdown mysqld
    mysqladmin shutdown
	
	echo "Created MySQL user 'structman' and database 'structman' and its password is set to 'structman_rocks'"
}

# Initializing MySQL server configuration
configure_mysql
configure_mysql_user_and_database

echo "*** Container configuration done, starting  $@ on $HOSTNAME ***"

exec "$@"
