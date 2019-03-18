#!/bin/bash
# Description: StructMAn docker container setup script, that will setup a structman container with a default folder structure and a compose file.
# Author: Sanjay kumar Srikakulam (sanjaysrikakulam@gmail.com)
# Changelog:
# 13-07-2018 Init

usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-p <path>                 Provide a path to create the container"
        echo "-d <directory_name>       Provide a name for the container directory"
        echo "Optional Parameters:"
        echo "-c <container_name>       Provide a name for the container"
        echo "Example:"
        echo "./setup_docker_container.sh -p <> -d <> -c <>"
        exit 1
}

while getopts ":p:d:c:" i; do
        case "${i}" in
        p)
                path=$OPTARG
        ;;
        d)
                directory_name=$OPTARG
        ;;
        c)
                container_name=$OPTARG
        ;;
        esac
done

if [[ "$path" == "" ]] ; then
        usage
fi

if [[ "$directory_name" == "" ]] ; then
        usage
fi

if [[ "$container_name" == "" ]] ; then
        container_name="StructMAn"
fi

# Creates the default folder structure
if [[ -d "$path" ]] ; then
    mkdir -p $path/$directory_name/{mysql_lib,mysql_logs,structman/{input_data,results,mysql_custom_conf.d}}
fi

# Creates a default docker-compose file
echo "version: '2'
services:
    structman_db:
        image: docker.io/sanjaysrikakulam/structman
        restart: unless-stopped
        container_name: $container_name
        hostname: structman
        ports:
            - 3306
        cap_add:
            - NET_ADMIN
            - NET_RAW
        volumes:
            - ./structman/input_data/:/structman/input_data/:Z
            - ./structman/results/:/structman/results/:Z
            - ./structman/mysql_custom_conf.d/:/etc/mysql/mysql_custom_conf.d/:Z
            - ./mysql_lib/:/var/lib/mysql/:Z
            - ./mysql_logs/:/var/log/mysql/:Z
        environment:
            MYSQL_STRUCTMAN_USER_NAME: \"structman\"
            MYSQL_STRUCTMAN_USER_PASSWORD: \"structman_rocks\" " > $path/$directory_name/docker-compose.yml
