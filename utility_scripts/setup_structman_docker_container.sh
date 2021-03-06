#!/bin/bash
# Description: StructMAn docker container setup script, that will setup a structman container with a default folder structure and a compose file.
# Author: Sanjay Kumar Srikakulam
# Changelog:
# 13-07-2018 Init

usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-p <path>                 Provide a path to create the container"
        echo "Optional Parameters:"
        echo "-c <container_name>       Provide a name for the container"
        echo "Example:"
        echo "./setup_structman_docker_container.sh -p <path/to/create/the/container> -c <container_name>"
        exit 1
}

while getopts ":p:c:" i; do
        case "${i}" in
        p)
                path=$OPTARG
        ;;
        c)
                container_name=$OPTARG
        ;;
        esac
done

if [[ "$path" == "" ]] ; then
        usage
fi

if [[ "$container_name" == "" ]] ; then
        container_name="structman"
fi

# Creates the default folder structure
mkdir -p $path/{mysql_lib,mysql_logs,structman/{input_data,results,resources}}


# Creates a default docker-compose file
echo "version: '3.5'
services:
    structman_db:
        shm_size: 8GB
        image: docker.io/alexgress/d-structman:latest
        container_name: $container_name
        hostname: structman-container
        ports:
            - 3306
        volumes:
            - ./structman/input_data/:/structman/input_data/:Z
            - ./structman/results/:/structman/results/:Z
            - ./structman/resources/:/structman/resources/:Z
            - ./mysql_lib/:/var/lib/mysql/:Z
            - ./mysql_logs/:/var/log/mysql/:Z
        environment:
            MYSQL_STRUCTMAN_USER_NAME: \"structman\"
            MYSQL_STRUCTMAN_USER_PASSWORD: \"structman_rocks\" " > $path/docker-compose.yml
