#!/bin/bash
# Description: StructMAn podman container setup script, that will setup a structman container with a default folder structure
# Based on: setup_structman_docker_container.sh 

usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-p <path>                 Provide a path to create the container"
        echo "Optional Parameters:"
        echo "-c <container_name>       Provide a name for the container"
        echo "Example:"
        echo "./setup_structman_podman_container.sh -p <path/to/create/the/container> -c <container_name>"
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

# podman bind mount paths must be absolute
path="$(realpath "$path")"

# Creates the default folder structure
mkdir -p $path/{mysql_lib,mysql_logs,structman/{input_data,results,resources}}

podman run -d --shm-size 8g -v "$path/mysql_lib/:/var/lib/mysql/:Z" -v "$path/mysql_logs/:/var/log/mysql/:Z" -v "$path/structman/input_data/:/structman/input_data/:Z" -v "$path/structman/results/:/structman/results/:Z" -v "$path/structman/resources/:/structman/resources/:Z" -e "MYSQL_STRUCTMAN_USER_NAME=structman" -e "MYSQL_STRUCTMAN_USER_PASSWORD=structman_rocks" --hostname "structman" --name "$container_name" structman
