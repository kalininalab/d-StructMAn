# StructMAn

# What is StructMAn?

> The Structural Mutation Annotation (StructMAn) software provides the annotation of non-synonymous single-nucleotide polymorphisms (nsSNPs) in the context of the structural neighbourhood of the resulting amino acid variations in the protein. Its rationale is that if a mutation is located on an interaction interface between the protein and another protein, DNA, RNA or a small molecule, it is likely to interfere with this interaction.

# Why docker container for StructMAn

* A docker container will help a user to setup their own local installation of StructMAn 
* Alternatively, you can also use the webserver version of [StructMAn](http://structman.mpi-inf.mpg.de/)

# Requirements

* Install docker and docker-compose

*Example CentOS (execute as root user(denoted by #))*
```bash
# yum install docker
# yum install docker-compose
```
*NOTE: Replace "yum" with the respective package manager of your operating system*

# To get this image

The recommended way to get this StructMAn docker image is to pull the prebuilt image from the [Docker Hub Registry](https://hub.docker.com/r/sanjaysrikakulam/structman/).

```bash
$ docker pull docker.io/sanjaysrikakulam/structman:latest
```

## To use/run this image in a container

Once you pulled the image, use this utility script [setup_structman_docker_container.sh](https://github.com/sanjaysrikakulam/structman/blob/master/utility_scripts/setup_structman_docker_container.sh) to create base setup and configuration for your StructMAn to run succesfully with docker-compose

*Example CentOS (execute as root user(denoted by #))*
```bash
# ./setup_structman_docker_container.sh -p <> -d <> -c <>
```
```
Where: 
-p => A path to create the container
-d => A name for the container directory
-c => A name for the container (This is an optional parameter, default value is "StructMAn")
```
- This script will create a folder structure, docker-compose.yml file and will change the selinux context of them

```
└── docker-compose.yml
└── mysql_logs
└── structman
    ├── input_data
    ├── results
    ├── mysql_data_backup
    └── mysql_custom_conf.d
```
## Directory setup

- **docker-compose.yml:** This compose file contains default setup and configuration for your StructMAn instance to run successfully
- **mysql_logs:** This is where MySQL will write its error log
- **structman:** This is where we organize all the files, so that we can keep them under one roof
- **input_data:** A user has to upload their input files in this directory
- **reuslts:** StructMAn by default will write all the generated output to this directory
- **mysql_data_backup:** If this ENV "MYSQL_FREQUENT_DATABASE_DUMPS" in docker-compose.yml file is set to "true" then a cron job will be initiated to make MySQL structman database dumps frequently as a backup to this directory
- **mysql_custom_conf.d:** Before starting the container, a user can add any additional MySQL server configuration, if it wasn't already implemented/enabled
- **NOTE: All these volumes are bind mounted and for more information on [bind mounts](https://docs.docker.com/storage/bind-mounts/)**

## Docker Compose

- Once this utility script [setup_structman_docker_container.sh](https://github.com/sanjaysrikakulam/structman/blob/master/utility_scripts/setup_structman_docker_container.sh) has been executed and the above mentioned directory structure is created, a user can change from the current working directory to the container directory where docker-compose.yml file is created, then the user can start the container in the detach mode using the following command

*Example CentOS (execute as root user(denoted by #))*
```bash
# docker-compose up -d
```

# Issues

If you encountered a problem running this container, you can file an [issue](https://github.com/sanjaysrikakulam/structman/issues). For us to provide any form of support, be sure to include the following information in your issue:

- Host OS and version
- Docker version ('docker version')
- Output of `docker info` and 'docker logs <container-name>'
