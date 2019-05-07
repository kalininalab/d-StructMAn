# StructMAn

# What is StructMAn?

> The Structural Mutation Annotation (StructMAn) software provides the annotation of non-synonymous single-nucleotide polymorphisms (nsSNPs) in the context of the structural neighbourhood of the resulting amino acid variations in the protein. Its rationale is that if a mutation is located on an interaction interface between the protein and another protein, DNA, RNA or a small molecule, it is likely to interfere with this interaction.

# Why docker container for StructMAn

* A docker container will help a user to setup their own local installation of StructMAn 
* Alternatively, you can also use the webserver version of [StructMAn](http://structman.mpi-inf.mpg.de/)

# Requirements

* Install docker and docker-compose (https://docs.docker.com/install/)

*NOTE: Enable/install epel repo if by default it is not available in your system*

*Example CentOS (execute as root user(denoted by #))*
```bash
# yum install epel-release
```

*Example CentOS (execute as root user(denoted by #))*
```bash
# yum install docker
# yum install docker-compose
```
**NOTE 1: Replace "yum" with the respective package manager of your operating system**

**NOTE 2: After installing docker, start and enable docker by doing systemctl start docker and systemctl enable docker (This is the case for CentOS only, check the docker documentation if your OS requires you to do so)**



# To get this image

The recommended way to get this StructMAn docker image is to pull the prebuilt image from the [Docker Hub Registry](https://hub.docker.com/r/sanjaysrikakulam/structman/).

```bash
# docker pull docker.io/sanjaysrikakulam/structman:latest
```

## To use/run this image

Once you pulled the image, use this utility script [setup_structman_docker_container.sh](https://github.com/sanjaysrikakulam/structman/blob/master/utility_scripts/setup_structman_docker_container.sh) by downloading the script to your local machine in order to create the base setup and configuration for your StructMAn to run succesfully using docker-compose client

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
└── mysql_lib
└── mysql_logs
└── structman
    ├── input_data
    └── results
```
## Directory setup

- **docker-compose.yml:** This compose file contains default setup and configuration for your StructMAn instance to run successfully
- **mysql_lib:** MySQL will use this folder as its data directory
- **mysql_logs:** MySQL will use this folder to store error logs
- **structman:** This is where we organize all the files, so that we can keep them under one roof
- **input_data:** A user has to upload their input files into this directory
- **reuslts:** StructMAn by default will write all the generated output to this directory
- **NOTE 1: All these volumes are bind mounted and for more information on [bind mounts](https://docs.docker.com/storage/bind-mounts/)**
- **NOTE 2: Since the mysql_lib directory is bind mounted, every time when the container starts with either new or old image the database will not be reset, therefore at any time if you want to have a clean setup, please delete all the files located under mysql_lib directory that was created using the utility script. Similarly, delete all the files located under input_data and results directories as well.**


## Docker Compose

- Once this utility script [setup_structman_docker_container.sh](https://github.com/sanjaysrikakulam/structman/blob/master/utility_scripts/setup_structman_docker_container.sh) has been executed and the above mentioned directory structure is created, a user can change from the current working directory to the container directory where docker-compose.yml file is created, then the user can start the container in the detach mode using the following command

*Example CentOS (execute as root user(denoted by #))*
```bash
# docker-compose up -d
```
# How to access the structman tool running in the container
* To list the running docker containers
 
 ```bash
 # docker ps
```

* Type 1: To run StructMAn from within the running container

 ```bash
 # docker exec -it <container_name> /bin/bash
```
*This will allow you to access the shell of the running container*

* Accessing the StructMAn application

 ```bash
 # structman.py -i </structman/input_data/<input_file_name>> -o </structman/results/> 
```

* Type 2: To access StructMAn from your localhost itself

 ```bash
 # docker exec -it StructMAn structman.py -i </structman/input_data/<input_file_name>> -o </structman/results/>
```
*This command will run the application (StructMAn) inside the container from your localhost*

**NOTE: By default StructMAn will look for the input here </structman/input_data/> and stores the output here </structman/results/>. These paths are located inside the docker container and these folders are bind mounted to your localhost**

## To update the existing StructMAn docker image

* If there is a new image of StructMAn available, one can easily update their container to use this latest image by doing the following


 ```bash
 - Stop the current running container by providing the path of the docker-compose.yml to the docker-compose client
 # docker-compose -f <docker-compose.yml> stop

 - Remove the container that we stopped just now
 # docker-compose -f <docker-compose.yml> rm -f

 - Pull the latest image of StructMAn from Docker Hub Registry (The docker compose.yml file already contains the image location which is pointing to the docker hub's registry, so just a pull command will download the latest image from the registry)
 # docker-compose -f <docker-compose.yml> pull

 - Finally start the container in the detach mode
 # docker-compose -f <docker-compose.yml> up -d
```


## Useful docker commands

* To start a container

```bash
 # docker-compose up -d
 ```
 
 * To start one or more stopped containers
 
 ```bash
 # docker start <container_name>
 ```
 
 * To start a container by specifying a path to a compose file
 
  ```bash
 # docker-compose -f <path_to_your_compose_file/docker-compose.yml> up -d
  ```
  
 * To stop a container

```bash
 # docker stop <container_name>
 ```

* To read the logs of a running docker container

```bash
# docker logs <container_name>
```

* To access the shell of a running docker container

```bash
# docker exec -it <container_name> /bin/bash
```

* To remove dangling images (images with <none> tag or repo name)

```bash
# docker images -q -f dangling=true | xargs --no-run-if-empty docker rmi
```

* To list all the docker images

```bash
# docker images -a
```

* To simply remove an image

```bash
# docker rmi <image_name or image_id>
```

* To list all exited containers

```bash
# docker ps -a -f status=exited
```

* To remove all exited containers

```bash
# docker rm $(docker ps -a -f status=exited -q)
```

* To remove a single container

```bash
# docker rm <container_name>
```

# For issues

If you encountered a problem running this container, you can file an [issue](https://github.com/sanjaysrikakulam/structman/issues). For us to provide any form of support, be sure to include the following information in your issue:

- Host OS and version
- Docker version ('docker version')
- Output of `docker info` and `docker logs <container_name>`
