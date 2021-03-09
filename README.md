# StructMAn

## What is StructMAn?

> The Structural Mutation Annotation (StructMAn) software provides annotation of non-synonymous single-nucleotide variants (nsSNVs) in the context of the structural properties of the resulting amino acid changes in the corresponding proteins. Its rationale is that if a mutation is located on an interaction interface between a protein and another protein, DNA, RNA or a small molecule, it is likely to interfere with this interaction, and mutations location in the protein core are likely to influence its stability. StructMAn was first published in [Nucl Acids Res, 2016](https://academic.oup.com/nar/article-abstract/44/W1/W463/2499349) where we showed that such structural annotation correlate well with established tools for predicting damaging effect of missense mutations. Later we have [applied StructMAn to a large collection of disease-associated and neutral mutations](https://www.nature.com/articles/oncsis201779) and discovered distinct trends of their spatial distribution.

## Why containers for StructMAn

* A container will help a user to setup their own local installation of StructMAn
* Here we present [docker](https://www.docker.com/) and [podman](https://podman.io/) versions of containers for StructMAn

# Docker version of StructMAn

## Requirements

**Before you start, please make sure you have internet connection to use SructMAn**

**NOTE: All the example code/commands here are based on CentOS 7 (everything is executed as root user(denoted by #))**

* Install [docker engine](https://docs.docker.com/install/) and docker-compose

*NOTE: Enable/install epel repo if by default it is not available in your system (for CentOS it is not available by default).*

```bash
# yum install epel-release
```

```bash
# yum install docker-compose
```
**NOTE 1: Replace "yum install" with the respective command of your operating system to install**

**NOTE 2: After installing docker, start and enable docker by doing the following (this is the case for CentOS only, check the [docker documentation](https://docs.docker.com/install/) if your OS requires you to do so)**

**NOTE 3: d-StructMAn should be provided with at least 4Gb of RAM. MAC users have to be aware of the memory given to the container ([see](https://docs.docker.com/docker-for-mac/) for more information).**

```bash
# systemctl start docker
# systemctl enable docker
```

## To use this image

Download this utility script [setup_structman_docker_container.sh](https://raw.githubusercontent.com/kalininalab/d-StructMAn/master/utility_scripts/setup_structman_docker_container.sh) to your local machine in order to create the base setup and configuration for your StructMAn to run succesfully using docker-compose client

```bash

- Download the script using cURL

# curl -O https://raw.githubusercontent.com/kalininalab/d-StructMAn/master/utility_scripts/setup_structman_docker_container.sh

- Set the execution flag to the script once it is downloaded from github

#  chmod +x <path>/setup_structman_docker_container.sh

- Now use the script to create default setup and configuration for StructMAn (Please change the values in the angular brackets "<>"),
make sure you always use the latest version of the sript, as there might be some updates from time to time.

# ./setup_structman_docker_container.sh -p <path/to/create/the/container_directory_name> -c <container_name>

Where: 
-p => A path to create the container direcotry
-c => A name for the container (This is an optional parameter, and its default value is "structman")
```
- This script will create the following folder structure along with the docker-compose.yml file in the path you provided above with "-p" option

```
<container_directory>
    └── docker-compose.yml
    └── mysql_lib
    └── mysql_logs
    └── structman
        ├── input_data
        └── results
        └── resources
        
```
## Directory setup

- **docker-compose.yml:** This compose file contains default setup and configuration for your StructMAn instance to run successfully
- **mysql_lib:** MySQL will use this folder as its data directory
- **mysql_logs:** MySQL will use this folder to store error logs
- **structman:** This is where we organize all the files, so that we can keep them under one roof
- **input_data:** A user has to upload their input files into this directory
- **results:** StructMAn by default will write all the generated output to this directory
- **resources:** Useful resources such as local PDB or other executables required by StructMAn should be put here before starting the container. For more on resources read [here](https://github.com/kalininalab/d-StructMAn/wiki/Local-resources)
- **NOTE 1: All these volumes are bind mounted and for more information on bind mounts read [here](https://docs.docker.com/storage/bind-mounts/)**
- **NOTE 2: Since the mysql_lib directory is bind mounted, every time when the container starts with either new or old image the database will not be reset, therefore at any time if you want [to have a clean setup](https://github.com/kalininalab/d-StructMAn#to-reset-structman-database-to-have-a-clean-structman-container), please delete all the files located under mysql_lib directory that was created using the utility script before starting the container. Similarly, delete all the files located under input_data and results directories as well.**

## Docker Compose

- After the script has created the above mentioned directory structure, you can change to the container directory where docker-compose.yml file is created, then you can start the container in the detach mode using the following command

```bash
# cd <container_directory>

# docker-compose up -d
```
- Wait for a while (minimum 5 minutes) before you start running the commands. **To check if the container has properly started or not, check the logs**

```bash
# docker-compose logs
```

- If you see the line ''***** Container setup and configuration is done, starting < mysqld > on the container <structman> '' at the end of the output from the above command, then it means the container has started successfully and that it is fully functional now.

## How to access the StructMAn tool running in the container
* To list the running docker containers
 
 ```bash
 # docker ps
```
Use the following command to run StructMAn, make sure you have some input file located under <container_directory>/structman/input_data

 ```bash
 # docker exec -it <container_name> structman.py -i /structman/input_data/<input_file_name>
```
For ease of use, lets setup an alias for the `docker exec` command either in the current terminal or by adding it to the end of the `.bashrc` file like below,
 ```bash
 # alias structman='docker exec -it <container_name> structman.py'
```
Once an alias is set, you can easily access the structman.py from your current terminal simply by,
 ```bash
 # structman -i /structman/input_data/<input_file_name>
```

**Refer this [tutorial](https://github.com/kalininalab/d-StructMAn/wiki/Tutorial) for more details on how to use StructMAn**

**More information about d-StructMAn can be found in its [wiki](https://github.com/kalininalab/d-StructMAn/wiki)** 

**NOTE: If you do not provide an input file using the "-i" option StructMAn will by default use all input files found in the <container_directory>/structman/input_data directory and store the output in the <container_directory>/structman/results directory. These paths are bind mounted to the container.**

## To update the existing StructMAn docker image

* If there is a new image of StructMAn available, one can easily update their container to use this latest image by doing the following

 ```bash
 - Stop the current running container by providing the path of the docker-compose.yml to the docker-compose client
 
 # docker-compose -f <docker-compose.yml> stop

 - Remove the container that we stopped just now
 
 # docker-compose -f <docker-compose.yml> rm -f
 
 - Delete the old image
 
 # docker rmi <image_id>

 - Finally pull the latest image and start the container in the detach mode
 
 # docker-compose -f <docker-compose.yml> up -d
```

## To reset StructMAn database/ To have a clean StructMAn container 

* If you encounter some issues with the database/container while executing structman.py with some input or you simply want a clean setup of the container, you can do the following

```bash
 - Stop the current running container by providing the path of the docker-compose.yml to the docker-compose client
 
 # docker-compose -f <docker-compose.yml> stop

 - Remove the container that we stopped just now
 
 # docker-compose -f <docker-compose.yml> rm -f

 - Remove all the files located under mysql_lib directory (Remember that this directory was created with the help of the utility script and that this direcotry is bind mounted into the container. Removing everything inside this folder will reset the database and reinitializes the mysql server setup)
 
 # rm -r <path_to_the_mysql_lib_directory> 

 - Finally start the container in the detach mode
 
 # docker-compose -f <docker-compose.yml> up -d
```

## To have more than one StructMAn container

* If you want to have more than one StructMAn container

1) Use the utility script like mentioned above in [To use/run this image](https://github.com/kalininalab/d-StructMAn#to-userun-this-image) section to create the new directory structure for your n'th instance of the StructMAn by exclusively changig the -p and -c option to some other value than the ones you used for your previous instances

# Podman version of StructMAn (for users with no root access)

## Requirements

**NOTE: All the example code/commands here are based on Fedora 30**

* Install [podman](https://github.com/containers/libpod/blob/master/install.md) **only for this step you require to be root**

```bash
# dnf install podman
```
**NOTE: Replace "dnf install" with the respective command of your operating system to install**

## To get this image

Pull the image, just like with docker

```bash
$ podman pull "docker.io/alexgress/d-structman:latest"
```

## To use/run this image

Like in the docker version we provide a setup script that creates the required directory structure and sets up the container with the right bind mounts. Download [setup_structman_podman_container.sh](https://raw.githubusercontent.com/kalininalab/d-StructMAn/master/utility_scripts/setup_structman_podman_container.sh) instead of `setup_structman_docker_container.sh` and follow the [instruction given for the Docker version of the script](#to-use-this-image).


Now Wait for a while (minimum 5 minutes) before you start running the commands. To check if the container has properly started or not, check the logs

```bash
# podman logs <container_name>
```

- If you see the line ''***** Container setup and configuration is done, starting < mysqld > on the container <structman> '' at the end of the output from the above command, then it means the container has started successfully and that it is fully functional now.

**That's it, now you can enjoy structman even without bothering your Admin** (*Remember with great power comes great responsibility, so!!!*)

**Almost all docker commands work for podman, just replace docker with podman and if you are not sure use this [commands manual](https://github.com/containers/podman/blob/master/commands-demo.md)**

## How to access the StructMAn tool running in the container
* To list the running podman containers
 
 ```bash
 $ podman ps
```
Use the following command to run StructMAn, make sure you have some input file located under <container_directory>/structman/input_data

 ```bash
 $ podman exec -it <container_name> structman.py -i /structman/input_data/<input_file_name>
```
For ease of use, lets setup an alias for the `podman exec` command either in the current terminal or by adding it to the end of the `.bashrc` file like below,
 ```bash
 # alias structman='podman exec -it <container_name> structman.py'
```
Once an alias is set, you can easily access the structman.py from your current terminal simply by,
 ```bash
 # structman -i /structman/input_data/<input_file_name>
```

**Refer this [tutorial](https://github.com/kalininalab/d-StructMAn/wiki/Tutorial) for more details on how to use StructMAn**

**More information about d-StructMAn can be found in its [wiki](https://github.com/kalininalab/d-StructMAn/wiki)** 

**NOTE: If you do not provide an input file using the "-i" option StructMAn will by default use all input files found in the <container_directory>/structman/input_data directory and store the output in the <container_directory>/structman/results directory. These paths are bind mounted to the container.**

## How to cite container tools

If you find these tools useful, please cite:

```Gress A, Srikakulam SK, Keller S, Kalinina OV. d-StructMAn: containerized structural annotation on the scale from genetic variants to whole proteomes. submitted.```

# Useful commands

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

If you encountered a problem running this container, you can file an [issue](https://github.com/kalininalab/d-StructMAn/issues). For us to provide any form of support, be sure to include the following information in your issue:

- Host OS and version
- Docker version ('docker version')
- Output of `docker info` and `docker logs <container_name>`
