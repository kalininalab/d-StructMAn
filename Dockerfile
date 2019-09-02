# StructMAn docker container

# Using Ubuntu 18.04 as the base image
FROM ubuntu:18.04

# Meta-data
LABEL maintainer="Alexander Gress (agress@mpi-inf.mpg.de)" \
      description="A docker container for the Structural Mutation Annotation (StructMAn) software which provides \
          the annotation of non-synonymous single-nucleotide polymorphisms (nsSNPs) in the \
          context of the structural neighbourhood of the resulting amino acid variations in \
          the protein. Its rationale is that if a mutation is located on an interaction \
          interface between the protein and another protein, DNA, RNA or a small molecule, \
          it is likely to interfere with this interaction."

# Install and update the required dependencies for StructMAn
RUN apt-get update \
&& DEBIAN_FRONTEND=noninteractive apt-get install -y mysql-server \
wget \
curl \
vim \
less \
gcc \
python-pip \
openbabel \
dssp \
ncbi-blast+ \
python-mysqldb \
mysql-client && \
rm -rf /var/lib/apt/lists/* && \
rm -rf /var/lib/mysql

RUN pip install numpy biopython matplotlib multiprocess
RUN wget -O /opt/mmseqs-linux-sse41.tar.gz https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz; tar xvfz /opt/mmseqs-linux-sse41.tar.gz -C /opt/; echo 'export PATH=$PATH:/opt/mmseqs/bin/"' >> ~/.profile

# Adding the StructMAn source
ADD ./structman_source /usr/structman_library/sources/

# Copying the entrypoint script to implement all the configuration
COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

# Default volumes to organize all the files under one roof and to allow backup
VOLUME ["/structman/input_data/", "/structman/results/", "/var/log/mysql/", "/var/lib/mysql/"]

# Ports
EXPOSE 3306/tcp

# Initialization and setup
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

CMD ["mysqld"]
