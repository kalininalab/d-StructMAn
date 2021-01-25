# StructMAn docker container

# Using Ubuntu 20.04 as the base image
FROM ubuntu:20.04

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
gzip \
python3 \
python3-dev \
python3-distutils \
python3-pip \
libboost-all-dev \
libzeep-dev \
libbz2-dev \
libz-dev \
autoconf \
automake \
autotools-dev \
rsync \
openbabel \
ncbi-blast+ \
python3-mysqldb \
mysql-client && \
rm -rf /var/lib/apt/lists/*

# Install StructMAn Python dependencies
RUN pip3 install numpy biopython matplotlib multiprocess pymysql python-igraph "pickle5>=0.0.10" psutil "chardet>=3.0.2,<4.0" ray

# Install and setup MMseqs2
RUN wget -O /opt/mmseqs-linux-sse41.tar.gz https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz && tar xvfz /opt/mmseqs-linux-sse41.tar.gz -C /opt/ && ln -s /opt/mmseqs/bin/mmseqs /usr/local/bin/

# Install xssp-3.0.7-mkdssp
RUN wget -O /opt/xssp-3.0.7.tar.gz https://github.com/cmbi/hssp/releases/download/3.0.7/xssp-3.0.7.tar.gz && tar xvzf /opt/xssp-3.0.7.tar.gz -C /opt/ && rm /opt/xssp-3.0.7.tar.gz
RUN cd /opt/xssp-3.0.7/ && ./autogen.sh && ./configure && make mkdssp && make install

# Add the StructMAn source
ADD ./structman_source /usr/structman_library/sources/

# Copy the entrypoint script to implement all the configuration
COPY docker-entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

# Default volumes to organize all the files under one roof and to allow backup
VOLUME ["/structman/input_data/", "/structman/results/", "/structman/resources/", "/var/log/mysql/", "/var/lib/mysql/"]

# Ports
EXPOSE 3306/tcp

# Initialization and setup
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

CMD ["mysqld"]
