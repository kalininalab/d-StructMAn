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
gzip \
python3.8 \
python3.8-dev \
python3.8-distutils \
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
python-mysqldb \
mysql-client && \
rm -rf /var/lib/apt/lists/* && \
rm -rf /var/lib/mysql

# Register the python version 3.8 in alternatives
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1

# Set python 3.8 as the default python3
RUN update-alternatives --set python3 /usr/bin/python3.8

# Upgrade pip to the latest version
RUN curl -s https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
python3 get-pip.py --force-reinstall && \
rm get-pip.py

# Clear pip cache
RUN rm -rf /root/.cache/pip

# Increase TCP backlog connection
RUN sysctl -w net.core.somaxconn=1024

# Install and setup MMseqs2
RUN wget -O /opt/mmseqs-linux-sse41.tar.gz https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz && tar xvfz /opt/mmseqs-linux-sse41.tar.gz -C /opt/ && ln -s /opt/mmseqs/bin/mmseqs /usr/local/bin/ && rm /opt/mmseqs-linux-sse41.tar.gz

# Custom XSSP
ADD ./tools/xssp /opt/
RUN cd /opt/xssp/ && ./autogen.sh && ./configure && make mkdssp && make install

# Add the StructMAn source
ADD ./structman_source /usr/structman_library/sources/

# Install StructMAn
RUN pip3 install /usr/structman_library/sources/

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
