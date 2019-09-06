FROM openjdk:8-jdk

# might need to install a particular version of perl
RUN apt-get update && apt-get install -y curl \
    gcc \
    git \
    tar \
    unzip \
    perl \
    wget \
    make \
    automake \
    curl \
    gzip \
    g++ \
    vim \
    python

# create the directories
RUN mkdir /software
RUN mkdir /rfam
RUN mkdir /rfam/local

RUN cd /rfam && git clone -b pipeline https://github.com/mb1069/rfam-production.git

# SOFTWARE INSTALLATION
# Infernal installation
RUN cd /software && \
curl -OL http://eddylab.org/infernal/infernal-1.1.2.tar.gz && \
tar -xvzf infernal-1.1.2.tar.gz && \
cd infernal-1.1.2 && \
./configure && \
make && \
make install && \
cd /software/infernal-1.1.2/easel && \
make install && \
cd miniapps

#COPY . /rfam/rfam-production
COPY config/rfam_local_cwl.py /rfam/rfam-production/config/rfam_local.py

ENV PYTHONPATH=/rfam/rfam-production
ENV PATH=/usr/bin:$PATH:/software/infernal-1.1.2/src:/software/infernal-1.1.2/src/miniapps
