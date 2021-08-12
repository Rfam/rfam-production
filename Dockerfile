FROM python:2.7

RUN apt-get update
RUN apt-get upgrade -y

# might need to install a particular version of perl
RUN apt-get install -y \
    automake \
    curl \
    g++ \
    gcc \
    git \
    gzip \
    make \
    perl \
    tar \
    unzip \
    vim \
    wget

# create the directories
RUN mkdir /software
RUN mkdir /rfam
RUN mkdir /rfam/local

ENV RFAM /rfam

WORKDIR $RFAM

RUN cd $RFAM && git clone -b pipeline https://github.com/mb1069/rfam-production.git

# SOFTWARE INSTALLATION
# Infernal installation
RUN cd /software && \
curl -OL http://eddylab.org/infernal/infernal-1.1.4.tar.gz && \
tar -xvzf infernal-1.1.4.tar.gz && \
cd infernal-1.1.4 && \
./configure && \
make && \
make install && \
cd /software/infernal-1.1.4/easel && \
make install && \
cd miniapps

#COPY . /rfam/rfam-production
#COPY config/rfam_local_cwl.py /rfam/rfam-production/config/rfam_local.py

# install requirements 
ENV RFAM_PRODUCTION "$RFAM/rfam-production"

ADD requirements.txt $RFAM_PRODUCTION/requirements.txt
RUN python -m pip install --upgrade pip
RUN python -m pip install -r $RFAM_PRODUCTION/requirements.txt

ENV PYTHONPATH=$RFAM_PRODUCTION
ENV PATH=/usr/bin:$PATH:/software/infernal-1.1.4/src:/software/infernal-1.1.4/src/miniapps
