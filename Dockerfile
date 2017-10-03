FROM centos:6.6

RUN yum install -y \
    curl \
    gcc \
    git \
    httpd \
    httpd-devel \
    libaio \
    mysql-devel \
    nc.x86_64 \
    openssl \
    openssl-devel \
    tar \
    unzip \
    zlib-devel

RUN mkdir /rfam
RUN mkdir /rfam/local

ENV LOC /rfam/local

# Install Python
# NOTE: Python-2.7.11 and python-2.7.11 are DIFFERENT folders, the former contains the sources, the later - binaries
RUN \
    cd $LOC && \
    curl -OL http://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz && \
    tar -zxvf Python-2.7.11.tgz && \
    cd Python-2.7.11 && \
    PREFIX=$LOC/python-2.7.11/ && \
    export LD_RUN_PATH=$PREFIX/lib && \
    ./configure --prefix=$PREFIX  --enable-shared && \
    make && \
    make install && \
    cd $LOC && \
    rm -Rf Python-2.7.11 && \
    rm Python-2.7.11.tgz

# Install virtualenv
RUN \
    cd $LOC && \
    curl -OL  https://pypi.python.org/packages/source/v/virtualenv/virtualenv-15.0.0.tar.gz && \
    tar -zxvf virtualenv-15.0.0.tar.gz && \
    cd virtualenv-15.0.0 && \
    $LOC/python-2.7.11/bin/python setup.py install && \
    cd $LOC && \
    rm -Rf virtualenv-15.0.0.tar.gz && \
    rm -Rf virtualenv-15.0.0

# Create virtual environment
RUN \
    cd $LOC && \
    mkdir virtualenvs && \
    cd virtualenvs && \
    $LOC/python-2.7.11/bin/virtualenv rfam-production --python=$LOC/python-2.7.11/bin/python

# Define container environment variables
ENV PYTHONPATH /rfam/rfam-production

ENV DJANGO_SETTINGS_MODULE rfam_schemas.rfam_schemas.settings

# Install Python requirements
ADD requirements.txt $RFAM_PRODUCTION
RUN \
    source $LOC/virtualenvs/rfam-production/bin/activate && \
    pip install -r $RFAM_PRODUCTION/requirements.txt

COPY entrypoint.sh /entrypoint.sh
RUN chmod 700 entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
