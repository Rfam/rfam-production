FROM perl:5.30-slim
RUN apt-get update
RUN apt-get install -y curl build-essential

# SOFTWARE INSTALLATION
# Infernal installation
RUN mkdir -p /Rfam/software && cd /Rfam/software && \
curl -OL http://eddylab.org/infernal/infernal-1.1.2.tar.gz && \
tar -xvzf infernal-1.1.2.tar.gz && \
cd infernal-1.1.2 && \
./configure && \
make && \
make install && \
cd /Rfam/software/infernal-1.1.2/easel && \
make install && \
cd miniapps

RUN apt-get install -y git

RUN cpan -f install File::ShareDir::Install && \
cpan -f install Inline::C

# install Bio-Easel
RUN cd /Rfam && \
git clone https://github.com/nawrockie/Bio-Easel.git && \
cd Bio-Easel && \
mkdir src && cd src && \
git clone https://github.com/EddyRivasLab/easel.git easel && \
cd easel && git checkout tags/Bio-Easel-0.06 && rm -rf .git && \
cd /Rfam/Bio-Easel && perl Makefile.PL && \
make && \
make install
