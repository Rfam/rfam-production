FROM gcc:10

ENV RNA=/rna

WORKDIR $RNA

RUN apt-get -y update

RUN apt-get install -y git


# Install Bio-Easel
RUN \
    cpan install Inline && \
    cpan install Inline::C
RUN \
    git clone https://github.com/nawrockie/Bio-Easel.git && \
    cd Bio-Easel && \
    git checkout e7fae0ab43fc47058183b71ff498ee0d0d2de6a7 && \
    mkdir src && \
    cd src && \
    curl -k -L -o easel-Bio-Easel-0.09.zip https://github.com/EddyRivasLab/easel/archive/Bio-Easel-0.09.zip && \
    unzip easel-Bio-Easel-0.09.zip && \
    mv easel-Bio-Easel-0.09 easel && \
    rm easel-Bio-Easel-0.09.zip && \
    cd .. && \
    perl Makefile.PL; make; make test; make install

RUN cd $RNA && git clone https://github.com/Rfam/rfam-family-pipeline.git

ENV PATH="$RNA/rfam-family-pipeline:$PATH"
ENV BIOEASELDIR="$RNA/Bio-Easel/blib/lib:$RNA/Bio-Easel/blib/arch:$RNA/Bio-Easel:$RNA/Bio-Easel/lib"
ENV PERL5LIB="$BIOEASELDIR:$PERL5LIB"

ENTRYPOINT ["/bin/bash"]