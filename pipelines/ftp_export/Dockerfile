FROM gcc:10

RUN apt-get update && apt-get install -y git libapr1-dev libaprutil1-dev

ENV RNA /rna

RUN mkdir $RNA

WORKDIR $RNA

# Install Bio-Easel
RUN \
    cpan install Inline && \
    cpan install Inline::C && \
    cpan install Data::Printer && \
    cpan install Moose && \
    cpan install Config::General && \
    cpan install Catalyst::Utils && \
    cpan install DBIx::Class::Schema && \
    cpan install JSON && \
    cpan install DateTime && \
    cpan install File::Slurp && \
    cpan install DateTime::Format::MySQL && \
    cpan install XML::LibXML && \
    cpan install Mail::Mailer && \
    cpan install Data::UUID && \
    cpan install MooseX::NonMoose && \
    cpan install SVN::Client
RUN \
    git clone https://github.com/nawrockie/Bio-Easel.git && \
    cd Bio-Easel && \
    mkdir src && \
    cd src && \
    curl -k -L -o easel-Bio-Easel-0.09.zip https://github.com/EddyRivasLab/easel/archive/Bio-Easel-0.09.zip && \
    unzip easel-Bio-Easel-0.09.zip && \
    mv easel-Bio-Easel-0.09 easel && \
    rm easel-Bio-Easel-0.09.zip && \
    cd .. && \
    perl Makefile.PL; make; make test; make install

RUN git clone https://github.com/Rfam/rfam-family-pipeline.git

ENV PATH="$RNA/rfam-family-pipeline/Rfam/Scripts/jiffies:$PATH"
ENV BIOEASELDIR="$RNA/Bio-Easel/blib/lib:$RNA/Bio-Easel/blib/arch:$RNA/Bio-Easel:$RNA/Bio-Easel/lib"
ENV PERL5LIB="$BIOEASELDIR:rfam-family-pipeline/Rfam/Lib:rfam-family-pipeline/Rfam/Schemata:$PERL5LIB"
ENV RFAM_CONFIG="$RNA/rfam-family-pipeline/Rfam/Conf/rfam.conf"
ENTRYPOINT ["/bin/bash"]