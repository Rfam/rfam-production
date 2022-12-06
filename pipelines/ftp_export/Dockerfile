FROM ubuntu:latest

ENV RNA=/rna

WORKDIR $RNA

RUN apt-get -y update

RUN apt-get install -y git

RUN cd $RNA && git clone https://github.com/Rfam/rfam-family-pipeline.git

ENV PATH="$RNA/rfam-family-pipeline:$PATH"

ENTRYPOINT ["/bin/bash"]