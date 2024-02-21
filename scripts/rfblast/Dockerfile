FROM python:3.9

ENV RNA /usr/src/rfam
ENV INFERNAL /usr/src/infernal

RUN mkdir $INFERNAL

WORKDIR $RNA
COPY . .
COPY pyproject.toml $RNA

ENV PYTHONPATH=${PYTHONPATH}:${PWD} 

# Install Infernal
RUN \
    cd $INFERNAL && \
    curl -OL http://eddylab.org/infernal/infernal-1.1.4.tar.gz && \
    tar -xvzf infernal-1.1.4.tar.gz && \
    cd infernal-1.1.4 && \
    ./configure --prefix=$INFERNAL/infernal-1.1.4 && \
    make && \
    make install && \
    cd easel && \
    make install && \
    cd $INFERNAL && \
    rm infernal-1.1.4.tar.gz

ENV PATH="$INFERNAL/infernal-1.1.4/bin:$RNA:$PATH"

RUN pip3 install poetry
RUN poetry config virtualenvs.create false
RUN poetry install --no-dev

ENTRYPOINT ["/bin/bash"]
