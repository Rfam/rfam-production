from python:3.11.8

# Modified from: https://github.com/python-poetry/poetry/discussions/1879
RUN apt-get update && \
    apt-get install -y \
    curl \
    gzip  \
    jq \
    rsync \
    wget

ENV RNA=/rna

WORKDIR $RNA

# Install Infernal/Easel
RUN \
    curl -OL http://eddylab.org/infernal/infernal-1.1.2.tar.gz && \
    tar -xvzf infernal-1.1.2.tar.gz && \
    rm infernal-1.1.2.tar.gz && \
    cd infernal-1.1.2 && \
    ./configure && \
    make && \
    make install && \
    cd easel && \
    make install

# Install seqkit, taken from:
# https://github.com/nextstrain/docker-base/blob/master/Dockerfile
RUN curl -L https://github.com/shenwei356/seqkit/releases/download/v2.2.0/seqkit_linux_amd64.tar.gz | tar xz --no-same-owner -C /usr/local/bin

ENV PYTHONFAULTHANDLER=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=random \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=off \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    POETRY_VERSION=1.8.1 \
    POETRY_NO_INTERACTION=1 \
    PATH="$PATH:/root/.local/bin"

RUN pip install "poetry==$POETRY_VERSION"

COPY pyproject.toml poetry.lock ./
RUN poetry config virtualenvs.create false
RUN poetry install --without=dev --no-interaction --no-ansi

ENTRYPOINT ["/bin/bash"]
