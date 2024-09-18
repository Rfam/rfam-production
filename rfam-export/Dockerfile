FROM python:3.11.8

# Modified from: https://github.com/python-poetry/poetry/discussions/1879
RUN apt-get update && \
    apt-get install -y \
    curl \
    gzip

ENV RNA=/rna

WORKDIR $RNA

ENV PYTHONFAULTHANDLER=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=random \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=off \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    POETRY_VERSION=1.8.2 \
    POETRY_NO_INTERACTION=1 \
    PATH="$PATH:/root/.local/bin"

RUN pip install "poetry==$POETRY_VERSION"

COPY pyproject.toml poetry.lock ./
RUN poetry config virtualenvs.create false
RUN poetry install --without=dev --no-interaction --no-ansi

ENTRYPOINT ["/bin/bash"]
