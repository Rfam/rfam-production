FROM python:3.9

ENV RNA /usr/src/rfam

WORKDIR $RNA
COPY . .
COPY pyproject.toml $RNA

ENV PYTHONPATH=${PYTHONPATH}:${PWD}

RUN pip3 install poetry
RUN poetry config virtualenvs.create false
RUN poetry install --only main

ENTRYPOINT ["/bin/bash"]
