#!/usr/bin/env python

import json
from io import StringIO
import tempfile
import typing as ty
from pathlib import Path

import click
from Bio import SeqIO

ENA_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=true"

import subprocess as sp

def generate_records(handle: ty.IO) -> ty.Iterable[SeqIO.SeqRecord]:
    for record in SeqIO.parse(handle, 'fasta'):
        record.id = record.id.split('|')[2]
        yield record


def fetch_ena(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    with tempfile.NamedTemporaryFile('w+') as tmp:
        sp.check_call(["wget", "-O", tmp.name, ENA_URL.format(accession=accession)])
        tmp.flush()

        info = sp.check_output(["file", '--mime-type', tmp.name])
        if b'application/gzip' in info:
            text = sp.check_output(['zcat', '-f', tmp.name])
            yield from generate_records(StringIO(text.decode()))
        elif b'text/plain' in info:
            yield from generate_records(tmp)
        else:
            raise ValueError(f"Do not know how handle {accession}")


def fetch(accessions: ty.Iterable[str]) -> ty.Iterable[SeqIO.SeqRecord]:
    for accession in accessions:
        yield from fetch_ena(accession)


@click.command()
@click.argument('ena_file', type=click.File('r'))
@click.argument('output', default='.', type=click.Path())
def main(ena_file, output):
    base = Path(output)
    for line in ena_file:
        data = json.loads(line)
        upi = data['upi']
        path = base / f"{upi}.fa"
        assert data['accession'] is None
        with path.open('w') as out:
            SeqIO.write(fetch(data['ids']), out, 'fasta')


if __name__ == '__main__':
    main()
