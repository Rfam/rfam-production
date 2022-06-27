#!/usr/bin/env python

import json
from io import StringIO
import tempfile
import typing as ty
from pathlib import Path
import subprocess as sp
import requests

import click
from Bio import SeqIO

ENA_FASTA_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=true"
ENA_EMBL_URL = "https://www.ebi.ac.uk/ena/browser/api/embl/{accession}?download=true"


class EnaIssue(Exception):
    pass


def generate_records(handle: ty.IO) -> ty.Iterable[SeqIO.SeqRecord]:
    seen = False
    for record in SeqIO.parse(handle, 'fasta'):
        record.id = record.id.split('|')[2]
        yield record
        seen = True
    if not seen:
        raise ValueError("Did not parse any sequences")


def fetch_ena_fasta(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    with tempfile.NamedTemporaryFile('w+') as tmp:
        try:
            sp.check_call(["wget", "--no-verbose", "-O", tmp.name, ENA_FASTA_URL.format(accession=accession)])
        except:
            raise EnaIssue
        tmp.flush()

        info = sp.check_output(["file", '--mime-type', tmp.name])
        if b'application/gzip' in info:
            with tempfile.NamedTemporaryFile('w+') as decomp:
                sp.run(['zcat', '-f', tmp.name], check=True, stdout=decomp)
                decomp.flush()
                yield from generate_records(decomp)
        elif b'text/plain' in info:
            yield from generate_records(tmp)
        else:
            raise ValueError(f"Do not know how handle {accession}")


def fetch_using_ena_embl(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    response = requests.get(ENA_EMBL_URL.format(accession=accession))
    response.raise_for_status()
    handle = StringIO(response.text)
    for line in handle:
        if not line.startswith('CON'):
            continue
        raw_contigs = line[3:].strip()
        yield from fetch_ena_fasta(raw_contigs)


def fetch_ena(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    try:
        yield from fetch_ena_fasta(accession)
    except EnaIssue:
        yield from fetch_using_ena_embl(accession)


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
