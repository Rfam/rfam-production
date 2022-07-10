#!/usr/bin/env python

import os
import json
from io import StringIO
import tempfile
import typing as ty
from pathlib import Path
import subprocess as sp
import logging
from functools import lru_cache

import requests
import click
from Bio import SeqIO

ENA_FASTA_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=true"
ENA_EMBL_URL = "https://www.ebi.ac.uk/ena/browser/api/embl/{accession}?download=true"

LOGGER = logging.getLogger(__name__)


class EnaIssue(Exception):
    pass


def generate_records(handle: ty.Union[ty.IO, str]) -> ty.Iterable[SeqIO.SeqRecord]:
    seen_ids = set()
    for record in SeqIO.parse(handle, 'fasta'):
        if '|' in record.id:
            record.id = record.id.split('|')[2]
        if record.id in seen_ids:
            raise ValueError(f"Duplicate id {record.id}")
        yield record
        seen_ids.add(record.id)
    if not seen_ids:
        raise ValueError("Did not parse any sequences")


def fetch_ena_fasta(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    LOGGER.info("Attempt to get ena fasta for %s", accession)
    with tempfile.NamedTemporaryFile('w+', dir=os.curdir) as tmp:
        try:
            sp.check_call(["wget", "--no-verbose", "-O", tmp.name, ENA_FASTA_URL.format(accession=accession)])
        except:
            LOGGER.warn("Fasta lookup failed for %s", accession)
            raise EnaIssue
        tmp.flush()

        info = sp.check_output(["file", '--mime-type', tmp.name])
        if b'application/gzip' in info:
            LOGGER.debug("Decompressing file for %s", accession)
            with tempfile.NamedTemporaryFile('w+', dir=os.curdir) as decomp:
                sp.run(['zcat', '-f', tmp.name], check=True, stdout=decomp)
                decomp.flush()
                yield from generate_records(decomp.name)
        elif b'text/plain' in info:
            LOGGER.info("Parsing fasta file for %s", accession)
            yield from generate_records(tmp)
        else:
            raise ValueError(f"Do not know how handle {accession}")


@lru_cache
def fetch_contigs(accession: str) -> ty.Iterable[str]:
    LOGGER.info("Fetching EMBL formatted file for %s", accession)
    response = requests.get(ENA_EMBL_URL.format(accession=accession))
    response.raise_for_status()
    handle = StringIO(response.text)
    for line in handle:
        if not line.startswith('CON'):
            continue
        yield line[3:].strip()


def fetch_using_ena_embl(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    for raw_contig in fetch_contigs(accession):
        LOGGER.debug("Found contigs %s", raw_contig)
        yield from fetch_ena_fasta(raw_contig)


def fetch_accession(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    LOGGER.info("Fetching ENA data for %s", accession)
    try:
        yield from fetch_ena_fasta(accession)
    except EnaIssue:
        LOGGER.info("Using backup EMBL format based lookup for %s", accession)
        yield from fetch_using_ena_embl(accession)


def fetch_accessions(accessions: ty.Iterable[str]) -> ty.Iterable[SeqIO.SeqRecord]:
    for accession in accessions:
        LOGGER.info("Fetching %s from ENA", accession)
        yield from fetch_accession(accession)


@click.command()
@click.argument('ena_file', type=click.File('r'))
@click.argument('output', default='.', type=click.Path())
def main(ena_file, output):
    logging.basicConfig(level=logging.INFO)
    base = Path(output)
    for line in ena_file:
        data = json.loads(line)
        upi = data['upi']
        path = base / f"{upi}.fa"
        assert data['accession'] is None
        with path.open('w') as out:
            SeqIO.write(fetch_accessions(data['ids']), out, 'fasta')


if __name__ == '__main__':
    main()
