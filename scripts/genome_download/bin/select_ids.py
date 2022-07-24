#!/usr/bin/env python

import os
import shutil
import json
from pathlib import Path
import typing as ty
import logging
import tempfile
import subprocess as sp

from Bio import SeqIO
from sqlitedict import SqliteDict

import click

import download_ena

LOGGER = logging.getLogger(__name__)

NCBI_SEQ_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"


def versionless(record: SeqIO.SeqRecord) -> str:
    return record.id.split(".", 1)[0]


def ncbi_fasta(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    with tempfile.NamedTemporaryFile("w+", dir=os.curdir) as tmp:
        sp.check_call(
            [
                "wget",
                "--no-verbose",
                "-O",
                tmp.name,
                NCBI_SEQ_URL.format(accession=accession),
            ]
        )
        yield from download_ena.generate_records(tmp.name)


def ncbi_lookup(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    try:
        LOGGER.info("Trying to fetch %s from NCBI", accession)
        yield from ncbi_fasta(accession)
    except:
        LOGGER.info("Trying to find contigs from %s", accession)
        contigs = list(download_ena.fetch_contigs(accession))
        for contig in contigs:
            LOGGER.info("Fetching contig %s from NCBI", contig)
            yield from ncbi_fasta(contig)


def select_ids(path: Path, allowed: ty.Set[str]) -> ty.Iterable[SeqIO.SeqRecord]:
    seen = set()
    with path.open("r") as raw:
        for record in SeqIO.parse(raw, "fasta"):
            LOGGER.info("Checking if %s is allowed", record.id)
            for key in [record.id, versionless(record)]:
                LOGGER.debug("Trying key %s", key)
                if key not in allowed:
                    LOGGER.debug("Skipping %s", key)
                    continue
                if key in seen:
                    raise ValueError(f"Duplicate ids in fasta file at {path}")
                LOGGER.debug("Keeping %s", key)
                seen.add(key)
                yield record
                break

    if len(seen) != len(allowed):
        missing = allowed - seen
        LOGGER.info("Missing %i ids: %s, using lookup", len(missing), missing)
        for accession in missing:
            LOGGER.info("Looking up %s with ENA", accession)
            try:
                yield from download_ena.fetch_accession(accession)
            except:
                LOGGER.info("Failed to fetch %s from ENA, trying NCBI", accession)
                yield from ncbi_lookup(accession)


def load_ignore(handle: ty.IO) -> ty.Set[str]:
    ignore = set()
    for line in handle:
        data = json.loads(line)
        ignore.add(data["upi"])
    return ignore


@click.command()
@click.option("--ignore-file", type=click.Path())
@click.argument("gca-file", type=click.File("r"))
@click.argument("directory", type=click.Path())
@click.argument("mapping_file", type=click.Path())
@click.argument("output", default=".", type=click.Path())
def main(
    gca_file: ty.IO,
    directory: str,
    mapping_file: str,
    output: str,
    ignore_file: ty.Optional[str] = None,
):
    logging.basicConfig(level=logging.INFO)
    base = Path(directory)
    out = Path(output)
    ignore = set()
    if ignore_file and Path(ignore_file).exists():
        with open(ignore_file, "r") as raw:
            ignore = load_ignore(raw)

    for line in gca_file:
        row = json.loads(line)
        upi = row["upi"]
        LOGGER.info("Checking %s", upi)
        if upi in ignore:
            LOGGER.debug("Skipping %s", upi)
            continue

        path = base / f"{upi}.fa"
        if not path.exists():
            raise ValueError(f"Failed to find genome for {upi} at {path}")

        final_path = out / f"{upi}.fa"
        if row["ids"] == ["*"]:
            LOGGER.info("Using all ids for %s", upi)
            shutil.copyfile(path, final_path)
            return

        allowed = {id for id in row["ids"]}
        if row["accession"].split(".", 1)[0] in allowed:
            LOGGER.info("Versionless match, using all ids for %s", upi)
            shutil.copyfile(path, final_path)
            return

        with open(final_path, "w") as f:
            SeqIO.write(select_ids(path, allowed), f, "fasta")


if __name__ == "__main__":
    main()
