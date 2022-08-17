# -*- coding: utf-8 -*-

"""
Copyright [2009-2022] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import csv
import json
import logging
from pathlib import Path
import typing as ty

import cattrs
import click
from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import download, metadata, ncbi, uniprot

LOGGER = logging.getLogger(__name__)


@click.group()
@click.option("--log-file", default=None, type=click.Path())
@click.option(
    "--log-level",
    default="info",
    type=click.Choice(
        ["critical", "error", "warn", "info", "debug"], case_sensitive=False
    ),
)
def cli(log_level="info", log_file=None):
    """
    Main entry point for updating rfamseq for Rfam.
    """
    logging.basicConfig(filename=log_file, level=getattr(logging, log_level.upper()))


@cli.command("build-metadata")
@click.argument("version")
@click.argument("ncbi-info", type=click.Path())
@click.argument("proteome-file", default="-", type=click.File("r"))
@click.argument("fasta-directory", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
def build_metadata_cmd(
    version: str,
    ncbi_info: str,
    proteome_file: ty.IO,
    fasta_directory,
    str,
    output: str,
):
    out = Path(output)
    fa_dir = Path(fasta_directory)
    with SqliteDict(ncbi_info, flag="r") as db:
        for line in proteome_file:
            raw = json.loads(line)
            proteome = cattrs.structure(raw, uniprot.ProteomeInfo)
            LOGGER.info("Building metadata for %s", proteome.upi)

            fetched = []
            fasta = fa_dir / f"{proteome.upi}.fa"
            for record in SeqIO.parse(fasta, "fasta"):
                fetched.append(metadata.FromFasta.from_record(record))

            metadata_out = out / f"{proteome.upi}.jsonl"
            genome = download.GenomeDownload.build(db, proteome)
            with metadata_out.open("w") as meta:
                info = metadata.build(version, proteome, genome.assembly_info, fetched)
                json.dump(cattrs.unstructure(info), meta)
                meta.write("\n")


@cli.command("download")
@click.argument("version")
@click.argument("ncbi-info", type=click.Path())
@click.argument("proteome-file", default="-", type=click.File("r"))
@click.argument("output", type=click.Path())
def download_cmd(version: str, ncbi_info: str, proteome_file: str, output: str):
    out = Path(output)
    with SqliteDict(ncbi_info, flag="r") as db:
        for line in proteome_file:
            raw = json.loads(line)
            proteome = cattrs.structure(raw, uniprot.ProteomeInfo)
            LOGGER.info("Downloading %s", proteome.upi)

            fetched = []
            fasta_out = out / f"{proteome.upi}.fa"
            genome = download.GenomeDownload.build(db, proteome)
            with fasta_out.open("w") as fasta:
                for sequence in genome.records(db):
                    fetched.append(metadata.FromFasta.from_record(sequence))
                    SeqIO.write(sequence, fasta, "fasta")

            metadata_out = out / f"{proteome.upi}.jsonl"
            with metadata_out.open("w") as meta:
                info = metadata.build(version, proteome, genome.assembly_info, fetched)
                json.dump(cattrs.unstructure(info), meta)
                meta.write("\n")


@cli.command("proteomes2genomes")
@click.option("--ignore", default=None, type=click.File("r"))
@click.argument("xml", type=click.Path())
@click.argument("output", default="-", type=click.File("w"))
def p2g_cmd(xml, output, ignore=None):
    to_skip = set()
    if ignore:
        to_skip.update(l.strip() for l in ignore)

    for proteome in uniprot.proteomes(Path(xml), to_skip):
        LOGGER.info("Working on %s", proteome.upi)
        data = cattrs.unstructure(proteome)
        output.write(json.dumps(data))
        output.write("\n")


@cli.command("parse-assembly-summary")
@click.argument("filename", default="-", type=click.File("r"))
@click.argument("output", type=click.Path())
def parse_assembly_info(filename, output):
    reader = csv.DictReader(filename, delimiter="\t")
    with SqliteDict(output, flag="r") as db:
        count = 0
        for index, row in enumerate(reader):
            count += 1
            summary = ncbi.NcbiAssemblySummary.from_ncbi_row(row)
            if summary.assembly_accession in db:
                raise ValueError("Multiple mappings to %s" % summary.assembly_accession)
            db[summary.assembly_accession] = summary
            if index % 100000 == 0:
                db.commit()
        db.commit()

        if count == 0:
            raise ValueError("Did not load any assemblies")

    with SqliteDict(output) as db:
        assert len(db) == count, "Did not load all assembiles"


if __name__ == "__main__":
    cli()
