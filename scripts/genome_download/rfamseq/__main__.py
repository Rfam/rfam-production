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

import json
import logging
from pathlib import Path

import attrs
import cattrs
import click
from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import download, uniprot, metadata

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
                    fetched.append((sequence.id, sequence.description))
                    SeqIO.write(sequence, fasta, "fasta")

            metadata_out = out / f"{proteome.upi}.json"
            with metadata_out.open("w") as meta:
                info = metadata.build(version, proteome, genome.assembly_info, fetched)
                json.dump(attrs.asdict(info), meta)
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
        data = attrs.asdict(proteome)
        output.write(json.dumps(data))
        output.write("\n")


if __name__ == "__main__":
    cli()
