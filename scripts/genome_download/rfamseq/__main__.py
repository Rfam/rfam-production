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
import sys
import typing as ty
from pathlib import Path

import cattrs
import click
from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import download, metadata, ncbi, uniprot
from rfamseq.utils import serialize

LOGGER = logging.getLogger(__name__)


@click.group()
@click.option(
    "--log-file", default=None, type=click.Path(), help="File to log to, default stdout"
)
@click.option(
    "--log-level",
    default="info",
    type=click.Choice(
        ["critical", "error", "warn", "info", "debug"], case_sensitive=False
    ),
    help="Set the log level",
)
def cli(log_level="info", log_file=None):
    """
    Main entry point for updating rfamseq for Rfam.
    """
    logging.basicConfig(filename=log_file, level=getattr(logging, log_level.upper()))


@cli.command("build-metadata")
@click.argument("version")
@click.argument("ncbi-info", type=click.Path())
@click.argument(
    "proteome-file",
    type=click.File("r"),
)
@click.argument(
    "fasta-directory",
    type=click.Path(exists=True),
)
@click.argument(
    "output",
    default=".",
    type=click.Path(),
)
def build_metadata_cmd(
    version: str,
    ncbi_info: str,
    proteome_file: ty.IO,
    fasta_directory: str,
    output: str,
):
    """
    Build metadata for downloaded fasta files. This is normally done when
    downloading files, but we are still tweaking the metadata generation so this
    will build metadata from fetched fasta files. This assumes the fasta files
    are named with UP* and stored in one directory.

    Arguments:
      VERSION          Version of the rfamseq database, eg 15.0
      NCBI-INFO        Path to file produced by parse-assembly-summary
      PROTEOME-FILE    Path to the file produced by: proteomes2genomes
      FASTA-DIRECTORY  Path to a directory with all fasta files produced by download
      OUTPUT           Path to write the up*.jsonl files to
    """
    out = Path(output)
    fa_dir = Path(fasta_directory)
    failed = False
    with SqliteDict(ncbi_info, flag="r") as db:
        for line in proteome_file:
            raw = json.loads(line)
            proteome = cattrs.structure(raw, uniprot.ProteomeInfo)
            LOGGER.info("Building metadata for %s", proteome.upi)

            fetched = []
            fasta = fa_dir / f"{proteome.upi}.fa"
            if not fasta.exists():
                LOGGER.error(
                    "Missing fasta file %s for proteome %s", fasta, proteome.upi
                )
                failed = True
                continue

            for record in SeqIO.parse(fasta, "fasta"):
                fetched.append(metadata.FromFasta.from_record(record))

            metadata_out = out / f"{proteome.upi}.jsonl"
            genome = download.GenomeDownload.build(db, proteome)
            with metadata_out.open("w") as meta:
                info = metadata.Metadata.build(
                    version, proteome, genome.assembly_info, fetched
                )
                json.dump(cattrs.unstructure(info), meta, default=serialize)
                meta.write("\n")

    if failed:
        sys.exit(1)


@cli.command("download")
@click.argument("version")
@click.argument(
    "ncbi-info",
    type=click.Path(),
)
@click.argument(
    "proteome-file",
    default="-",
    type=click.File("r"),
)
@click.argument(
    "output",
    default=".",
    type=click.Path(),
)
def download_cmd(version: str, ncbi_info: str, proteome_file: str, output: str):
    """
    Download the genomes specified in proteome-file. The file is the result of
    proteomes2genomes or a chunk of that file. This will download the genome and
    limit it to the specified components. This will produce a UP*.fa and
    UP*.jsonl file for all proteomes in the file. The .jsonl contains metadata
    and .fa contains the specified metadata.

    \b
    Arguments:
      VERSION        Version of the rfamseq database, eg 15.0
      NCBI-INFO      Path to file produced by parse-assembly-summary
      PROTEOME-FILE  Path to the file produced by: proteomes2genomes
      OUTPUT         Path to write the sequence fasta and metadata jsonl files to
    """
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
                info = metadata.Metadata.build(
                    version, proteome, genome.assembly_info, fetched
                )
                json.dump(cattrs.unstructure(info), meta, default=serialize)
                meta.write("\n")


@cli.command("proteomes2genomes")
@click.option(
    "--ignore",
    default=None,
    type=click.File("r"),
    help="A file with one UPI per line of UPIs to ignore when parsing.",
)
@click.argument(
    "xml",
    type=click.Path(),
)
@click.argument(
    "output",
    default="-",
    type=click.File("w"),
)
def p2g_cmd(xml, output, ignore=None):
    """
    Parse the XML file provided by uniprot to a jsonl file. The jsonl file
    contains one JSON encoded object per line where the object contains
    everything needed to download a genome from NCBI/ENA.

    \b
    Arguments:
      XML     Path to the proteome.xml file fetched from uniprot.
      OUTPUT  Path to write a summary jsonl to
    """
    to_skip = set()
    if ignore:
        to_skip.update(l.strip() for l in ignore)

    for proteome in uniprot.proteomes(Path(xml), to_skip):
        LOGGER.info("Working on %s", proteome.upi)
        data = cattrs.unstructure(proteome)
        output.write(json.dumps(data))
        output.write("\n")


@cli.command("parse-assembly-summary")
@click.argument(
    "filename",
    default="-",
    type=click.File("r"),
)
@click.argument("output", type=click.Path())
def parse_assembly_info(filename, output):
    """
    Parse the assembly summaries provided by NCBI to generate an sqlite database
    that can be used to quickly lookup assembly summary information. The
    assembly summary information is used in different parts of the downloading
    pipeline.

    \b
    Arguments
      FILENAME  Path to the TSV file with NCBI assembly summaries to parse
      OUTPUT    Path to write an sqlite database to
    """

    reader = csv.DictReader(filename, delimiter="\t")
    with SqliteDict(output) as db:
        count = 0
        for index, row in enumerate(reader):
            try:
                summary = ncbi.NcbiAssemblySummary.from_ncbi_row(row)
            except Exception:
                LOGGER.error("Could not parse row %s", row)
                continue
            if summary.assembly_accession in db:
                if db[summary.assembly_accession] != summary:
                    LOGGER.debug(
                        "Found duplicate entry for %s", summary.assembly_accession
                    )
                    continue
                raise ValueError("Multiple mappings to %s" % summary.assembly_accession)

            count += 1
            db[summary.assembly_accession] = summary
            if index % 100000 == 0:
                db.commit()
        LOGGER.info("Done parsing assemblies")
        db.commit()

        if count == 0:
            raise ValueError("Did not load any assemblies")

    with SqliteDict(output, flag="r") as db:
        assert len(db) == count, "Did not load all assemblies"


if __name__ == "__main__":
    cli()
