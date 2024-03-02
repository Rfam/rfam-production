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
import typing as ty
from pathlib import Path

import cattrs
import click
from boltons.fileutils import atomic_save
from boltons.jsonutils import JSONLIterator
from sqlitedict import SqliteDict

from rfamseq import download, ncbi
from rfamseq import uniprot as uni
from rfamseq.converter import camel_case_converter

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
    """Main entry point for updating rfamseq for Rfam. This is where the
    main logic for fetching sequences and building the metadata for the genome
    live.
    """
    if log_file:
        path = Path(log_file)
        if not path.parent.exists():
            path.parent.mkdir(parents=True)

    logging.basicConfig(
        filename=log_file,
        level=getattr(logging, log_level.upper()),
        # format="%(asctime)s %(levelname)-8s %(name)-15s %(message)s",
    )


@cli.group("uniprot")
def uniprot():
    """Commands dealing with fetching and processing data from UniProt. This is
    generally focused around getting proteomes, downloading their genomes and
    then building metadata for the Rfam database from this.
    """


@uniprot.command("fetch-proteomes")
@click.argument("requested", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def fetch_cmd(requested: ty.IO, output: ty.IO):
    """Fetch the proteomes specificed in the requested file. This file should
    be a list of UPID to fetch, one per line with no other content on each
    line. This is written to the given output file which is overwritten.

    \b
    Arguments:
      REQUESTED       The file of UPID's to fetch
      OUTPUT          The output file to write to, defaults stdout
    """
    ids = set(l.strip() for l in requested)
    for proteome in uni.proteome.fetch(ids):
        LOGGER.debug("Saving %s", proteome.id)
        data = cattrs.unstructure(proteome)
        output.write(json.dumps(data))
        output.write("\n")


@uniprot.command("deduplicate")
@click.argument("files", nargs=-1)
@click.argument("output", default="-", type=click.File("w"))
def dedup_cmd(files: ty.List[str], output: ty.IO):
    """Parse all given files to remove all duplicate proteomes. Proteomes are
    duplicated if they have the same id.
    """

    seen: ty.Set[str] = set()
    converter = camel_case_converter()
    for file in files:
        with open(file, "r") as raw:
            for raw in JSONLIterator(raw):
                proteome = converter.structure(raw, uni.proteome.Proteome)
                if proteome.id in seen:
                    continue
                seen.add(proteome.id)
                json.dump(raw, output)
                output.write("\n")


@uniprot.command("download-genomes")
@click.option(
    "--failed-filename",
    default="failed.jsonl",
    help="Name of the file to write failed entries to",
)
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
def download_cmd(
    version: str,
    ncbi_info: str,
    proteome_file: str,
    output: str,
    failed_filename="failed.jsonl",
):
    """Download the genomes specified in proteome-file. The file is the result
    of parse-proteomes, fetch-proteomes or a chunk of one of those files. This
    will download the entire genome for each proteome. This produces a UP*.fa
    and UP*.jsonl file for all proteomes in the file. The .jsonl contains
    metadata suitable for import to the Rfam database and .fa contains the
    sequences for the genome. The output directory will be created if needed.

    \b
    Arguments:
      VERSION        Version of the rfamseq database, eg 15.0
      NCBI-INFO      Path to file produced by `ncib parse-assembly-summary`
      PROTEOME-FILE  JSONL file of proteomes to download.
      OUTPUT         Path to write the sequence fasta and metadata jsonl files to.

    \b
    Options:
      --failed-filename  Name of the file to write the proteomes which were not
                         downloaded to. This will be a JSONL file with a 'status'
                         field showing the status and then the data itself.
    """
    out = Path(output)
    if not out.exists():
        out.mkdir(parents=True)

    converter = camel_case_converter()
    failed = out / failed_filename
    with SqliteDict(ncbi_info, flag="r") as db, failed.open("w") as fail:
        downloader = download.GenomeDownloader.build(version, db)
        for line in proteome_file:
            raw = json.loads(line)
            proteome = converter.structure(raw, uni.proteome.Proteome)
            LOGGER.info("Downloading %s", proteome.id)

            fasta_out = out / f"{proteome.id}.fa"
            metadata_out = out / f"{proteome.id}.jsonl"
            try:
                with atomic_save(str(fasta_out), text_mode=True) as fasta, atomic_save(
                    str(metadata_out), text_mode=True
                ) as meta:
                    downloader.download_proteome(proteome, fasta, meta)
            except download.SuppressedGenome:
                LOGGER.warning("Skipping proteome with suppressed genome %s", proteome)
                json.dump({"status": "suppressed", "data": raw}, fail)
                fail.write("\n")
                continue
            except (ncbi.ftp.UnknownGCA, ncbi.ftp.UnknownGCF):
                LOGGER.warning("Skipping proteome with unknown genome %s", proteome)
                json.dump({"status": "unknown-genome", "data": raw}, fail)
                fail.write("\n")
                continue
            except Exception as err:
                LOGGER.exception(err)
                LOGGER.error("Failed to download %s", proteome)
                json.dump({"status": "failed", "data": raw}, fail)
                fail.write("\n")


@uniprot.command("parse-proteomes")
@click.option(
    "--ignore",
    default=None,
    type=click.File("r"),
    help="A file with one UPI per line of UPIs to ignore when parsing.",
)
@click.argument(
    "jsonl",
    type=click.Path(),
)
@click.argument(
    "output",
    default="-",
    type=click.File("w"),
)
def p2g_cmd(jsonl, output, ignore=None):
    """Parse the JSON file provided by uniprot to a jsonl file. The jsonl file
    contains one JSON encoded object per line where the object contains
    everything needed to download a genome from NCBI/ENA.

    \b
    Arguments:
      JSONL   Path to the proteome.xml file fetched from uniprot.
      OUTPUT  Path to write a summary jsonl to

    \b
    Options:
      --ignore   Filename of UPIs to ignore, one per line
    """
    to_skip = set()
    if ignore:
        to_skip.update(l.strip() for l in ignore)

    converter = camel_case_converter()
    for proteome in uni.proteome.parse(Path(jsonl), to_skip):
        LOGGER.debug("Saving on %s", proteome.id)
        data = converter.unstructure(proteome)
        output.write(json.dumps(data))
        output.write("\n")


@cli.group("ncbi")
def ncbi_cmd():
    """This is a group of commands dealing with parsing data from NCBI."""


@ncbi_cmd.command("parse-assembly-summary")
@click.option("--fail-on-duplicate", is_flag=True)
@click.argument("filenames", nargs=-1)
@click.argument("output", type=click.Path())
def parse_assembly_info(filenames, output, fail_on_duplicate=False):
    """Parse the assembly summaries provided by NCBI to generate an sqlite
    database that can be used to quickly lookup assembly summary information.
    The assembly summary information is used in different parts of the
    downloading pipeline. The files are parsed in the order provided and if
    there are duplicate ids earlier ones take priority.

    \b
    Arguments
      FILENAMEs  Path to the TSV files with NCBI assembly summaries to parse
      OUTPUT     Path to write an sqlite database to

    \b
    Options
        --fail-on-duplicate    Fail if a duplicate assembly is found
    """

    count = 0
    with SqliteDict(output) as db:
        for index, summary in enumerate(ncbi.parse_assembly_files(filenames)):
            if summary.assembly_accession in db:
                LOGGER.warning(
                    "Found duplicate entry for %s", summary.assembly_accession
                )
                if db[summary.assembly_accession] != summary:
                    LOGGER.debug("Have: %s", db[summary.assembly_accession])
                    LOGGER.debug("Found: %s", summary.assembly_accession)
                if fail_on_duplicate:
                    raise ValueError(
                        "Duplicate accession for {summary.assembly_accession}"
                    )
                continue

            db[summary.assembly_accession] = summary
            count += 1
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
