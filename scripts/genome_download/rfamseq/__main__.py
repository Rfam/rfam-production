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
import math
import typing as ty
from pathlib import Path
from string import Template

import cattrs
import click
from boltons.fileutils import atomic_save
from boltons.jsonutils import JSONLIterator
from humanfriendly import parse_size
from sqlitedict import SqliteDict

from rfamseq import download, mgnify, ncbi, seqstat
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
    converter = camel_case_converter()
    for proteome in uni.proteome.fetch(ids):
        LOGGER.debug("Saving %s", proteome.id)
        data = converter.unstructure(proteome)
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
    count = 0
    converter = camel_case_converter()
    for file in files:
        LOGGER.info("Deduplicating file %s", file)
        with open(file, "r") as raw:
            for raw in JSONLIterator(raw):
                try:
                    proteome = converter.structure(raw, uni.proteome.Proteome)
                except Exception as err:
                    LOGGER.error("Failed to handle %s", raw)
                    raise err

                LOGGER.debug("Validating %s", proteome.id)
                if proteome.id in seen:
                    LOGGER.debug("Skipping duplicate genome %s", proteome.id)
                    count += 1
                    continue
                seen.add(proteome.id)
                json.dump(raw, output)
                output.write("\n")
    LOGGER.info("Saw %i duplicates", count)


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
                    if not fasta:
                        raise ValueError(f"Failed to open fasta file {fasta_out}")
                    if not meta:
                        raise ValueError(f"Failed to open metadata file {metadata_out}")
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
    "reference",
    type=click.Path(),
)
@click.argument(
    "output",
    default="-",
    type=click.File("w"),
)
def p2g_cmd(reference, output, ignore=None):
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
    for proteome in uni.proteome.parse(Path(reference), to_skip):
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


@cli.group("mgnify")
def mgnify_cmd():
    """This is a set of commands dealing with selecting and downloading mgnify
    genomes."""


@mgnify_cmd.command("select-mags")
@click.argument("summary-file", type=click.File("r"))
@click.argument("metric-file", type=click.File("r"))
@click.argument("output-file", default="-", type=click.File("w"))
def select_mags(summary_file, metric_file, output_file):
    """This is a command to select mags matching some quality cutoffs and
    completeness. This will load all MAG information from the given summary-file
    and then filter it to only those which match the criteria in the
    metric-file. This writes a JSONL formatted file out.

    ARGUMENTS:
      summary-file   A TSV file that summarises the status of all MAGs. This is
                     provided by MGnify.
      metric-file    A JSON file which lists the quality cutoffs to use.
      output         The file to write to, defaults to stdout.
    """

    metric = cattrs.structure(json.load(metric_file), mgnify.Selector)
    reader = csv.DictReader(summary_file, delimiter="\t")
    count = 0
    accepted = 0
    for row in reader:
        count += 1
        mag = mgnify.MAGInfo.from_mgnify(row)
        if metric.accepts(mag):
            accepted += 1
            json.dump(cattrs.unstructure(mag), output_file)
            output_file.write("\n")
    LOGGER.info("Accepted %i of %i sequences", accepted, count)


@mgnify_cmd.command("download")
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
    "mag-file",
    default="-",
    type=click.File("r"),
)
@click.argument(
    "output",
    default=".",
    type=click.Path(),
)
def mgnify_download(
    version, ncbi_info, mag_file, output, failed_filename="failed.jsonl"
):
    out = Path(output)
    if not out.exists():
        out.mkdir(parents=True)

    failed = out / failed_filename
    with SqliteDict(ncbi_info, flag="r") as db, failed.open("w") as fail:
        downloader = download.GenomeDownloader.build(version, db)
        for line in mag_file:
            raw = json.loads(line)
            mag = cattrs.structure(raw, mgnify.MAGInfo)
            LOGGER.info("Downloading %s", mag.accession)

            fasta_out = out / f"{mag.accession}.fa"
            metadata_out = out / f"{mag.accession}.jsonl"
            try:
                with atomic_save(str(fasta_out), text_mode=True) as fasta, atomic_save(
                    str(metadata_out), text_mode=True
                ) as meta:
                    if not fasta:
                        raise ValueError(f"Failed to open fasta file {fasta_out}")
                    if not meta:
                        raise ValueError(f"Failed to open metadata file {metadata_out}")
                    downloader.download_mag(mag, fasta, meta)
            except Exception as err:
                LOGGER.exception(err)
                LOGGER.error("Failed to download %s", mag)
                json.dump({"status": "failed", "data": raw}, fail)
                fail.write("\n")


@cli.group("seqstat")
def seqstat_group():
    """A series of commands dealing with seqstat files. These are produced by
    esl-seqstat and summarize the sequences in a fasta file.
    """
    pass


@seqstat_group.command("chunk-ids")
@click.option(
    "--prefix",
    default="chunk_",
    help="The prefix for each chunk filename to use",
)
@click.option(
    "--extension",
    default=".txt",
    help="The file extension for each chunk to use",
)
@click.argument("seqstat-file", type=click.File("r"))
@click.argument("chunk-size")
@click.argument("output", default=".", type=click.Path())
def chunk_ids_cmd(seqstat_file, chunk_size, output, prefix="chunk_", extension=".fa"):
    """Chunk the given seqstat file into several files where each file only
    contains the ids from the file and the size of the file that will be
    produced by fetching those sequences is roughly chunk-size. Each chunk will
    be named: ${prefix}${index}${extension}. The indexes are 1 indexed. This
    assumes that each sequence is about 1 bytes per nucleotide and uses that to
    compute the expected size of the resulting file. The output directory will
    be created if needed.

    Arguments:
      seqstat-file    The file produced by esl-seqstat.
      chunk-size      The size of the chunk, may be values like '1GB', '10MB'.
      output          The output directory to write chunks to
    """

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)

    byte_size = parse_size(chunk_size)
    chunks = seqstat.chunk(seqstat.parse(seqstat_file), byte_size)
    for index, chunk in enumerate(chunks):
        out_path = output / f"{prefix}{index + 1}{extension}"
        with open(out_path, "w") as out:
            for info in chunk.entries:
                out.write(info.sequence_id)
                out.write("\n")


@seqstat_group.command("take-fraction")
@click.argument("seqstat-file", type=click.Path())
@click.argument("fraction", type=float)
@click.argument("output", default="-", type=click.File("w"))
def take_fraction_cmd(seqstat_file, fraction, output):
    """This takes the given fraction of entries from the seqstat file and writes
    them out. The entries are taken from the start of the file and

    Arguments:
      seqstat-file    The file produced by esl-seqstat.
      fraction        The fraction of the file to take, between 0 and 1
      output          The output file to write to, default stdout.
    """
    assert (
        0 < fraction < 1
    ), f"Fraction ${fraction} must be between 0 and 1, exclusively"

    count = 0
    with open(seqstat_file, "r") as raw:
        for line in raw:
            if line.startswith("="):
                count += 1

    if not count:
        raise ValueError("Cannot work with an empty file")
    LOGGER.info("Found %i seqstat entries", count)

    take = round(count * fraction)
    if not take:
        raise ValueError(f"Asked to take ${count} * ${fraction} = ${take}")
    LOGGER.info("Will take %i lines", take)

    with open(seqstat_file, "r") as raw:
        written = 0
        for index, line in enumerate(raw):
            if index >= take:
                break
            output.write(line)
            written += 1

        if written != take:
            raise ValueError(f"Only wrote {written} of {take}")


@cli.command("database-size")
@click.argument("seqstat-files", nargs=-1, type=click.File("r"))
def db_size_cmd(seqstat_files):
    """Compute the database size based upon the given seqstat files."""
    total = 0
    for file in seqstat_files:
        total += 2 * seqstat.total_residues(file)
    size = total / 1_000_000.0
    print(math.ceil(size))


@cli.command("build-config")
@click.option("--define", multiple=True, help="Define a key=value pair")
@click.argument("version")
@click.argument("template-file", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def build_command(version, template_file, output, define=None):
    """Build a new Rfamseq database. This will fill out the template file to
    create the required structure. If there are keys from the template not
    present in the --define options then this will fail.

    Arguments:
      version        Rfam version
      template-file  Path to the template file to fill out

    Options:
      --define key=value   Define some key value pairs for the template. This
                           can be done as many times as needed.
    """

    data = {}
    if define:
        for entry in define:
            key, value = entry.split("=", 1)
            data[key] = value

    data["version"] = version

    template = Template(template_file.read())
    assert template.is_valid(), "Invalid template"
    output.write(template.substitute(data))


if __name__ == "__main__":
    cli()
