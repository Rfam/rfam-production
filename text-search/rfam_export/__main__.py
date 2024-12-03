# -*- coding: utf-8 -*-

# Copyright [2009-2024] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import typing as ty
from pathlib import Path

import click
from loguru import logger
from rfam_export import families, full_region, genomes
from rfam_export.context import Context


@click.group()
@click.option("--rfam-version", envvar="RFAM_VERSION", help="Version of Rfam to use")
@click.option(
    "--rfamseq-version", envvar="RFAMSEQ_VERSION", help="Version of Rfamseq to use"
)
@click.option("--host", envvar="DB_HOST", help="Database host to use")
@click.option("--user", envvar="DB_USER", help="Database user to use")
@click.option("--password", envvar="DB_PASSWORD", help="Database password to use")
@click.option("--database-name", envvar="DB_DATABASE", help="Database name to use")
@click.option("--port", type=int, envvar="DB_PORT", help="Database port to use")
@click.option(
    "--log-level",
    default="INFO",
    help="Log level to use",
    type=click.Choice(
        ["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"],
        case_sensitive=False,
    ),
)
@click.pass_context
def main(
    ctx,
    host: str,
    user: str,
    password: str,
    database_name: str,
    port: int,
    rfam_version: str,
    rfamseq_version: str,
    log_level="INFO",
):
    """This is the main entry point to the text search export for Rfam.

    This requires defining what the Rfam and Rfamseq versions being used are.
    This can be done either by setting `--rfam-version` or the `$RFAM_VERSION`
    env var. Similiary, `--rfamseq-version` or `$RFAMSEQ_VERSION` must be set.

    Additionally, there must be DB_HOST,USER, etc variables to define what
    MySQL database to use to fetch data from.
    """
    logger.remove()
    logger.add(sys.stderr, level=log_level)

    if not rfam_version:
        logger.critical("Must specify the Rfam release version")
        ctx.exit(2)

    if not rfamseq_version:
        logger.critical("Must specify the rfamseq version")
        ctx.exit(2)

    ctx.ensure_object(dict)
    ctx.obj["context"] = Context(
        rfam_version=rfam_version,
        rfamseq_version=rfamseq_version,
        mysql_host=host,
        mysql_user=user,
        mysql_password=password,
        mysql_port=port,
        mysql_database=database_name,
    )


@main.command("list-genomes")
@click.argument("output", default="-", type=click.File("w"))
@click.pass_context
def list_genomes_cmd(ctx, output: ty.TextIO):
    """Create a listing of all genomes to export for the given Rfam and Rfamseq
    versions.

    This produces a JSONL file which can be used by `full-regions`.

    Arguments

    ---------

      output - The file to output to, default stdout.
    """
    context = ctx.obj["context"]
    with context.cursor() as cursor:
        for accession in genomes.all_genome_accessions(cursor, context):
            output.write(accession)
            output.write("\n")


@main.command("full-regions")
@click.option("--failed-file", default="failed.txt", type=click.File("w"))
@click.argument("genome-file", type=click.File("r"))
@click.argument(
    "output",
    default=".",
    type=click.Path(dir_okay=True, file_okay=False, path_type=Path),
)
@click.pass_context
def full_regions_cmd(ctx, genome_file: ty.TextIO, output: Path, failed_file: ty.TextIO):
    """Get all sequence matches for all genomes in the given file.

    This will produce a XML file suitable for import by the EBeye search
    service. The input file should be the reuslt of `list-genomes` above. The
    output files will be named like: `accession.xml` and placed in the given
    directory.

    Arguments

    ---------

      genome-file - The list of genomes to import.

      output - The output directory to write to, defaults current directory.
    """
    context: Context = ctx.obj["context"]
    logger.info("Starting to export all full regions in file {}", genome_file.name)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    with context.cursor() as cursor:
        context.cache("families", families.family_mapping(cursor))
        count = 0
        attempted = 0
        for line in genome_file:
            accession = line.strip()
            logger.info("Starting to export all full regions from {}", accession)
            attempted += 1
            out = output / f"{accession}.xml"
            try:
                genome = genomes.genome_info(cursor, context, accession)
                logger.info("Writing genome {} to {}", genome.accession, out)
                with out.open("w") as handle:
                    full_region.export_region(cursor, context, genome, handle)
                    count += 1
            except Exception as err:
                out.unlink(missing_ok=True)
                match err:
                    case full_region.NoHits():
                        logger.info("Genome {} has no full_regions", accession)
                        failed_file.write(f"{accession} no_hits\n")
                    case _:
                        logger.exception(err)
                        logger.info("Failed to export {}", accession)
                        failed_file.write(f"{accession} other\n")
                failed_file.flush()
    logger.info("Exported {} of {} genomes", count, attempted)
    if count == 0:
        logger.error("Failing, since no genomes exported")
        sys.exit(70)


@main.command("export-families")
@click.argument(
    "output",
    default=".",
    type=click.Path(dir_okay=True, file_okay=False, path_type=Path),
)
@click.pass_context
def export_families_cmd(ctx, output: Path):
    """Export all families into the given directory."""
    context: Context = ctx.obj["context"]
    logger.info("Starting to export all families to {}", output)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)

    attempted = 0
    count = 0
    with context.cursor() as cursor:
        all_families = families.families(cursor)
        exporter = families.Exporter.build(context, cursor, output)
        for family in all_families:
            attempted += 1
            try:
                exporter.export_family(family)
            except Exception as err:
                logger.info("Failed to export {}", family)
                logger.exception(err)
            count += 1

    logger.info("Exported {} of {} families", count, attempted)
    if count == 0:
        logger.error("Failing, since no families exported")
        sys.exit(70)


if __name__ == "__main__":
    main()
