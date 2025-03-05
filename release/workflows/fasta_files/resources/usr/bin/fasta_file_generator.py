#!/usr/bin/env python3

# Copyright [2009-2025] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#      http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Script to generate fasta files for all family regions in full_region.

Easel tools should be installed and added to PATH. Easel tools can be installed
along with the Infernal suite
"""

import os
import re
import subprocess
import typing as ty
from pathlib import Path

import click
import pymysql.cursors

BAD_SEQUENCE = re.compile(r"[.-@|\s| -)|z-~|Z-`|EFIJLOPQX|efijlopqx+,]+")


def sfetch(database, accession, substring=None) -> str:
    cmd = ["esl-sfetch", database, accession]
    if substring:
        cmd = ["esl-sfetch", "-c", substring, database, accession]
    output = subprocess.check_output(cmd, text=True)
    parts = output.split()
    return "\n".join(parts[1:])


def extract_full_sequences(
    cursor, seq_db: Path, rfam_acc: str
) -> ty.Iterator[ty.Tuple[str, str, str]]:
    query = """
        SELECT
            fr.rfam_acc,
            fr.rfamseq_acc,
            fr.seq_start,
            fr.seq_end,
            rf.description
        FROM full_region fr JOIN rfamseq rf
        ON
            fr.rfamseq_acc=rf.rfamseq_acc
        WHERE
            fr.is_significant = 1
            AND fr.type = 'full'
            AND fr.rfam_acc = %s
    """

    cursor.execute(query, (rfam_acc,))
    for region in cursor:
        accession = region["rfamseq_acc"]
        substring = "{start}-{stop}".format(
            start=region["seq_start"],
            stop=region["seq_end"],
        )

        yield (
            f"{accession}/{substring}",
            region["description"],
            sfetch(seq_db, accession, substring=substring),
        )


def extract_seed_sequences(
    cursor, rfam_seed: Path, rfam_acc: str
) -> ty.Iterator[ty.Tuple[str, str, str]]:
    query = """
        SELECT
            sr.rfamseq_acc,
            sr.seq_start,
            sr.seq_end,
            rf.description
        FROM seed_region sr JOIN rfamseq rf
        ON
            sr.rfamseq_acc = rf.rfamseq_acc
        WHERE
            sr.rfam_acc = %s
    """

    cursor.execute(query, (rfam_acc,))
    for region in cursor:
        accession = "{accession}/{start}-{stop}".format(
            accession=region["rfamseq_acc"],
            start=region["seq_start"],
            stop=region["seq_end"],
        )
        yield accession, region["description"], sfetch(rfam_seed, accession)


def extract_family_sequences(
    connection, seq_db: Path, rfam_seed: Path, rfam_acc: str, output: ty.TextIO
):
    with connection.cursor() as cursor:
        for id, name, sequence in extract_seed_sequences(cursor, rfam_seed, rfam_acc):
            if not sequence or BAD_SEQUENCE.search(sequence) is not None:
                raise ValueError(f"Invalid SEED sequence for {id}")
            output.write(">%s %s\n%s\n" % (id, name, sequence))

    with connection.cursor() as cursor:
        for id, name, sequence in extract_full_sequences(cursor, seq_db, rfam_acc):
            if not sequence or BAD_SEQUENCE.search(sequence) is not None:
                raise ValueError(f"Invalid FULL sequence for {id}")
            output.write(">%s %s\n%s\n" % (id, name, sequence))


@click.command("")
@click.option("--host")
@click.option("--port", type=int)
@click.option("--user")
@click.option("--password")
@click.option("--database")
@click.argument("sequence-database", type=click.Path())
@click.argument("rfam-seed", type=click.Path())
@click.argument("rfam-accession")
@click.argument("output", type=click.File("w"), default="-")
def main(
    sequence_database,
    rfam_seed,
    rfam_accession,
    output,
    host=None,
    port=None,
    user=None,
    password=None,
    database=None,
):
    """Generate a fasta file for the given accessions and write it to output.
    If any sequence cannot be exported or is invalid this will crash. This may
    produce incomplete files in the case of a crash.

    Arguments\b
    ----------\b
    sequence-database: The path to the combined Rfam.fa in Rfamseq.\b
    rfam-seed: The path to Rfam.seed for this release.\b
    rfam-accession: The Rfam accession to write sequences for.\b
    output: The output file to write to.\b
    """
    connection = pymysql.connect(
        host=host,
        port=port,
        user=user,
        password=password,
        database=database,
        cursorclass=pymysql.cursors.DictCursor,
    )
    with connection:
        extract_family_sequences(
            connection, Path(sequence_database), Path(rfam_seed), rfam_accession, output
        )
    if not os.path.exists(output) or not os.path.getsize(output) > 0:
        raise ValueError(f"File {output} is empty")


if __name__ == "__main__":
    main()
