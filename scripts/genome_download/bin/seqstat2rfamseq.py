#!/usr/bin/env python3

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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
import typing as ty

import click
from sqlitedict import SqliteDict


def extract_metadata_from_fasta(
    seqstat_file: ty.IO, info: SqliteDict
) -> ty.Iterable[ty.List[str]]:
    """
    Parses a esl-seqstat file and generates rfamseq like matadata
    """

    previous_acc = ""
    mol_type = "genomic DNA"

    for line in seqstat_file:
        seq_elements = [x.strip() for x in line.strip().split(" ") if x != ""]

        rfamseq_acc = seq_elements[0]
        accession = rfamseq_acc.partition(".")[0]
        version = rfamseq_acc.partition(".")[2]
        length = seq_elements[1]
        description = " ".join(seq_elements[2:])

        yield [
            rfamseq_acc,
            accession,
            "00000" + str(version),
            info[rfamseq_acc]["taxid"],
            mol_type,
            length,
            description,
            previous_acc,
            info[rfamseq_acc]["source"],
        ]


@click.command()
@click.argument("info_file", type=click.Path())
@click.argument("seqstat-file", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def main(info_file: str, seqstat_file: ty.IO, output: ty.IO):
    writer = csv.writer(output, delimiter="\t")
    with SqliteDict(info_file) as db:
        metadata = extract_metadata_from_fasta(seqstat_file, db)
        writer.writerows(metadata)


if __name__ == "__main__":
    main()
