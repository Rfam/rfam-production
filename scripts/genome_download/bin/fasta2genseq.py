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


def build_genseq(
    version: str, db: SqliteDict, handle: ty.IO
) -> ty.Iterable[ty.List[str]]:
    for line in handle:
        if not line.startswith(">"):
            continue
        parts = line.split(" ")
        accession = parts[0][1:]
        info = db[accession]
        yield [
            info["upid"],
            accession,
            "",
            "",
            version,
        ]


@click.command()
@click.argument("version")
@click.argument("info-file", type=click.Path())
@click.argument("header-file", default="-", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def main(version, info_file, header_file, output):
    writer = csv.writer(output, delimiter="\t")
    with SqliteDict(info_file) as db:
        writer.writerows(build_genseq(version, db, header_file))


if __name__ == "__main__":
    main()
