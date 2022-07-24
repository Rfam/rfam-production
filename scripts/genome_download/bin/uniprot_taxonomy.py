#!/usr/bin/env python3

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

import click

NS = {
    "uni": "http://uniprot.org/uniprot",
    "pro": "http://uniprot.org/proteomes",
}

FIELDS = [
    "id",
    "common_name",
    "science_name",
    "lineage",
]


@click.command()
@click.argument("taxonomy", type=click.File("r"))
@click.argument("summary", type=click.File("r"))
@click.argument("output", type=click.File("w"))
def main(taxonomy, summary, output):
    logging.basicConfig(level=logging.INFO)
    taxids = set()
    for line in summary:
        info = json.loads(line)
        taxids.add(info["taxid"])

    reader = csv.DictReader(taxonomy, fieldnames=FIELDS, delim="\t")
    writer = csv.writer(output)
    for row in reader:
        if row["id"] not in taxids:
            continue
        tree_name = row["Scientific name"].replace(" ", "_")
        lineage = row["lineage"]
        lineage = "; ".join(reversed(lineage.split(", ")))
        writer.writerow(
            [
                row["id"],
                row["science_name"],
                lineage,
                tree_name,
                f"{tree_name}[{row['id']}]",
            ]
        )


if __name__ == "__main__":
    main()
