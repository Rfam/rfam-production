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

import io
import logging
import subprocess as sp
import typing as ty
from pathlib import Path

from Bio import SeqIO

LOGGER = logging.getLogger(__name__)


def extract_ids(path: Path) -> ty.Set[str]:
    """
    Extract all sequence ids from the fasta file at the given path. This uses
    seqkit to parse the file and extract the ids quickly. This will fail of no
    ids are extracted.
    """

    ids = set()
    text = sp.check_output(
        ["seqkit", "seq", "--name", "--only-id", str(path)], text=True
    )
    for line in io.StringIO(text):
        ids.add(line.strip())
    if not ids:
        raise ValueError(f"Somehow read no ids from {path}")
    return ids


def parse(handle: ty.IO) -> ty.Iterable[SeqIO.SeqRecord]:
    """
    Parses a fasta file and removes any entries with duplicate ids, handle ids
    that use the '|' style and use the third section as the id and ensure that
    the handle contains fasta entries.
    """

    seen_ids = set()
    for record in SeqIO.parse(handle, "fasta"):
        if "|" in record.id:
            record.id = record.id.split("|")[2]
        if record.id in seen_ids:
            LOGGER.error(f"Fasta file contains duplicate id %s", record.id)
            continue
        if not record.seq:
            LOGGER.error(f"Fasta file contains empty sequence %s", record.id)
            continue
        yield record
        seen_ids.add(record.id)
    if not seen_ids:
        raise ValueError("Did not parse any sequences")
