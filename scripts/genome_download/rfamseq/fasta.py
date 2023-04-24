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

import typing as ty

from Bio import SeqIO


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
            raise ValueError(f"Duplicate id {record.id}")
        if not record.seq:
            raise ValueError(f"Empty sequence is invalid {record.id}")
        yield record
        seen_ids.add(record.id)
    if not seen_ids:
        raise ValueError("Did not parse any sequences")
