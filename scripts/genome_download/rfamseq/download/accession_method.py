# -*- coding: utf-8 -*-

"""
Copyright [2009-${2023}] EMBL-European Bioinformatics Institute
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

from __future__ import annotations

import logging
import typing as ty

from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import ena, ncbi, wget

Records = ty.Iterable[SeqIO.SeqRecord]

LOGGER = logging.getLogger(__name__)


def lookup_fasta(info: SqliteDict, accession: str) -> Records:
    try:
        LOGGER.info("Trying to fetch %s from NCBI", accession)
        yield from ncbi.fetch_fasta(info, accession)
    except wget.FetchError:
        LOGGER.info("Failed to fetch %s from NCBI, trying ENA", accession)
        yield from ena.fetch_fasta(accession)


def fetch(info: SqliteDict, accession: str) -> Records:
    LOGGER.info("Trying to fetch %s from NCBI", accession)
    try:
        yield from lookup_fasta(info, accession)
    except wget.FetchError:
        LOGGER.info("Trying to find contigs from %s", accession)
        contigs = list(ena.fetch_contigs(accession))
        for contig in contigs:
            LOGGER.info("Fetching contig %s from NCBI", contig)
            try:
                yield from ncbi.efetch_fasta(contig)
            except wget.FetchError:
                LOGGER.info("Failed to get %s from NCBI", contig)
                try:
                    yield from ena.fetch_fasta(contig)
                except wget.FetchError:
                    LOGGER.info("Failed to get %s from ENA", contig)
