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

import logging
import typing as ty
from functools import lru_cache

from Bio import SeqIO

from rfamseq import fasta, wget

LOGGER = logging.getLogger(__name__)

ENA_FASTA_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=true"
ENA_EMBL_URL = "https://www.ebi.ac.uk/ena/browser/api/embl/{accession}?download=true"
ENA_WGS_FASTA_URL = (
    "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/{prefix}/{name}.fasta.gz"
)


class MissingContigs(Exception):
    """
    Raised if there are no contigs in the EMBL file.
    """


def fetch(template: str, **data) -> ty.Iterable[SeqIO.SeqRecord]:
    url = template.format(**data)
    LOGGER.debug("Fetching %s", url)
    with wget.wget(url) as handle:
        yield from fasta.parse(handle)


@lru_cache
def fetch_contigs(accession: str) -> ty.List[str]:
    LOGGER.info("Fetching EMBL formatted file for %s", accession)
    url = ENA_EMBL_URL.format(accession=accession)
    contigs = []
    with wget.wget(url) as handle:
        for line in handle:
            if not line.startswith("CON"):
                continue
            contigs.append(line[3:].strip())
    if not contigs:
        raise MissingContigs(f"No contigs for {accession}")
    return contigs


def fetch_fasta(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    LOGGER.info("Fetching %s fasta from ENA", accession)
    yield from fetch(ENA_FASTA_URL, accession=accession)


def wgs_fasta_url(accession: str) -> str:
    prefix = accession[0:3].lower()
    name = accession[0:6].upper()
    return ENA_WGS_FASTA_URL.format(prefix=prefix, name=name)


def fetch_wgs_sequences(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    LOGGER.info("Fetching the wgs fasta set for %s", accession)
    yield from fetch(wgs_fasta_url(accession))


def lookup(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    LOGGER.info("Fetching %s from ENA", accession)
    try:
        yield from fetch_fasta(accession)
    except wget.FetchError:
        LOGGER.info("Failed to get directly, will try via contigs")
        contigs = list(fetch_contigs(accession))
        if not contigs:
            raise ValueError(f"Could not find contigs for {accession}")
        for contig in contigs:
            LOGGER.info("Fetching contig %s for %s", contig, accession)
            yield from fetch_fasta(contig)