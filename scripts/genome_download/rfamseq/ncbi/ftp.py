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

from __future__ import annotations

import csv
import logging
import re
import typing as ty
from io import StringIO

import requests
from Bio import SeqIO
from ratelimit import limits, sleep_and_retry
from sqlitedict import SqliteDict

from rfamseq import fasta, wget

LOGGER = logging.getLogger(__name__)

NCBI_SEQ_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"

NCBI_SUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={accessions}"

NCBI_WGS_URL = "https://www.ncbi.nlm.nih.gov/Traces/wgs/{accession}/contigs/tsv"


class UnknownGCF(Exception):
    """
    Raised if an Unknown GCF id is given
    """


class InvalidGenomeId(Exception):
    """
    Raised if given a non GCA/GCF id.
    """


class UnknownGenomeId(Exception):
    """
    Raised if the genome id looks valid but NCBI assembly data does not know
    about it.
    """


def add_version_if_missing(info: SqliteDict, id: str) -> str:
    if "." in id:
        return id
    possible = {}
    pattern = re.compile(f"^{id}.(\\d+)$")
    for key in info.iterkeys():
        if match := re.match(pattern, key):
            index = int(match.group(1))
            possible[index] = key
    to_use = max(possible.keys())
    return possible[to_use]


def ftp_path(info: SqliteDict, accession: str, suffix: str) -> ty.Optional[str]:
    versioned = add_version_if_missing(info, accession)
    if not versioned or versioned not in info:
        LOGGER.info("Accession %s not found in ncbi db", versioned)
        return None

    path = info[versioned].ftp_path
    if path == "na" or path is None:
        LOGGER.info("Accession %s has no path", versioned)
        return None
    parts = path.split("/")
    name = parts[-1]
    return f"{path}/{name}_{suffix}"


def genome_ftp_path(info: SqliteDict, accession: str) -> ty.Optional[str]:
    return ftp_path(info, accession, "genomic.fna.gz")


@sleep_and_retry
@limits(3, period=1)
def efetch_fasta(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    LOGGER.info("Trying efetch for %s", accession)
    url = NCBI_SEQ_URL.format(accession=accession)
    try:
        with wget.wget(url) as handle:
            yield from fasta.parse(handle)
    except wget.FetchError as err:
        LOGGER.debug(err)


def ftp_fasta(info: SqliteDict, accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    LOGGER.info("Trying FTP access to %s", accession)
    prefix = accession[0:4]
    if prefix not in {"GCA_", "GCF_"}:
        raise InvalidGenomeId(accession)
    url = genome_ftp_path(info, accession)
    if url is None and prefix == "GCF_":
        raise UnknownGCF(accession)
    if url is None:
        raise Exception("Not yet implemented")
    with wget.wget(url) as handle:
        yield from fasta.parse(handle)


def fetch_fasta(info: SqliteDict, accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    if accession.startswith("GCA_") or accession.startswith("GCF_"):
        yield from ftp_fasta(info, accession)
    else:
        yield from efetch_fasta(accession)


def resolve_wgs(accession: str) -> ty.Optional[ty.List[str]]:
    LOGGER.info("Trying to resolve WGS set %s", accession)
    url = NCBI_WGS_URL.format(accession=accession)
    LOGGER.debug("Fetching %s", url)
    response = requests.get(url)
    try:
        response.raise_for_status()
    except:
        LOGGER.info("Request to resolve wgs set failed")
        return None

    reader = csv.DictReader(StringIO(response.text), delimiter="\t")
    accessions = [r["accession"] for r in reader]
    if not accessions:
        LOGGER.info("Failed to get load any accessions from response")
        return None
    return accessions
