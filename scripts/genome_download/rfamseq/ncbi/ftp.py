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
from contextlib import closing, contextmanager
from ftplib import FTP
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


class UnknownGCA(Exception):
    """
    Raised if an Unknown GCA id is given
    """


class InvalidGenomeId(Exception):
    """
    Raised if given a non GCA/GCF id.
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


def increment_version(id: str) -> str:
    matches = re.match(r"^(.+)\.(\d+)$", id)
    if not matches:
        raise ValueError(f"Failed to parse {id} to increment")
    prefix = matches.group(1)
    version = matches.group(2)
    next_version = int(version) + 1
    return f"{prefix}.{next_version}"


def ftp_path(
    info: SqliteDict, accession: str, suffix: str, try_next_version=True
) -> ty.Optional[str]:
    versioned = add_version_if_missing(info, accession)
    if not versioned or versioned not in info:
        LOGGER.warning("Accession %s not found in ncbi db", versioned)
        return None

    path = info[versioned].ftp_path
    if path == "na" or path is None:
        LOGGER.warning("Accession %s has no path", versioned)
        LOGGER.debug("Have %s", info[versioned])
        if not try_next_version:
            return None
        next_version = increment_version(versioned)
        return ftp_path(info, next_version, suffix, try_next_version=False)
    parts = path.split("/")
    name = parts[-1]
    return f"{path}/{name}_{suffix}"


def guess_genome_url(accession: str, suffix: str) -> ty.Optional[str]:
    LOGGER.info("Trying to fetch NCBI FTP url for %s", accession)
    if accession[0:4] not in {"GCA_", "GCF_"}:
        LOGGER.debug("Accession %s is not a genome accession", accession)
        return None

    if len(accession) not in {13, 15}:
        LOGGER.debug("Accession %s is not expected size", accession)
        return None

    host = "ftp.ncbi.nlm.nih.gov"
    try:
        with closing(FTP(host)) as ftp:
            ftp.login()
            ftp_dir = f"genomes/all/{accession[0:3]}/{accession[4:7]}/{accession[7:10]}/{accession[10:13]}"
            ftp.cwd(ftp_dir)
            possible = ftp.nlst()
            matching = [p for p in possible if p.startswith(accession)]
            if len(matching) == 1:
                return f"ftp://{host}/{ftp_dir}/{matching[0]}/{matching[0]}{suffix}"
            else:
                LOGGER.info("Cannot determine which of %s matches", possible)
    except Exception as err:
        LOGGER.debug("Failed to access NCBI FTP")
        LOGGER.debug(err)
        return None
    return None


def ftp_fasta_url(info: SqliteDict, accession: str) -> str:
    prefix = accession[0:4]
    if prefix not in {"GCA_", "GCF_"}:
        raise InvalidGenomeId(accession)

    if (url := ftp_path(info, accession, "genomic.fna.gz")) is not None:
        return url

    if (url := guess_genome_url(accession, "_genomic.fna.gz")) is not None:
        return url

    if prefix == "GCF_":
        raise UnknownGCF(accession)
    raise UnknownGCA(accession)


def efetch_fasta_url(accession: str) -> str:
    return NCBI_SEQ_URL.format(accession=accession)


def fasta_url(info: SqliteDict, accession: str) -> str:
    if accession.startswith("GCA_") or accession.startswith("GCF_"):
        LOGGER.info("Trying FTP access to %s", accession)
        return ftp_fasta_url(info, accession)
    return efetch_fasta_url(accession)


@contextmanager
def fetch_fasta(info: SqliteDict, accession: str) -> ty.Iterator[ty.IO]:
    with wget.wget(fasta_url(info, accession)) as handle:
        yield handle


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
