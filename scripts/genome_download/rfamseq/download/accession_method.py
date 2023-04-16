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
from contextlib import contextmanager

from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import ena, fasta, ncbi, wget

Records = ty.Iterable[SeqIO.SeqRecord]

LOGGER = logging.getLogger(__name__)


class FastaLookupFailed(Exception):
    """
    Raised when if all fasta lookup methods failed to fetch some data.
    """


class AccessionLookupFailed(Exception):
    """
    Raised when if all accession lookup methods failed to fetch some id.
    """


@contextmanager
def fetch(info: SqliteDict, accession: str) -> ty.Iterator[ty.IO]:
    try:
        LOGGER.info("Trying to fetch %s from NCBI", accession)
        with ncbi.fetch_fasta(info, accession) as handle:
            yield handle
            return
    except Exception as err:
        LOGGER.debug(err)
        LOGGER.debug("Failed to fetch from NCBI")

    try:
        LOGGER.info("Trying to fetch %s from ENA", accession)
        with ena.fetch_fasta(accession) as handle:
            yield handle
            return
    except Exception as err:
        LOGGER.debug(err)
        LOGGER.debug("Failed to fetch from ENA")

    raise AccessionLookupFailed(f"Failed to lookup {accession}")


def records(info: SqliteDict, accession: str) -> ty.Iterator[SeqIO.SeqRecord]:
    with fetch(info, accession) as handle:
        yield from fasta.parse(handle)
