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

import logging
import typing as ty

from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import ncbi, uniprot

from .accession_method import fetch as fetch_accessions

Records = ty.Iterable[SeqIO.SeqRecord]

LOGGER = logging.getLogger(__name__)


def fetch(info: SqliteDict, proteome: uniprot.ProteomeInfo) -> Records:
    LOGGER.info("Looking up each component for %s", proteome.upi)
    genome = proteome.genome_info
    assert isinstance(
        genome.components, uniprot.SelectedComponents
    ), f"Invalid components for {proteome}"
    ids = []
    for component in genome.components:
        assert isinstance(component, str), f"Invalid component in {proteome}"
        ids.append(component)

    ids = ",".join(ids)
    LOGGER.info("Trying to lookup all ids as a batch: %s", ids)
    try:
        fetched = ncbi.efetch_fasta(ids)
        if fetched:
            yield from fetched
        else:
            raise ValueError("Failed to fetch using efetch")
    except:
        LOGGER.info("Failed to efetch all ids, will try individual lookup")
        for component in genome.components:
            assert isinstance(component, str), f"Invalid component in {proteome}"
            yield from fetch_accessions(info, component)
