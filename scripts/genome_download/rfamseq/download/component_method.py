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

from rfamseq import fasta, ncbi, uniprot, wget, wgs
from rfamseq.accession import Accession

from .accession_method import records as accession_records

LOGGER = logging.getLogger(__name__)


def components(proteome: uniprot.ProteomeInfo) -> ty.List[str]:
    LOGGER.info("Looking up each component for %s", proteome.upi)
    genome = proteome.genome_info
    if not isinstance(genome.components, uniprot.SelectedComponents):
        raise ValueError(f"Invalid components for {proteome}")

    ids = []
    for component in genome.components:
        if isinstance(component, Accession):
            ids.append(str(component))
        elif isinstance(component, (wgs.WgsPrefix, wgs.WgsSequenceId)):
            ids.append(component.to_wgs_string())
        else:
            raise ValueError(f"Cannot fetch by component for {genome}")

    return ids


def records(
    info: SqliteDict, proteome: uniprot.ProteomeInfo
) -> ty.Iterator[SeqIO.SeqRecord]:
    comp_ids = components(proteome)
    ids = ",".join(comp_ids)
    LOGGER.info("Trying to lookup all ids as a batch: %s", ids)
    try:
        url = ncbi.efetch_fasta_url(ids)
        with wget.wget(url) as handle:
            yield from fasta.parse(handle)
    except:
        LOGGER.info("Failed to efetch all ids, will try individual lookup")
        for component in ids:
            LOGGER.debug("Fetching %s", component)
            yield from accession_records(info, component)
