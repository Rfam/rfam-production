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
import re
import typing as ty

from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import ena, fasta_filter, ncbi, uniprot, wgs

from .accession_method import fetch as fetch_accessions

Records = ty.Iterable[SeqIO.SeqRecord]

LOGGER = logging.getLogger(__name__)


def fetch_ena_genome(accession: str) -> Records:
    try:
        yield from ena.fetch_fasta(accession)
    except Exception:
        LOGGER.info("Failed to fetch from ENA")
        return None


def fetch_fasta(info: SqliteDict, genome: uniprot.GenomeInfo) -> Records:
    assert genome.accession, f"Missing genome accession"
    if genome.source is uniprot.GenomeSource.ENA:
        records = fetch_ena_genome(genome.accession)
        if records:
            LOGGER.info("Fetching genome %s from ENA", genome.accession)
            yield from records
            return
    LOGGER.info("Fetching genome %s from NCBI", genome.accession)
    yield from fetch_accessions(info, genome.accession)


def is_wgs(genome: uniprot.GenomeInfo, ncbi_info: ncbi.NcbiAssemblyReport) -> bool:
    found = (
        ncbi_info.wgs_project
        and isinstance(genome.components, uniprot.SelectedComponents)
        and any(
            ncbi_info.wgs_project[0:4] in acc for acc in genome.components.accessions
        )
    )
    return bool(found)


def fetch(
    info: SqliteDict, proteome: uniprot.ProteomeInfo, ncbi_info: ncbi.NcbiAssemblyReport
) -> Records:
    genome = proteome.genome_info
    assert genome.accession, f"Missing genome accession in {proteome}"
    wgs_accessions = None

    try:
        if ncbi_info.wgs_project and is_wgs(genome, ncbi_info):
            wgs_accessions = wgs.resolve_wgs(ncbi_info.wgs_project)
            if wgs_accessions:
                LOGGER.info("Resolved WGS set %s", ncbi_info.wgs_project)
                LOGGER.debug("Found %s", wgs_accessions)
            else:
                LOGGER.info("Did not resolve WGS set %s", ncbi_info.wgs_project)
    except TypeError as e:
        LOGGER.debug(
            "Catching TypeError to ignore type Unplaced items in accession list", e
        )

    assert isinstance(
        genome.components, uniprot.SelectedComponents
    ), f"Invalid components type {genome}"
    comp_set = fasta_filter.ComponentSelector.from_selected(
        ncbi_info,
        genome.components,
        wgs_accessions,
    )

    records = fetch_fasta(info, genome)
    for classification in comp_set.filter(records):
        if isinstance(classification, fasta_filter.Extra):
            LOGGER.info(
                "Accession %s contains extra record %s",
                genome.accession,
                classification.extra.id,
            )
            continue
        elif isinstance(classification, fasta_filter.Found):
            LOGGER.info(
                "Found expected record %s in %s",
                classification.record.id,
                genome.accession,
            )
            yield classification.record
        elif isinstance(classification, fasta_filter.MissingAccession):
            LOGGER.info(
                "Accession %s did not contain %s",
                genome.accession,
                classification.accession,
            )
            yield from fetch_accessions(info, classification.accession)
        elif isinstance(classification, fasta_filter.MissingWgsSet):
            LOGGER.info(
                "Accession %s did not contain any known sequences from WGS set %s",
                genome.accession,
                classification.prefix,
            )
            yield from ena.fetch_wgs_sequences(classification.prefix)
        else:
            raise ValueError(f"Impossible state with {classification} for {proteome}")
