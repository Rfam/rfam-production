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

from rfamseq import ena, fasta, fasta_filter, ncbi, uniprot, wgs

from .accession_method import records as accession_records

LOGGER = logging.getLogger(__name__)


@contextmanager
def fetch_fasta(info: SqliteDict, genome: uniprot.GenomeInfo) -> ty.Iterator[ty.IO]:
    assert genome.accession, f"Missing genome accession"
    did_ena = False
    if genome.source is uniprot.GenomeSource.ENA:
        try:
            LOGGER.info("Trying to fetch %s from ENA", genome.accession)
            with ena.fetch_fasta(genome.accession) as handle:
                yield handle
                did_ena = True
                return
        except:
            LOGGER.info("Fetching from ENA failed")

    if did_ena:
        raise ValueError("impossible state")

    LOGGER.info("Fetching genome %s from NCBI", genome.accession)
    with ncbi.fetch_fasta(info, genome.accession) as handle:
        yield handle


def wgs_requested(
    genome: uniprot.GenomeInfo, ncbi_info: ncbi.NcbiAssemblyReport
) -> bool:
    if not ncbi_info.wgs_project:
        return False

    if not isinstance(genome.components, uniprot.SelectedComponents):
        return False

    start = ncbi_info.wgs_project[0:4]
    for accession in genome.components.accessions:
        if isinstance(accession, str) and start in accession:
            return True
    return False


def wgs_summary(
    ncbi_info: ncbi.NcbiAssemblyReport, genome: uniprot.GenomeInfo
) -> ty.Optional[wgs.WgsSummary]:
    if not ncbi_info.wgs_project or not wgs_requested(genome, ncbi_info):
        return None

    summary = wgs.resolve_wgs(wgs.WgsPrefix.build(ncbi_info.wgs_project))
    if not summary:
        LOGGER.info("Did not resolve WGS set %s", ncbi_info.wgs_project)
        return None

    LOGGER.info("Resolved WGS set %s", ncbi_info.wgs_project)
    LOGGER.debug("Found %s", summary)
    return summary


def records(
    info: SqliteDict, proteome: uniprot.ProteomeInfo, ncbi_info: ncbi.NcbiAssemblyReport
) -> ty.Iterator[SeqIO.SeqRecord]:
    genome = proteome.genome_info
    assert genome.accession, f"Missing genome accession in {proteome}"

    if not isinstance(genome.components, uniprot.SelectedComponents):
        raise ValueError(f"Invalid components type {genome}")

    comp_set = fasta_filter.ComponentSelector.from_selected(
        ncbi_info,
        genome.components,
        wgs_summary(ncbi_info, genome),
    )

    with fetch_fasta(info, genome) as records:
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
                yield from accession_records(info, classification.accession)

            elif isinstance(classification, fasta_filter.MissingWgsSet):
                LOGGER.info(
                    "Accession %s did not contain any known sequences from WGS set %s",
                    genome.accession,
                    classification.prefix,
                )
                with ena.wgs_fasta(classification.prefix) as handle:
                    yield from fasta.parse(handle)

            else:
                raise ValueError(
                    f"Impossible state with {classification} for {proteome}"
                )
