# -*- coding: utf-8 -*-

"""
Copyright [2009-2023] EMBL-European Bioinformatics Institute
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

from attrs import define
from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import ena, fasta, ncbi, uniprot, wget, wgs
from rfamseq.accession import Accession
from rfamseq.fasta_filter import (
    FastaFilter,
    Found,
    MissingAccession,
    MissingWgsSequence,
    MissingWgsSet,
)
from rfamseq.utils import assert_never, batched

Records = ty.Iterable[SeqIO.SeqRecord]

LOGGER = logging.getLogger(__name__)


class AccessionLookupFailed(Exception):
    """
    Raised when if all accession lookup methods failed to fetch some id.
    """


@define
class Missing:
    accessions: ty.Set[Accession]
    wgs_sets: ty.Set[wgs.WgsPrefix]
    wgs_sequences: ty.Set[wgs.WgsSequenceId]

    @classmethod
    def empty(cls) -> Missing:
        return Missing(accessions=set(), wgs_sets=set(), wgs_sequences=set())

    def add_components(self, comps: uniprot.SelectedComponents):
        if comps.unplaced:
            raise ValueError("Cannot add unplaced to Missing")
        self.accessions.update(comps.accessions)
        self.wgs_sets.update(comps.wgs_sets)
        self.wgs_sequences.update(comps.wgs_sequences)


@contextmanager
def fetch_genome(info: SqliteDict, genome: uniprot.GenomeInfo) -> ty.Iterator[ty.IO]:
    assert genome.accession, "Genome must have a primary accession"

    if genome.source is uniprot.GenomeSource.ENA:
        try:
            LOGGER.info("Trying to fetch %s from ENA", genome.accession)
            with ena.fetch_fasta(genome.accession) as handle:
                yield handle
            return
        except Exception as err:
            LOGGER.info("Fetching from ENA failed")
            LOGGER.exception(err)

    LOGGER.info("Fetching genome %s from NCBI", genome.accession)
    with ncbi.fetch_fasta(info, genome.accession) as handle:
        yield handle


def genomic_records(
    info: SqliteDict,
    genome: uniprot.GenomeInfo,
    assembly_report: ty.Optional[ncbi.NcbiAssemblyReport],
    missing: Missing,
) -> Records:
    assert genome.accession, "Genome must have a primary accession"

    with fetch_genome(info, genome) as records:
        selector = FastaFilter.from_selected(assembly_report, genome.components)
        for classification in selector.filter(records):
            match classification:
                case Found():
                    LOGGER.info(
                        "Found expected record %s in %s",
                        classification.record.id,
                        genome.accession,
                    )
                    yield classification.record

                case MissingAccession():
                    LOGGER.info(
                        "Accession %s did not contain %s",
                        genome.accession,
                        classification.accession,
                    )
                    missing.accessions.add(classification.accession)

                case MissingWgsSet():
                    LOGGER.info(
                        "Genome %s did not contain any known sequences from WGS set %s",
                        genome.accession,
                        classification.prefix,
                    )
                    missing.wgs_sets.add(classification.prefix)

                case MissingWgsSequence():
                    LOGGER.info(
                        "Genome %s did not contain WGS sequence %s",
                        genome.accession,
                        classification.sequence_id,
                    )
                    missing.wgs_sequences.add(classification.sequence_id)

                case _:
                    assert_never(classification)


def accession_fetch(accessions: ty.List[str]) -> Records:
    try:
        LOGGER.info("Trying to fetch %s from NCBI", accessions)
        accession = ",".join(accessions)
        url = ncbi.efetch_fasta_url(accession)
        with wget.wget(url) as handle:
            yield from fasta.parse(handle)
            return
    except Exception as err:
        LOGGER.debug(err)
        LOGGER.debug("Failed to fetch from NCBI")

    try:
        for accession in accessions:
            LOGGER.info("Trying to fetch %s from ENA", accession)
            with ena.fetch_fasta(accession) as handle:
                yield from fasta.parse(handle)
                return
    except Exception as err:
        LOGGER.debug(err)
        LOGGER.debug("Failed to fetch from ENA")

    raise AccessionLookupFailed(f"Failed to lookup {accessions}")


def missing_records(missing: Missing) -> Records:
    for chunk in batched(missing.accessions, 3):
        ids = [str(c) for c in chunk]
        yield from accession_fetch(ids)

    for chunk in batched(missing.wgs_sequences, 5):
        ids = [c.to_wgs_string() for c in chunk]
        yield from accession_fetch(ids)

    for wgs_set in missing.wgs_sets:
        with ena.wgs_fasta(wgs_set) as handle:
            yield from fasta.parse(handle)


@define
class GenomeDownloader:
    info: SqliteDict
    proteome: uniprot.ProteomeInfo
    assembly_report: ty.Optional[ncbi.NcbiAssemblyReport]

    @classmethod
    def build(
        cls, info: SqliteDict, proteome: uniprot.ProteomeInfo
    ) -> GenomeDownloader:
        genome = proteome.genome_info
        assembly_report = None
        try:
            if genome.accession:
                assembly_report = ncbi.fetch_assembly_report(info, genome.accession)
        except Exception as err:
            LOGGER.debug("Failed to load assembly info for %s", genome.accession)
            LOGGER.debug(err)
        return cls(
            info=info,
            proteome=proteome,
            assembly_report=assembly_report,
        )

    def fetch_records(self) -> Records:
        genome = self.proteome.genome_info
        missing = Missing.empty()
        if not genome.accession:
            match genome.components:
                case uniprot.All():
                    raise ValueError("Must have a primary accession of all requested")
                case uniprot.SelectedComponents():
                    missing.add_components(genome.components)
                case _:
                    assert_never(genome.components)
        else:
            yield from genomic_records(self.info, genome, self.assembly_report, missing)

        yield from missing_records(missing)

    def records(self) -> Records:
        seen = set()
        for record in self.fetch_records():
            if not record.seq:
                raise ValueError(f"Empty sequence is invalid {record.id}")
            if record.id in seen:
                LOGGER.error("Somehow got duplicate record id %s", record.id)
                continue
            yield record
            seen.add(record.id)
