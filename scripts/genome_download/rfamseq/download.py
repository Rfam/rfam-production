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
from pathlib import Path

from attrs import define
from Bio import SeqIO
from ratelimit import limits, sleep_and_retry
from sqlitedict import SqliteDict

from rfamseq import easel, ena, fasta, ncbi, uniprot, wget
from rfamseq.fasta_filter import FastaFilter
from rfamseq.metadata import FromFasta, Metadata
from rfamseq.missing import Missing
from rfamseq.utils import assert_never, batched

Records = ty.Iterable[SeqIO.SeqRecord]

LOGGER = logging.getLogger(__name__)


class AccessionLookupFailed(Exception):
    """
    Raised when if all accession lookup methods failed to fetch some id.
    """


@contextmanager
def fetch_genome(
    info: SqliteDict, genome: uniprot.GenomeInfo, fetched: ty.Optional[Path]
) -> ty.Iterator[ty.IO]:
    assert genome.accession, "Genome must have a primary accession"

    LOGGER.info("Fetching the requested genome %s", genome)
    if fetched and fetched.exists():
        LOGGER.info("Using given file %s", fetched)
        with fetched.open("rb") as raw:
            yield raw
        return

    if genome.source and genome.source.from_ebi():
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
    fetched: ty.Optional[Path],
) -> Records:
    assert genome.accession, "Genome must have a primary accession"

    selector = FastaFilter.from_selected(assembly_report, genome.components)
    with fetch_genome(info, genome, fetched) as handle:
        ids = fasta.extract_ids(Path(handle.name))
        selected, missed = selector.filter_ids(ids)
        missing.update(missed)
        if not selected:
            raise ValueError("Selected no sequences")
        elif selected == ids:
            yield from SeqIO.parse(handle.name, "fasta")
        else:
            LOGGER.debug("Filtering with easel")
            with easel.filtered(Path(handle.name), selected) as filtered:
                yield from SeqIO.parse(filtered.name, "fasta")


@sleep_and_retry
@limits(calls=3, period=1)
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
    LOGGER.info("Will fetch missing records %s", missing)

    for chunk in batched(missing.accessions, 3):
        ids = [str(c) for c in chunk]
        yield from accession_fetch(ids)

    for chunk in batched(missing.wgs_sequences, 5):
        ids = [c.to_wgs_string() for c in chunk]
        yield from accession_fetch(ids)

    for wgs_set in missing.wgs_sets:
        try:
            with ena.wgs_fasta(wgs_set) as handle:
                yield from fasta.parse(handle)
                continue
        except wget.FetchError as err:
            LOGGER.debug("Failed to lookup %s, will try suppressed", wgs_set)

        with ena.wgs_fasta(wgs_set, max_increase=0, use_suppressed=True) as handle:
            yield from fasta.parse(handle)
            continue


@define
class GenomeDownloader:
    info: SqliteDict
    proteome: uniprot.ProteomeInfo
    assembly_report: ty.Optional[ncbi.NcbiAssemblyReport]
    from_fasta: ty.List[FromFasta]
    prefetched: ty.Optional[Path] = None

    @classmethod
    def build(
        cls,
        info: SqliteDict,
        proteome: uniprot.ProteomeInfo,
        prefetched=None,
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
            from_fasta=[],
            prefetched=prefetched,
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
            yield from genomic_records(
                self.info, genome, self.assembly_report, missing, self.prefetched
            )

        if missing:
            yield from missing_records(missing)

    def metadata(self, version: str) -> Metadata:
        return Metadata.build(
            version, self.proteome, self.assembly_report, self.from_fasta
        )

    def records(self) -> Records:
        seen = set()
        for record in self.fetch_records():
            if not record.seq:
                raise ValueError(f"Empty sequence is invalid {record.id}")
            if record.id in seen:
                LOGGER.error("Somehow got duplicate record id %s", record.id)
                continue
            self.from_fasta.append(FromFasta.from_record(record))
            yield record
            seen.add(record.id)
