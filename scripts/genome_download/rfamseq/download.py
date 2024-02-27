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

from attrs import Factory, define
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
            LOGGER.info("Fetching from ENA failed, falling back to NCBI")

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
            LOGGER.error("Selected no sequences")
            return
        elif selected == ids:
            LOGGER.debug("Using all sequences")
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
        LOGGER.info("Failed to fetch from NCBI")
        LOGGER.debug(err)

    try:
        for accession in accessions:
            LOGGER.info("Trying to fetch %s from ENA", accession)
            with ena.fetch_fasta(accession) as handle:
                yield from fasta.parse(handle)
                return
    except Exception as err:
        LOGGER.info("Failed to fetch from ENA")
        LOGGER.debug(err)

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
        except Exception:
            LOGGER.debug("Failed to lookup %s, will try suppressed", wgs_set)

        with ena.wgs_fasta(wgs_set, max_increase=0, use_suppressed=True) as handle:
            yield from fasta.parse(handle)
            continue


@define
class GenomeDownloader:
    """
    This class is handles downloading all sequences for a specific proteome
    from UniProt. This handles fetching the genomic sequences, filtering them
    and then fetching things from NCBI or ENA with fallbacks if needed. It then
    tracks all sequences found and produces some metadata files that can be
    used to update the Rfam database.

    :param SqliteDict info: An SqliteDict of the the assembly summary information
    as produced by ncbi.parse_assembly_files. This is assumed to cover all
    possible assemblies.
    :param ProteomeInfo proteome: The summary of the proteome to download.
    :param None | NcbiAssemblyReport: The assembly report, if any for this
    proteome.
    :param [FromFasta] _from_fasta: A list of all fasta entries seen when
    downloading the requested sequences. This is populated by the class itself.
    """

    info: SqliteDict
    proteome: uniprot.ProteomeInfo
    _assembly_report: ty.Optional[ncbi.NcbiAssemblyReport] = None
    _from_fasta: ty.List[FromFasta] = Factory(list)
    prefetched: ty.Optional[Path] = None

    @classmethod
    def build(
        cls,
        info: SqliteDict,
        proteome: uniprot.ProteomeInfo,
        prefetched=None,
    ) -> GenomeDownloader:
        """
        Create a new GenomeDownloader object. This will lookup
        """
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
            _assembly_report=assembly_report,
            _from_fasta=[],
            prefetched=prefetched,
        )

    def is_suppressed_proteome(self) -> bool:
        """
        Check if this genome is suppressed by looking at the assembly report.
        If no assembly report is assigned, this returns False.
        """
        assembly_summary = self.assembly_summary
        if not assembly_summary:
            return False
        return assembly_summary.is_suppressed

    @property
    def assembly_summary(self) -> None | ncbi.assembly_summary.NcbiAssemblySummary:
        """
        Get the assembly summary if any for this genome. This returns None if
        the proteome's genome does not have an accession. This will also return
        None if the accession is not a known genome. This assumes the genome
        accession is versioned.
        """
        accession = self.proteome.genome_info.accession
        if not accession:
            return None
        return self.info.get(accession, None)

    def __fetch_records(self) -> Records:
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
                self.info, genome, self._assembly_report, missing, self.prefetched
            )

        if missing:
            yield from missing_records(missing)

    def metadata(self, version: str) -> Metadata:
        """
        Build a metadata object of all seen records. This can only be called
        after `records`, otherwise there will be no records to build the
        metadata from.
        """
        return Metadata.build(
            version, self.proteome, self._assembly_report, self._from_fasta
        )

    def records(self) -> Records:
        """
        Download all individual sequences for this GenomeDownloader. This will
        fetch all sequences and produce an iterable of records. This will track
        all seen records and use them to build the metadata object.
        """
        seen = set()
        for record in self.__fetch_records():
            if not record.seq:
                raise ValueError(f"Empty sequence is invalid {record.id}")
            if record.id in seen:
                LOGGER.error("Somehow got duplicate record id %s", record.id)
                continue
            self._from_fasta.append(FromFasta.from_record(record))
            yield record
            seen.add(record.id)
