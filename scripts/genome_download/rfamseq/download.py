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

import enum
import logging
import re
import typing as ty

from attrs import define, frozen
from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import ena, fasta_filter, ncbi, uniprot, wget, wgs

LOGGER = logging.getLogger(__name__)

WGS_PATTERN = re.compile(r"^([A-Z]{4,6})\d{2,}")

Records = ty.Iterable[SeqIO.SeqRecord]


@define
class GenomeInfo:
    assembly_name: ty.Optional[str]
    study_ref: ty.Optional[str]
    description: str
    ciruclar: bool
    taxid: int
    kingdom: str


@define
class SequenceInfo:
    fasta_accession: str
    taxid: str
    mol_type: str
    description: str
    length: int
    source: str

    @property
    def accession(self) -> str:
        return self.fasta_accession.split(".", 1)[0]

    def version(self) -> str:
        parts = self.fasta_accession.split(".")
        version = "1"
        if len(parts) == 2:
            version = parts[1]
        return "%05i" % version


def lookup_fasta(info: SqliteDict, accession: str) -> Records:
    try:
        LOGGER.info("Trying to fetch %s from NCBI", accession)
        yield from ncbi.fetch_fasta(info, accession)
    except wget.FetchError:
        LOGGER.info("Failed to fetch %s from NCBI, trying ENA", accession)
        yield from ena.fetch_fasta(accession)


def lookup_by_accession(info: SqliteDict, accession: str) -> Records:
    LOGGER.info("Trying to fetch %s from NCBI", accession)
    try:
        yield from lookup_fasta(info, accession)
    except wget.FetchError:
        LOGGER.info("Trying to find contigs from %s", accession)
        contigs = list(ena.fetch_contigs(accession))
        for contig in contigs:
            LOGGER.info("Fetching contig %s from NCBI", contig)
            try:
                yield from ncbi.efetch_fasta(contig)
            except wget.FetchError:
                LOGGER.info("Failed to get %s from NCBI", contig)
                yield from ena.fetch_fasta(contig)


def sequences_by_components(
    info: SqliteDict, proteome: uniprot.ProteomeInfo
) -> Records:
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
        yield from ncbi.efetch_fasta(ids)
    except:
        LOGGER.info("Failed to efetch all ids, will try individual lookup")
        for component in genome.components:
            assert isinstance(component, str), f"Invalid component in {proteome}"
            yield from lookup_by_accession(info, component)


def wgs_sequences(wgs: str) -> Records:
    yield from ena.fetch_wgs_sequences(wgs)


def lookup_with_fallback(
    info: SqliteDict, proteome: uniprot.ProteomeInfo, ncbi_info: ncbi.NcbiAssemblyInfo
) -> Records:
    genome = proteome.genome_info
    assert genome.accession, f"Missing genome accession in {proteome}"

    wgs_accessions = None
    if ncbi_info.wgs_project:
        resolved = wgs.resolve_wgs(ncbi_info.wgs_project)
        if resolved:
            LOGGER.info("Resolved WGS set %s", ncbi_info.wgs_project)
            LOGGER.debug("Found %s", resolved)
            wgs_accessions = {ncbi_info.wgs_project: resolved.ids()}
        else:
            LOGGER.info("Did not resolve WGS set %s", ncbi_info.wgs_project)

    assert isinstance(
        genome.components, uniprot.SelectedComponents
    ), f"Invalid components type {genome}"
    comp_set = fasta_filter.ComponentSet.from_selected(
        ncbi_info.sequence_info,
        genome.components,
        wgs_accessions,
    )

    records = lookup_by_accession(info, genome.accession)

    for classification in fasta_filter.filter(records, comp_set):
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
        elif isinstance(classification, fasta_filter.Missing):
            LOGGER.info(
                "Accession %s did not contain %s",
                genome.accession,
                classification.accession,
            )
            missing = classification.accession
            if match := re.match(WGS_PATTERN, missing):
                LOGGER.info("Accession %s looks like a wgs set", missing)
                try:
                    yield from ena.fetch_wgs_sequences(match.group(1))
                    continue
                except:
                    LOGGER.info("Failed to download WGS sequences")
            yield from lookup_by_accession(info, classification.accession)
        else:
            raise ValueError(f"Impossible state with {classification} for {proteome}")


@enum.unique
class DownloadMethod(enum.Enum):
    BY_COMPONENT = "LOOKUP_EACH_COMPONENT"
    LOOKUP_GENOME_ACCESSION = "LOOKUP_GENOME_ACCESSION"
    FALLBACK = "LOOKUP_WITH_FALLBACK"

    def fetch(
        self,
        info: SqliteDict,
        proteome: uniprot.ProteomeInfo,
        ncbi_info: ty.Optional[ncbi.NcbiAssemblyInfo] = None,
    ) -> Records:
        if self is DownloadMethod.BY_COMPONENT:
            yield from sequences_by_components(info, proteome)
        elif self is DownloadMethod.LOOKUP_GENOME_ACCESSION:
            assert isinstance(
                proteome.genome_info.accession, str
            ), "Invalid genome_info"
            yield from lookup_by_accession(info, proteome.genome_info.accession)
        elif self is DownloadMethod.FALLBACK:
            assert isinstance(ncbi_info, ncbi.NcbiAssemblyInfo), "Must give ncbi_info"
            yield from lookup_with_fallback(info, proteome, ncbi_info)
        else:
            raise ValueError("Cannot fetch with method %s", self)


@frozen
class GenomeDownload:
    method: DownloadMethod
    proteome: uniprot.ProteomeInfo
    assembly_info: ty.Optional[ncbi.NcbiAssemblyInfo]

    @classmethod
    def by_components(cls, proteome: uniprot.ProteomeInfo) -> GenomeDownload:
        LOGGER.info("Fetching proteome %s by listed components", proteome)
        return cls(DownloadMethod.BY_COMPONENT, proteome, None)

    @classmethod
    def by_lookup(cls, proteome: uniprot.ProteomeInfo) -> GenomeDownload:
        LOGGER.info("Using all sequences listed in %s", proteome.genome_info.accession)
        return cls(DownloadMethod.LOOKUP_GENOME_ACCESSION, proteome, None)

    @classmethod
    def build(cls, info: SqliteDict, proteome: uniprot.ProteomeInfo) -> GenomeDownload:
        genome = proteome.genome_info
        if genome.accession is None:
            return cls.by_components(proteome)

        if not isinstance(genome.accession, str):
            raise ValueError(f"Unknown type of genome accession {genome}")

        LOGGER.info("Extracting based on genome %s", genome.accession)
        if isinstance(genome.components, uniprot.All):
            return cls.by_lookup(proteome)

        if not isinstance(genome.components, uniprot.SelectedComponents):
            raise ValueError(f"Unknown type of components {genome.components}")

        LOGGER.info("Extracting selected components for %s", genome.accession)

        try:
            ncbi_info = ncbi.assembly_info(info, genome.accession)
        except ncbi.UnknownGenomeId:
            return cls.by_components(proteome)

        return cls(DownloadMethod.FALLBACK, proteome, ncbi_info)

    def records(self, info: SqliteDict) -> Records:
        LOGGER.info("Fetching sequences method %s", self.method.name)
        yield from self.method.fetch(info, self.proteome, ncbi_info=self.assembly_info)


# def sequences(info: SqliteDict, proteome: uniprot.ProteomeInfo) -> Records:
#     genome = proteome.genome_info
#     if genome.accession is None:
#         yield from sequences_by_components(info, proteome)
#     elif isinstance(genome.accession, str):
#         LOGGER.info("Extracting based on genome %s", genome.accession)
#         if isinstance(genome.components, uniprot.All):
#             LOGGER.info("Using all components in %s", genome.accession)
#             yield from lookup_sequences(info, genome.accession)

#         elif isinstance(genome.components, uniprot.SelectedComponents):
#             LOGGER.info("Extracting selected components for %s", genome.accession)

#             try:
#                 ncbi_info = ncbi.assembly_info(info, genome.accession)
#             except ncbi.UnknownGenomeId:
#                 yield from sequences_by_components(info, proteome)
#                 return

#             wgs_accessions = None
#             if ncbi_info.wgs_project:
#                 resolved = ncbi.resolve_wgs(ncbi_info.wgs_project)
#                 if resolved:
#                     wgs_accessions = {ncbi_info.wgs_project: resolved}

#             comp_set = fasta_filter.ComponentSet.from_selected(
#                 ncbi_info.sequence_info,
#                 genome.components,
#                 wgs_accessions,
#             )
#             records = lookup_sequences(info, genome.accession)

#             for classification in fasta_filter.filter(records, comp_set):
#                 if isinstance(classification, fasta_filter.Extra):
#                     LOGGER.info(
#                         "Accession %s contains extra record %s",
#                         genome.accession,
#                         classification.extra.id,
#                     )
#                     continue
#                 elif isinstance(classification, fasta_filter.Found):
#                     LOGGER.info(
#                         "Found expected record %s in %s",
#                         classification.record.id,
#                         genome.accession,
#                     )
#                     yield classification.record
#                 elif isinstance(classification, fasta_filter.Missing):
#                     LOGGER.info(
#                         "Accession %s did not contain %s",
#                         genome.accession,
#                         classification.accession,
#                     )
#                     missing = classification.accession
#                     if match := re.match(WGS_PATTERN, missing):
#                         LOGGER.info("Accession %s looks like a wgs set", missing)
#                         try:
#                             yield from ena.fetch_wgs_sequences(match.group(1))
#                             continue
#                         except:
#                             LOGGER.info("Failed to download WGS sequences")
#                     yield from lookup_sequences(info, classification.accession)
#                 else:
#                     raise ValueError(
#                         f"Impossible state with {classification} for {proteome}"
#                     )
#         else:
#             raise ValueError("Impossible state")
#     else:
#         raise ValueError(f"Unknown type of accession {genome}")
