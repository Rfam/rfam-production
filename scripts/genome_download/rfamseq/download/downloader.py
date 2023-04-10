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

import enum
import logging
import typing as ty

from attrs import frozen
from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import ncbi, uniprot

from .accession_method import fetch as fetch_accessions
from .component_method import fetch as fetch_components
from .fallback_method import fetch as fallback_fetch

Records = ty.Iterable[SeqIO.SeqRecord]

LOGGER = logging.getLogger(__name__)


@enum.unique
class DownloadMethod(enum.Enum):
    BY_COMPONENT = "LOOKUP_EACH_COMPONENT"
    LOOKUP_GENOME_ACCESSION = "LOOKUP_GENOME_ACCESSION"
    FALLBACK = "LOOKUP_WITH_FALLBACK"

    def fetch(
        self,
        info: SqliteDict,
        proteome: uniprot.ProteomeInfo,
        ncbi_info: ty.Optional[ncbi.NcbiAssemblyReport] = None,
    ) -> Records:
        if self is DownloadMethod.BY_COMPONENT:
            yield from fetch_components(info, proteome)
        elif self is DownloadMethod.LOOKUP_GENOME_ACCESSION:
            assert isinstance(
                proteome.genome_info.accession, str
            ), "Invalid genome_info"
            yield from fetch_accessions(info, proteome.genome_info.accession)
        elif self is DownloadMethod.FALLBACK:
            assert isinstance(ncbi_info, ncbi.NcbiAssemblyReport), "Must give ncbi_info"
            yield from fallback_fetch(info, proteome, ncbi_info)
        else:
            raise ValueError("Cannot fetch with method %s", self)


@frozen
class GenomeDownload:
    method: DownloadMethod
    proteome: uniprot.ProteomeInfo
    assembly_info: ty.Optional[ncbi.NcbiAssemblyReport]

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
            ncbi_info = ncbi.fetch_assembly_report(info, genome.accession)
        except ncbi.UnknownGenomeId:
            return cls.by_components(proteome)

        # If we are fetching a genome which has a single component, and that
        # component is the WGS project assigned the genome in the NCBI assembly
        # info, we c
        return cls(DownloadMethod.FALLBACK, proteome, ncbi_info)

    def records(self, info: SqliteDict) -> Records:
        LOGGER.info("Fetching sequences method %s", self.method.name)
        yield from self.method.fetch(info, self.proteome, ncbi_info=self.assembly_info)
