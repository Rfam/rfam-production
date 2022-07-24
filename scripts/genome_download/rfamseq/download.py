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

import typing as ty
import logging

from attrs import define

from Bio import SeqIO
from sqlitedict import SqliteDict

from rfamseq import ncbi, ena, wget, uniprot, fasta_filter

LOGGER = logging.getLogger(__name__)

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


def lookup_sequences(info: SqliteDict, accession: str) -> Records:
    LOGGER.info("Trying to fetch %s from NCBI", accession)
    try:
        yield from ncbi.fetch_fasta(info, accession)
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


def sequences(info: SqliteDict, proteome: uniprot.ProteomeInfo) -> Records:
    genome = proteome.genome_info
    if genome.accession is None:
        LOGGER.info("Looking up each component for %s", proteome.upi)
        assert isinstance(
            genome.components, uniprot.SelectedComponents
        ), f"Invalid components for {proteome}"
        for component in genome.components:
            assert isinstance(component, str), f"Invalid component in {proteome}"
            LOGGER.info("Trying to lookup %s", component)
            yield from lookup_sequences(info, component)
    elif isinstance(genome.accession, str):
        ncbi_info = ncbi.assembly_info(info, genome.accession)

        if isinstance(genome.components, uniprot.All):
            LOGGER.info("Using all components in %s", genome.accession)
            yield from lookup_sequences(info, genome.accession)

        elif isinstance(genome.components, uniprot.SelectedComponents):
            comp_set = fasta_filter.ComponentSet.from_selected(
                ncbi_info.sequence_info, genome.components
            )
            records = lookup_sequences(info, genome.accession)
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
                    yield from lookup_sequences(info, classification.accession)
                else:
                    raise ValueError(
                        f"Impossible state with {classification} for {proteome}"
                    )
        else:
            raise ValueError("Impossible state")
    else:
        raise ValueError(f"Unknown type of accession {genome}")
