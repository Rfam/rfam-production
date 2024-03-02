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

import json
import logging
import typing as ty
from contextlib import contextmanager

from attrs import frozen
from Bio import SeqIO
from boltons import iterutils
from sqlitedict import SqliteDict

from rfamseq import ena, fasta, ncbi, uniprot, wget
from rfamseq.accession import Accession
from rfamseq.metadata import FromFasta, Metadata
from rfamseq.missing import Missing
from rfamseq.utils import assert_never, batched

Records = ty.Iterable[SeqIO.SeqRecord]

LOGGER = logging.getLogger(__name__)


class NoRecordsFetched(Exception):
    """Rasied when no records are fetched"""


class NoComponentSource(Exception):
    """Raised if there are no sources for a given component"""


class TooManyComponentSources(Exception):
    """Raised if there is more than one source for a given component"""


class NoPossibleSequences(Exception):
    """Raised if there is nothing to fetch"""


class SuppressedGenome(Exception):
    """Raised if the given genome is suppressed by NCBI"""


class AccessionLookupFailed(Exception):
    """
    Raised when if all accession lookup methods failed to fetch some id.
    """


# @contextmanager
# def fetch_genome(
#     info: SqliteDict, genome: uniprot.GenomeInfo, fetched: ty.Optional[Path]
# ) -> ty.Iterator[ty.IO]:
#     assert genome.accession, "Genome must have a primary accession"
#
#     LOGGER.info("Fetching the requested genome %s", genome)
#     if fetched and fetched.exists():
#         LOGGER.info("Using given file %s", fetched)
#         with fetched.open("rb") as raw:
#             yield raw
#         return
#
#     if genome.source and genome.source.from_ebi():
#         try:
#             LOGGER.info("Trying to fetch %s from ENA", genome.accession)
#             with ena.fetch_fasta(genome.accession) as handle:
#                 yield handle
#             return
#         except Exception as err:
#             LOGGER.info("Fetching from ENA failed, falling back to NCBI")
#
#     LOGGER.info("Fetching genome %s from NCBI", genome.accession)
#     with ncbi.fetch_fasta(info, genome.accession) as handle:
#         yield handle
#
#
# def genomic_records(
#     info: SqliteDict,
#     genome: uniprot.GenomeInfo,
#     assembly_report: ty.Optional[ncbi.NcbiAssemblyReport],
#     missing: Missing,
#     fetched: ty.Optional[Path],
# ) -> Records:
#     assert genome.accession, "Genome must have a primary accession"
#
#     selector = FastaFilter.from_selected(assembly_report, genome.components)
#     with fetch_genome(info, genome, fetched) as handle:
#         ids = fasta.extract_ids(Path(handle.name))
#         selected, missed = selector.filter_ids(ids)
#         missing.update(missed)
#         if not selected:
#             LOGGER.error("Selected no sequences")
#             return
#         elif selected == ids:
#             LOGGER.debug("Using all sequences")
#             yield from SeqIO.parse(handle.name, "fasta")
#         else:
#             LOGGER.debug("Filtering with easel")
#             with easel.filtered(Path(handle.name), selected) as filtered:
#                 yield from SeqIO.parse(filtered.name, "fasta")
#
#
# @sleep_and_retry
# @limits(calls=3, period=1)
# def accession_fetch(accessions: ty.List[str]) -> Records:
#     try:
#         LOGGER.info("Trying to fetch %s from NCBI", accessions)
#         accession = ",".join(accessions)
#         url = ncbi.efetch_fasta_url(accession)
#         with wget.wget(url) as handle:
#             yield from fasta.parse(handle)
#             return
#     except Exception as err:
#         LOGGER.info("Failed to fetch from NCBI")
#         LOGGER.debug(err)
#
#     try:
#         for accession in accessions:
#             LOGGER.info("Trying to fetch %s from ENA", accession)
#             with ena.fetch_fasta(accession) as handle:
#                 yield from fasta.parse(handle)
#                 return
#     except Exception as err:
#         LOGGER.info("Failed to fetch from ENA")
#         LOGGER.debug(err)
#
#     raise AccessionLookupFailed(f"Failed to lookup {accessions}")
#
#
# def missing_records(missing: Missing) -> Records:
#     LOGGER.info("Will fetch missing records %s", missing)
#
#     for chunk in batched(missing.accessions, 3):
#         ids = [str(c) for c in chunk]
#         yield from accession_fetch(ids)
#
#     for chunk in batched(missing.wgs_sequences, 5):
#         ids = [c.to_wgs_string() for c in chunk]
#         yield from accession_fetch(ids)
#
#     for wgs_set in missing.wgs_sets:
#         try:
#             with ena.wgs_fasta(wgs_set) as handle:
#                 yield from fasta.parse(handle)
#                 continue
#         except Exception:
#             LOGGER.debug("Failed to lookup %s, will try suppressed", wgs_set)
#
#         with ena.wgs_fasta(wgs_set, max_increase=0, use_suppressed=True) as handle:
#             yield from fasta.parse(handle)
#             continue


@frozen
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
    """

    version: str
    ncbi: ncbi.FtpWrapper
    ena: ena.EnaWrapper

    @classmethod
    def build(cls, version: str, info: SqliteDict) -> GenomeDownloader:
        """Create a new GenomeDownloader object."""
        return cls(
            version=version,
            ncbi=ncbi.FtpWrapper.build(info),
            ena=ena.EnaWrapper.build(),
        )

    def __ncbi_fetch__(self, accessions: ty.List[Accession], batch_size=3) -> Records:
        for batch in iterutils.chunked_iter(accessions, batch_size):
            LOGGER.info("NCBI Fetching accessions %s", batch)
            url = self.ncbi.efetch_url(batch)
            with wget.wget(url) as handle:
                yield from fasta.parse(handle)

    def __ebi_fetch__(self, accessions: ty.List[Accession], **kwargs) -> Records:
        for accession in accessions:
            LOGGER.info("ENA Fetching accession %s", accession)
            url = self.ena.fasta_url(accession)
            with wget.wget(url) as handle:
                yield from fasta.parse(handle)

    def __fetch_accessions__(
        self, accessions: ty.List[Accession], batch_size=3
    ) -> Records:
        assert accessions, "Must give at least one accession to fetch"
        errors = []
        LOGGER.info("Fetching components %s", ", ".join(str(a) for a in accessions))
        for fetcher in [self.__ncbi_fetch__, self.__ebi_fetch__]:
            LOGGER.debug("Trying %s", fetcher.__name__)
            try:
                # TODO Validate this finds all expected accessions. The issue
                # is that the requested accession and the found accessions may
                # or may not have versions, so the comparision should be
                # versionless. Additionally, if this is looking up something
                # like a genome or wgs set, somehow it may provide many
                # accessions, which do not easily match the requested one but
                # are correct.
                yield from fetcher(accessions, batch_size=batch_size)
                break
            except Exception as err:
                LOGGER.debug("Could not use %s for %s", fetcher.__name__, accessions)
                errors.append(err)
        else:
            for err in errors:
                LOGGER.exception(err)
            raise AccessionLookupFailed(f"Could not fetch {accessions}")

    def __fetch_genome__(self, accession: Accession) -> Records:
        assert accession.is_genomic()
        if assembly_summary := self.ncbi.assembly_summary(accession):
            if not assembly_summary.is_latest:
                if latest := self.ncbi.find_latest_version(accession):
                    LOGGER.info(
                        "Accession %s is not the latest, using newer: %s",
                        accession,
                        latest,
                    )
                    accession = latest
                else:
                    raise SuppressedGenome(accession)
        else:
            LOGGER.warning("No assembly summary for %s, is it real", str(accession))

        urls: ty.List[str] = []
        if url := self.ncbi.genome_url(accession):
            urls.append(url)
        else:
            LOGGER.debug("Could not get NCBI url for %s", str(accession))

        urls.append(self.ena.fasta_url(accession))

        errors = []
        for url in urls:
            try:
                with wget.wget(url) as tmp:
                    yield from fasta.parse(tmp)
                break
            except Exception as err:
                LOGGER.debug("Failed to get %s", url)
                errors.append(err)
                continue
        else:
            for err in errors:
                LOGGER.exception(err)
            raise AccessionLookupFailed(f"Could not get {accession}")

    def __fetch_proteome_records__(
        self, proteome: uniprot.proteome.Proteome
    ) -> Records:
        """Fetch all records for the given proteome."""

        genome = proteome.genome_assembly
        accessions: ty.List[Accession] = []
        if genome.assembly_id:
            accessions.append(Accession.build(genome.assembly_id))
        else:
            for component in proteome.components:
                if component.name.lower() == "unplaced":
                    continue

                if not component.proteome_cross_references:
                    raise NoComponentSource(component)

                component_sources = []
                for xref in component.proteome_cross_references:
                    match xref.database:
                        case uniprot.proteome.ComponentSource.GENOME_ACCESSION:
                            component_sources.append(Accession.build(xref.id))
                        case uniprot.proteome.ComponentSource.BIOSAMPLE:
                            pass
                        case _:
                            assert_never(xref.database)
                if not component_sources:
                    raise NoComponentSource(component)
                accessions.extend(component_sources)

        if not accessions:
            raise NoPossibleSequences(genome)

        if len(accessions) == 1 and accessions[0].is_genomic():
            LOGGER.info("Fetching via complete genome %s", str(accessions[0]))
            yield from self.__fetch_genome__(accessions[0])
        else:
            for accession in accessions:
                if accession.is_genomic():
                    LOGGER.error("Multiple accessions/genome %s", str(accession))
                    raise ValueError("Cannot fetch many accessions along side a genome")
            LOGGER.info(
                "Fetching via components %s", ", ".join(str(a) for a in accessions)
            )
            yield from self.__fetch_accessions__(accessions)

    def download_proteome(
        self, proteome: uniprot.proteome.Proteome, fasta: ty.IO, metadata: ty.IO
    ):
        """Download all individual sequences for the given proteome. This will
        fetch all sequences and write each entry to the given fasta handle. A
        metadata object of all sequences downloaded will be created and written
        to the metadata handle. It is an error if no sequences are downloaded.
        """

        from_fasta = []
        for record in self.__fetch_proteome_records__(proteome):
            from_fasta.append(FromFasta.from_record(record))
            SeqIO.write(record, fasta, "fasta")

        if not from_fasta:
            raise NoRecordsFetched(f"Found no records for {proteome.id}")

        # assembly_report = None
        # if genome.assembly_id:
        #     assembly_report = self.ncbi.assembly_summary(genome.assembly_id)
        # meta = Metadata.build(self.version, proteome, assembly_report, from_fasta)
        # json.dump(cattrs.unstructure(meta), metadata)
