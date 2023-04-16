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

import collections as coll
import logging
import typing as ty

from attrs import define, frozen
from Bio import SeqIO

from rfamseq import fasta, ncbi, uniprot, wgs
from rfamseq.accession import Accession

LOGGER = logging.getLogger(__name__)


@frozen
class MissingAccession:
    accession: str


@frozen
class MissingWgsSet:
    prefix: wgs.WgsPrefix


@frozen
class Found:
    matching_accession: str
    record: SeqIO.SeqRecord


@frozen
class Extra:
    extra: SeqIO.SeqRecord


Missing = ty.Union[MissingAccession, MissingWgsSet]
RecordTypes = ty.Union[Missing, Found, Extra]


@define
class WgsSequenceContainer:
    """
    This represents a way to refer to a collection of WGS sets.
    """

    any_member_of: ty.Set[wgs.WgsPrefix]
    ranges: ty.Set[wgs.WgsSequenceRange]

    @classmethod
    def empty(cls) -> WgsSequenceContainer:
        return cls(any_member_of=set(), ranges=set())

    def __add_sequence_range(self, range: wgs.WgsSequenceRange):
        if range.is_single_range() and range.start.is_wgs_set_reference():
            self.any_member_of.add(range.prefix)
        else:
            self.ranges.add(range)
        return None

    def add(
        self,
        accession: ty.Union[wgs.WgsPrefix, wgs.WgsSequenceId, wgs.WgsSequenceRange],
    ):
        if isinstance(accession, wgs.WgsPrefix):
            self.any_member_of.add(accession)
            return None
        elif isinstance(accession, wgs.WgsSequenceRange):
            return self.__add_sequence_range(accession)
        elif isinstance(accession, wgs.WgsSequenceId):
            seq_range = wgs.WgsSequenceRange.from_endpoint(
                wgs.WgsSequenceKind.SEQUENCE, accession
            )
            return self.__add_sequence_range(seq_range)
        raise ValueError(f"Cannot treat {accession} as a WGS accession")

    def is_empty(self) -> bool:
        return not (bool(self.ranges) or bool(self.any_member_of))

    def matching_sets(
        self,
        accession: ty.Union[str, Accession, wgs.WgsSequenceId],
        within_one_version=False,
    ) -> ty.Optional[ty.Set[wgs.WgsPrefix]]:
        if isinstance(accession, Accession):
            accession = wgs.WgsSequenceId.build(str(accession))
        if isinstance(accession, str):
            accession = wgs.WgsSequenceId.build(accession)

        for prefix in self.any_member_of:
            if prefix.matches(accession.prefix, within_one_version=within_one_version):
                return set([accession.prefix, prefix])

        for range in self.ranges:
            if range.includes(accession, within_one_version=within_one_version):
                return set([range.prefix])
        return None


@frozen
class RequestedAccessions:
    unplaced: bool
    standard_accessions: ty.Set[Accession]
    wgs_set: WgsSequenceContainer

    @classmethod
    def build(cls, selected: uniprot.SelectedComponents) -> RequestedAccessions:
        accessions = set()
        unplaced = False
        wgs_set = WgsSequenceContainer.empty()
        for accession in selected.accessions:
            if isinstance(accession, uniprot.Unplaced):
                unplaced = True
                continue
            elif isinstance(accession, (wgs.WgsPrefix, wgs.WgsSequenceId)):
                wgs_set.add(accession)
            else:
                accessions.add(accession)

        return cls(
            unplaced=unplaced,
            standard_accessions=accessions,
            wgs_set=wgs_set,
        )

    def includes_wgs(self):
        return not self.wgs_set.is_empty()

    def includes_unplaced(self):
        return self.unplaced


@define
class SeenAccessions:
    accessions: ty.Set[Accession]
    wgs_prefix: ty.Set[wgs.WgsPrefix]

    @classmethod
    def empty(cls) -> SeenAccessions:
        return cls(accessions=set(), wgs_prefix=set())

    def mark_wgs_prefixes(self, wgs_ids: ty.Set[wgs.WgsPrefix]):
        self.wgs_prefix.update(wgs_ids)

    def mark_accession(self, accession: Accession):
        self.accessions.add(accession)

    def seen_wgs(self) -> bool:
        return bool(self.wgs_prefix)

    def __contains__(self, accession: Accession) -> bool:
        return any(a.matches(accession) for a in self.accessions)


@frozen
class ComponentSelector:
    requested: RequestedAccessions
    assembly_report: ncbi.NcbiAssemblyReport
    summaries: ty.Dict[wgs.WgsSummary, wgs.WgsSummary]

    @classmethod
    def from_selected(
        cls,
        assembly_report: ncbi.NcbiAssemblyReport,
        selected: uniprot.SelectedComponents,
        wgs_summary: ty.Optional[wgs.WgsSummary],
    ) -> ComponentSelector:
        summaries = {}
        if wgs_summary:
            summaries[wgs_summary] = wgs_summary

        return cls(
            requested=RequestedAccessions.build(selected),
            assembly_report=assembly_report,
            summaries=summaries,
        )

    def matching_wgs_set(self, id: str) -> ty.Optional[ty.Set[wgs.WgsPrefix]]:
        if not wgs.looks_like_wgs_accession(id):
            return None

        # Check if this sequence is from a WGS set. This means the id is any one
        # of the ones known to be a part of that WGS set. We will treat seeing any
        # one sequence from a wgs set as seeing the whole WGS set. Slightly
        # risky but should be fine. Maybe.
        #
        # This should take the accession from the sequence and check if it
        # looks like any sequence id know about for the given WGS set
        if prefix := self.requested.wgs_set.matching_sets(id, within_one_version=True):
            return prefix
        return None

    def matching_accession(self, id: Accession) -> ty.Optional[Accession]:
        to_search = id
        # Build an accession that is has all aliases from the NCBI assembly
        # report.
        if info := self.assembly_report.info_for(id):
            to_search = id.alias(info.accession())

        # Check if any requested assembly matches the given or possible
        # aliases.
        for accession in self.requested.standard_accessions:
            if accession.matches(to_search):
                return accession
        return None

    def filter(self, handle: ty.IO) -> ty.Iterable[RecordTypes]:
        """
        Parse a fasta handle and compare each sequence to the requested set. If the
        sequence has been requested, yield a Found object, it has not then yield an
        Extra object. For any ids which are not present in the file yield a Missing
        object. Does not handle duplicate ids within a file.
        """

        count = 0
        seen = SeenAccessions.empty()
        for record in fasta.parse(handle):
            LOGGER.info("Checking if %s is allowed", record.id)
            accession = Accession.build(record.id)
            if self.requested.includes_unplaced() and self.assembly_report.is_unplaced(
                accession
            ):
                LOGGER.info("Found - Accession %s is an expected unplaced", record.id)
                if self.matching_accession(accession):
                    LOGGER.debug("Accession %s was also requested", record.id)
                    seen.mark_accession(accession)
                count += 1
                yield Found(matching_accession=record.id, record=record)

            elif wgs_ids := self.matching_wgs_set(record.id):
                LOGGER.info(
                    "Found - Accession %s is an expected member of WGS set %s",
                    record.id,
                    wgs_ids,
                )
                seen.mark_wgs_prefixes(wgs_ids)
                count += 1
                wgs_id = list(wgs_ids)[0]
                yield Found(matching_accession=wgs_id.to_wgs_string(), record=record)

            elif matching := self.matching_accession(accession):
                LOGGER.info(
                    "Found - Accession %s matches requested %s", record.id, matching
                )
                seen.mark_accession(accession)
                count += 1
                yield Found(matching_accession=record.id, record=record)

            else:
                LOGGER.info("Extra - Accession %s is extra", record.id)
                yield Extra(extra=record)

        LOGGER.info("Found %i requested sequences in the fasta file", count)
        LOGGER.debug("Requested %s", self.requested)
        LOGGER.debug("Seen %s", seen)

        for accession in self.requested.standard_accessions:
            if accession in seen:
                continue
            LOGGER.debug("Missing accession %s", accession)
            yield MissingAccession(accession=str(accession))

        if self.requested.includes_wgs() and not seen.seen_wgs():
            missing = self.requested.wgs_set.any_member_of - seen.wgs_prefix
            LOGGER.debug("Missing wgs sets: %s", missing)
            for wgs_prefix in missing:
                yield MissingWgsSet(prefix=wgs_prefix)
