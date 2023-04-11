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

from rfamseq import ncbi, uniprot, wgs
from rfamseq.accession import Accession

LOGGER = logging.getLogger(__name__)


@frozen
class Missing:
    accession: str


@frozen
class Found:
    matching_accession: str
    record: SeqIO.SeqRecord


@frozen
class Extra:
    extra: SeqIO.SeqRecord


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

    def add(self, accession: str):
        try:
            prefix = wgs.WgsPrefix.build(accession)
            self.any_member_of.add(prefix)
            return None
        except:
            LOGGER.debug("Accession %s is not a wgs prefix", accession)
            pass

        try:
            range = wgs.WgsSequenceRange.build(wgs.WgsSequenceKind.SEQUENCE, accession)
            return self.__add_sequence_range(range)
        except:
            LOGGER.debug("Accession %s is not a wgs sequence range", accession)
            pass

        raise ValueError(f"Cannot treat {accession} as a WGS accession")

    def is_empty(self) -> bool:
        return not (bool(self.ranges) or bool(self.any_member_of))

    def matching_set(
        self,
        accession: ty.Union[str, Accession, wgs.WgsSequenceId],
        within_one_version=False,
    ) -> ty.Optional[wgs.WgsPrefix]:
        if isinstance(accession, Accession):
            accession = wgs.WgsSequenceId.build(str(accession))
        if isinstance(accession, str):
            accession = wgs.WgsSequenceId.build(accession)

        for prefix in self.any_member_of:
            if prefix.matches(accession.prefix, within_one_version=within_one_version):
                return prefix

        for range in self.ranges:
            if range.includes(accession, within_one_version=within_one_version):
                return range.prefix
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
            elif wgs.looks_like_wgs_accession(accession):
                wgs_set.add(accession)
            else:
                accessions.add(Accession.build(accession))

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
    wgs_id: coll.Counter[str]

    @classmethod
    def empty(cls) -> SeenAccessions:
        return cls(accessions=set(), wgs_id=coll.Counter())

    def mark_wgs(self, wgs_id: str):
        self.wgs_id[wgs_id] += 1

    def mark_accession(self, accession: Accession):
        self.accessions.add(accession)

    def seen_wgs(self) -> bool:
        return bool(self.wgs_id)

    def __contains__(self, accession: Accession) -> bool:
        return any(a.matches(accession) for a in self.accessions)


@frozen
class ComponentSelector:
    requested: RequestedAccessions
    assembly_report: ncbi.NcbiAssemblyReport

    @classmethod
    def from_selected(
        cls,
        assembly_report: ncbi.NcbiAssemblyReport,
        selected: uniprot.SelectedComponents,
        wgs_accessions: ty.Optional[wgs.WgsSummary],
    ) -> ComponentSelector:
        return cls(
            requested=RequestedAccessions.build(selected),
            assembly_report=assembly_report,
        )

    def matching_wgs_set(self, id: str) -> ty.Optional[str]:
        if not wgs.looks_like_wgs_accession(id):
            return None

        # Check if this sequence is from a WGS set. This means the id is any one
        # of the ones known to be a part of that WGS set. We will treat seeing any
        # one sequence from a wgs set as seeing the whole WGS set. Slightly
        # risky but should be fine. Maybe.
        #
        # This should take the accession from the sequence and check if it
        # looks like any sequence id know about for the given WGS set
        if prefix := self.requested.wgs_set.matching_set(id, within_one_version=True):
            return prefix.to_wgs_string()
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

    def filter(self, records: ty.Iterable[SeqIO.SeqRecord]) -> ty.Iterable[RecordTypes]:
        """
        Parse a fasta handle and compare each sequence to the requested set. If the
        sequence has been requested, yield a Found object, it has not then yield an
        Extra object. For any ids which are not present in the file yield a Missing
        object. Does not handle duplicate ids within a file.
        """

        count = 0
        seen = SeenAccessions.empty()
        for record in records:
            LOGGER.info("Checking if %s is allowed", record.id)
            accession = Accession.build(record.id)
            if self.requested.includes_unplaced() and self.assembly_report.is_unplaced(
                accession
            ):
                LOGGER.info("Found - Accession %s is an expected unplaced", record.id)
                count += 1
                yield Found(matching_accession=record.id, record=record)

            elif wgs_id := self.matching_wgs_set(record.id):
                LOGGER.info(
                    "Found - Accession %s is an expected member of WGS set %s",
                    record.id,
                    wgs_id,
                )
                seen.mark_wgs(wgs_id)
                count += 1
                yield Found(matching_accession=wgs_id, record=record)

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
            yield Missing(accession=str(accession))

        if self.requested.includes_wgs() and not seen.seen_wgs():
            # for accession in self.requested.wgs_set.largest_ids():
            #     LOGGER.debug("Missing sequences from WGS set %s", accession)
            raise ValueError("Not reimplemented")
