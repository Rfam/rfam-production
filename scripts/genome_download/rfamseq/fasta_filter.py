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
import re
import typing as ty

from attrs import define, frozen
from Bio import SeqIO

from rfamseq import ncbi, uniprot, utils, wgs
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


def wgs_matches(
    wgs_accession: ty.Dict[str, ty.List[str]], component: str
) -> ty.List[str]:
    for short_wgs, accessions in wgs_accession.items():
        pattern = re.compile(f"{short_wgs}\\d+")
        if re.match(pattern, component):
            return accessions
    return []


@frozen
class RequestedAccessions:
    unplaced: bool
    standard_accessions: ty.Set[Accession]
    wgs_set: ty.Optional[wgs.WgsSummary]

    @classmethod
    def build(cls, selected: uniprot.SelectedComponents) -> RequestedAccessions:
        accessions = set()
        unplaced = False
        wgs_set = None
        for accession in selected.accessions:
            if isinstance(accession, uniprot.Unplaced):
                unplaced = True
                continue
            elif wgs.is_wgs_id(accession):
                # If the WGS accession is the same as what we have in the NCBI
                # report then use the already fetched WGS summary. If it is not
                # then we need to look it up
                pass

            else:
                # Use the NCBI Assembly Summary to build an accession that
                # includes all aliases we can find.
                accessions.add(accession)

        return cls(
            unplaced=unplaced,
            standard_accessions=accessions,
            wgs_set=wgs_set,
        )

    def includes_unplaced(self):
        return self.unplaced

    def __contains__(self, id: Accession) -> bool:
        if any(a.matches(id) for a in self.standard_accessions):
            return True
        return any(id in wgs for wgs in self.wgs_sets)


@define
class SeenAccessions:
    accessions: ty.Set[Accession]
    wgs_id: ty.Set[str]

    @classmethod
    def empty(cls) -> SeenAccessions:
        return cls(accessions=set(), wgs_id=set())

    def mark_wgs(self, wgs_id: str):
        self.wgs_id.add(wgs_id)

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
        accessions: ty.Set[Accession] = set()
        unplaced = False
        for component in selected:
            if isinstance(component, uniprot.Unplaced):
                unplaced = True
            elif wgs.looks_like_wgs_accession(component):
                if wgs_accessions:
                    if wgs_accessions.wgs_id in component:
                        continue
                    # If the component to fetch is a wgs record id which we already have
                    # stored in the wgs_accession object we do not search for it in the
                    # file.
                    # That id will never appear in the file, but the ids which
                    # compose a record may. By this point we have already resolved the
                    # wgs accession into records (hopefully) so we can ignore it as
                    # something to search for and just rely on the wgs ids we have
                    # determined.
                    # If the given component is single version difference from the
                    # record id we already see in the wgs_accessions, we accept it.
                    # This is technically wrong, however, there are cases where old
                    # versions disappear. This is painful but there is generally a new exists,
                    # so we can get a close enough sequence.
                    #
                    # Example:
                    # UP000077684 (GCA_001645045.2) wants LWDE01000000, which is gone,
                    # but LWDE02 exists.
                else:
                    raise ValueError(f"Not yet implemented {component}")
            else:
                accessions.add(Accession.build(component))

        return cls(
            requested=RequestedAccessions(
                unplaced=unplaced,
                standard_accessions=accessions,
                wgs_set=wgs_accessions,
            ),
            assembly_report=assembly_report,
        )

    def matching_wgs_set(self, id: str) -> ty.Optional[str]:
        # Check if this sequence is from a WGS set. This means the id is any one
        # of the ones known to be a part of that WGS set. We will treat seeing any
        # one sequence from a wgs set as seeing the whole WGS set. Slightly
        # risky but should be fine. Maybe.
        #
        # This should take the accession from the sequence and check if it
        # looks like any sequence id know about for the given WGS set
        if self.requested.wgs_set and self.requested.wgs_set.id_matches(
            id, within_one_version=True
        ):
            return self.requested.wgs_set.wgs_id
        return None

    def matching_accession(self, id: Accession) -> ty.Optional[Accession]:
        to_search = id
        if info := self.assembly_report.info_for(id):
            to_search = id.alias(info.accession())
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

        seen = SeenAccessions.empty()
        for record in records:
            LOGGER.info("Checking if %s is allowed", record.id)
            accession = Accession.build(record.id)
            if self.requested.includes_unplaced() and self.assembly_report.is_unplaced(
                accession
            ):
                LOGGER.info("Found - Accession %s is an expected unplaced", record.id)
                yield Found(matching_accession=record.id, record=record)

            elif wgs_id := self.matching_wgs_set(record.id):
                LOGGER.info(
                    "Found - Accession %s is an expected member of WGS set %s",
                    record.id,
                    wgs_id,
                )
                seen.mark_wgs(wgs_id)
                yield Found(matching_accession=wgs_id, record=record)
            elif matching := self.matching_accession(accession):
                LOGGER.info(
                    "Found - Accession %s matches requested %s", record.id, matching
                )
                seen.mark_accession(accession)
                yield Found(matching_accession=record.id, record=record)
            else:
                LOGGER.info("Extra - Accession %s is extra", record.id)
                yield Extra(extra=record)

        import sys

        sys.exit(1)

        for accession in self.requested.standard_accessions:
            if accession in seen:
                continue
            yield Missing(accession=str(accession))

        if self.requested.wgs_set and not seen.seen_wgs():
            for accession in self.requested.wgs_set.largest_ids():
                yield Missing(accession=accession)
