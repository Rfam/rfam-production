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

from attrs import frozen
from Bio import SeqIO

from rfamseq import ncbi, uniprot, utils, wgs

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
class ComponentSet:
    components: ty.Union[uniprot.All, ty.Dict[str, ty.Set[str]]]
    allow_unplaced: bool
    unplaced: ty.Set[str]
    wgs_sets: ty.Optional[wgs.WgsSummary]

    @classmethod
    def from_all(cls) -> ComponentSet:
        return cls(
            components=uniprot.ALL_CHROMOSOMES,
            allow_unplaced=False,
            unplaced=set(),
            wgs_sets=None,
        )

    @classmethod
    def from_selected(
        cls,
        sequence_info: ty.List[ncbi.NcbiSequenceInfo],
        selected: uniprot.SelectedComponents,
        wgs_accessions: ty.Optional[wgs.WgsSummary],
    ) -> ComponentSet:
        components = coll.defaultdict(set)
        allow_unplaced = False
        unplaced = set()

        for component in selected:
            if isinstance(component, uniprot.Unplaced):
                allow_unplaced = True
                continue

            # If the component to fetch is a wgs record id which we already have
            # stored in the wgs_accession object we do not search for it in the
            # file. That id will never appear in the file, but the ids which
            # compose a record may. By this point we have already resolved the
            # wgs accession into records (hopefully) so we can ignore it as
            # something to search for and just rely on the wgs ids we have
            # determined.
            if wgs_accessions and component in wgs_accessions.record_ids():
                continue

            # If the given component is single version difference from the
            # record id we already see in the wgs_accessions, we accept it.
            # This is techinically wrong, however, there are cases where old
            # versions dissappear. This is painful but there is generally a new exists,
            # so we can get a close enough sequence.
            #
            # Example:
            # UP000077684 (GCA_001645045.2) wants LWDE01000000, which is gone,
            # but LWDE02 exists.
            if wgs_accessions and wgs_accessions.within_one_version(component):
                continue

            components[component].add(component)
            components[component].add(utils.versionless(component))

        if allow_unplaced:
            for info in sequence_info:
                if info.role is not ncbi.SequenceRole.UNPLACED_SCAFFOLD:
                    continue
                unplaced.add(info.genbank_accession)
                unplaced.add(utils.versionless(info.genbank_accession))

        return cls(
            components=dict(components),
            allow_unplaced=allow_unplaced,
            unplaced=unplaced,
            wgs_sets=wgs_accessions,
        )

    def wgs_project(self, accession: str) -> ty.Optional[str]:
        if isinstance(self.components, uniprot.All):
            return None
        if not self.wgs_sets:
            return None

        if accession in self.wgs_sets:
            return self.wgs_sets.wgs_id

        return None

    def normalize(self, accession: str) -> ty.Optional[str]:
        if isinstance(self.components, uniprot.All):
            return accession

        versionless = utils.versionless(accession)
        for key, ids in self.components.items():
            if (
                key == accession
                or accession in ids
                or key == versionless
                or versionless in ids
            ):
                return key

        if self.allow_unplaced:
            if accession in self.unplaced:
                return accession

        if self.wgs_sets:
            if versionless in self.wgs_sets:
                return versionless
            if accession in self.wgs_sets:
                return accession

        return None

    def __iter__(self):
        if isinstance(self.components, uniprot.All):
            return iter(())

        yield from self.components.keys()
        if self.wgs_sets:
            yield from self.wgs_sets.ids()


def filter(
    records: ty.Iterable[SeqIO.SeqRecord], requested: ComponentSet
) -> ty.Iterable[RecordTypes]:
    """
    Parse a fasta handle and compare each sequence to the requested set. If the
    sequence has been requestd yeild a Found object, it has not then yeild and
    Extra object. For any ids which are not present in the file yield a Missing
    object. Does not handle duplicate ids within a file.
    """

    seen = set()
    seen_wgs = set()
    for record in records:
        LOGGER.info("Checking if %s is allowed", record.id)
        normalized_id = requested.normalize(record.id)
        if normalized_id:
            seen.add(normalized_id)
            if wgs_id := requested.wgs_project(normalized_id):
                seen_wgs.add(wgs_id)
            yield Found(matching_accession=record.id, record=record)
        else:
            yield Extra(extra=record)

        # import sys
        # sys.exit()

    for accession in requested:
        if accession in seen:
            continue
        wgs_project = requested.wgs_project(accession)
        if wgs_project in seen_wgs:
            continue
        yield Missing(accession=accession)
