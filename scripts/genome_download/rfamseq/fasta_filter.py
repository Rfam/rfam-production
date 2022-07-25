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

import re
import typing as ty
import logging
import collections as coll

from attrs import frozen, define
from Bio import SeqIO

from rfamseq import utils, uniprot, fasta, ncbi


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


def wgs_matches(wgs_accession: ty.Dict[str, ty.List[str]], component: str) -> ty.List[str]:
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

    @classmethod
    def from_all(cls) -> ComponentSet:
        return cls(
            components=uniprot.ALL_CHROMOSOMES, allow_unplaced=False, unplaced=set()
        )

    @classmethod
    def from_selected(
        cls,
        sequence_info: ty.List[ncbi.NcbiSequenceInfo],
        selected: uniprot.SelectedComponents,
        wgs_accessions: ty.Optional[ty.Dict[str, ty.List[str]]],
    ) -> ComponentSet:
        components = coll.defaultdict(set)
        allow_unplaced = False
        unplaced = set()

        for component in selected:
            if isinstance(component, uniprot.Unplaced):
                allow_unplaced = True
                continue
            components[component].add(component)
            components[component].add(utils.versionless(component))

            if wgs_accessions and (matches := wgs_matches(wgs_accessions, component)):
                components[component].update(matches)

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
        )

    def normalize(self, accession: str) -> ty.Optional[str]:
        if isinstance(self.components, uniprot.All):
            return None

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

        return None

    def __iter__(self):
        if isinstance(self.components, uniprot.All):
            return iter(())
        return iter(self.components.keys())


def filter(
    records: ty.Iterable[SeqIO.SeqRecord], requested: ComponentSet
) -> ty.Iterable[fasta.RecordTypes]:
    """
    Parse a fasta handle and compare each sequence to the requested set. If the
    sequence has been requestd yeild a Found object, it has not then yeild and
    Extra object. For any ids which are not present in the file yield a Missing
    object. Does not handle duplicate ids within a file.
    """

    seen = set()
    for record in records:
        LOGGER.info("Checking if %s is allowed", record.id)
        normalized_id = requested.normalize(record.id)
        if normalized_id:
            seen.add(normalized_id)
            yield Found(matching_accession=record.id, record=record)
        else:
            yield Extra(extra=record)

    for accession in requested:
        if accession not in seen:
            yield Missing(accession=accession)
