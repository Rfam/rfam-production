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

import logging
import typing as ty

from attrs import define, frozen
from Bio import SeqIO

from rfamseq import fasta, ncbi, uniprot, wgs
from rfamseq.accession import Accession
from rfamseq.utils import assert_never

LOGGER = logging.getLogger(__name__)

SequenceId = ty.Union[Accession, wgs.WgsSequenceId]


def parse_sequence_id(raw: str) -> SequenceId:
    if accession := wgs.parse_wgs_sequence_id(raw):
        return accession
    return Accession.build(raw)


def sequence_id_str(seq_id: SequenceId) -> str:
    match seq_id:
        case Accession():
            return str(seq_id)
        case wgs.WgsSequenceId():
            return seq_id.to_wgs_string()
        case _:
            assert_never(seq_id)


@frozen
class MissingAccession:
    accession: Accession


@frozen
class MissingWgsSet:
    prefix: wgs.WgsPrefix


@frozen
class MissingWgsSequence:
    sequence_id: wgs.WgsSequenceId


@frozen
class Found:
    record: SeqIO.SeqRecord


Missing = ty.Union[MissingAccession, MissingWgsSet, MissingWgsSequence]
RecordTypes = ty.Union[Missing, Found]


@define
class AccessionTracker:
    count: int
    unplaced: bool
    accessions: ty.Set[Accession]
    wgs_sequences: ty.Set[wgs.WgsSequenceId]
    wgs_prefix: ty.Set[wgs.WgsPrefix]

    @classmethod
    def empty(cls) -> AccessionTracker:
        return cls(
            count=0,
            unplaced=False,
            accessions=set(),
            wgs_prefix=set(),
            wgs_sequences=set(),
        )

    def seen(self, accession: ty.Union[SequenceId, wgs.WgsPrefix]) -> bool:
        match accession:
            case Accession():
                return any(a.matches(accession) for a in self.accessions)
            case wgs.WgsSequenceId():
                return any(
                    w.matches(accession, within_one_version=True)
                    for w in self.wgs_sequences
                ) or any(
                    w.matches(accession.prefix, within_one_version=True)
                    for w in self.wgs_prefix
                )
            case wgs.WgsPrefix():
                return any(
                    w.matches(accession, within_one_version=True)
                    for w in self.wgs_prefix
                )
            case _:
                assert_never(accession)

    def mark(self, seq_id: uniprot.All | ty.List[uniprot.Component]):
        self.count += 1
        if isinstance(seq_id, uniprot.All):
            return  # Taken care of above with += 1

        for sid in seq_id:
            match sid:
                case uniprot.Unplaced():
                    self.unplaced = True
                case Accession():
                    self.accessions.add(sid)
                case wgs.WgsPrefix():
                    self.wgs_prefix.add(sid)
                case wgs.WgsSequenceId():
                    self.wgs_sequences.add(sid)
                    self.wgs_prefix.add(sid.prefix)
                case _:
                    assert_never(sid)


@frozen
class FastaFilter:
    requested: uniprot.Components
    assembly_report: ty.Optional[ncbi.NcbiAssemblyReport]

    @classmethod
    def from_selected(
        cls,
        assembly_report: ty.Optional[ncbi.NcbiAssemblyReport],
        requested: uniprot.Components,
    ) -> FastaFilter:
        return cls(requested=requested, assembly_report=assembly_report)

    def matching_components(
        self, accession: SequenceId
    ) -> uniprot.All | ty.List[uniprot.Component]:
        match self.requested:
            case uniprot.All():
                return self.requested

            case uniprot.SelectedComponents():
                matching = set(self.requested.matching(accession))
                if self.assembly_report:
                    to_lookup = None
                    if isinstance(accession, wgs.WgsSequenceId):
                        text = accession.to_wgs_string()
                        to_lookup = Accession.build(text)
                    else:
                        to_lookup = accession

                    if info := self.assembly_report.info_for(to_lookup):
                        aliased = info.accession()
                        matching.update(self.requested.matching(aliased))
                        if info.is_unplaced():
                            matching.add(uniprot.UNPLACED)

                return list(matching)

            case _:
                assert_never(self.requested)

    def __missing__(self, seen: AccessionTracker) -> ty.Iterable[RecordTypes]:
        match self.requested:
            case uniprot.All():
                if not seen.count:
                    raise ValueError("Asked for all, saw nothing")

            case uniprot.SelectedComponents():
                if self.requested.unplaced and not seen.unplaced:
                    if not self.assembly_report:
                        LOGGER.error("Asked for unplaced but cannot infer any")
                    else:
                        for sequence_info in self.assembly_report.sequence_info:
                            if not sequence_info.is_unplaced():
                                continue

                            accession = parse_sequence_id(
                                str(sequence_info.accession())
                            )
                            if seen.seen(accession):
                                LOGGER.debug(
                                    "There appears to be confusion about %s", accession
                                )
                                continue
                            match accession:
                                case wgs.WgsSequenceId():
                                    yield MissingWgsSequence(sequence_id=accession)
                                case Accession():
                                    yield MissingAccession(accession=accession)
                                case _:
                                    assert_never(accession)

                for accession in self.requested.accessions:
                    if not seen.seen(accession):
                        yield MissingAccession(accession=accession)

                for wgs_sequence in self.requested.wgs_sequences:
                    if not seen.seen(wgs_sequence):
                        yield MissingWgsSequence(sequence_id=wgs_sequence)

                for wgs_prefix in self.requested.wgs_sets:
                    if not seen.seen(wgs_prefix):
                        yield MissingWgsSet(prefix=wgs_prefix)

            case _:
                assert_never(self.requested)

    def filter(self, handle: ty.IO) -> ty.Iterable[RecordTypes]:
        """
        Parse a fasta handle and compare each sequence to the requested set. If the
        sequence has been requested, yield a Found object, it has not then yield an
        Extra object. For any ids which are not present in the file yield a Missing
        object. Does not handle duplicate ids within a file.
        """

        seen = AccessionTracker.empty()
        for record in fasta.parse(handle):
            LOGGER.info("Checking if %s is allowed", record.id)
            accession = parse_sequence_id(record.id)
            if matching := self.matching_components(accession):
                LOGGER.info("Found - Accession %s was requested", accession)
                LOGGER.debug("Accession %s matches %s", accession, matching)
                seen.mark(matching)

                if not len(record.seq):
                    raise ValueError(f"Record {record.id} is empty")

                yield Found(record=record)

            else:
                LOGGER.info("Extra - Accession %s is extra", record.id)

        LOGGER.debug("Requested %s", self.requested)
        LOGGER.debug("Seen %s", seen)

        yield from self.__missing__(seen)
