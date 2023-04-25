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
    matching_accession: str
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

    def seen_accession(
        self, accession: ty.Union[Accession, wgs.WgsSequenceId, wgs.WgsPrefix]
    ) -> bool:
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

    def mark(self, seq_id: uniprot.All | uniprot.Component):
        self.count += 1
        match seq_id:
            case uniprot.All():
                pass  # Already handled above
            case uniprot.Unplaced():
                self.unplaced = True
            case Accession():
                self.accessions.add(seq_id)
            case wgs.WgsPrefix():
                self.wgs_prefix.add(seq_id)
            case wgs.WgsSequenceId():
                self.wgs_sequences.add(seq_id)
                self.wgs_prefix.add(seq_id.prefix)
            case _:
                assert_never(seq_id)


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

    def requested_component(
        self, accession: ty.Union[Accession, wgs.WgsSequenceId]
    ) -> ty.List[ty.Union[uniprot.All, uniprot.Component]]:
        match self.requested:
            case uniprot.All():
                return [self.requested]

            case uniprot.SelectedComponents():
                matches = []
                to_lookup = None
                match accession:
                    case Accession():
                        # Build an accession that is has all aliases from the NCBI assembly
                        # report.
                        to_search = accession
                        if self.assembly_report and (
                            info := self.assembly_report.info_for(accession)
                        ):
                            to_search = accession.alias(info.accession())
                        matches.extend(self.requested.matching_accessions(to_search))
                        to_lookup = accession

                    case wgs.WgsSequenceId():
                        matches.extend(
                            self.requested.matching_wgs_sets(accession.prefix)
                        )
                        matches.extend(self.requested.matching_wgs_sequences(accession))
                        to_lookup = Accession.build(accession.to_wgs_string())

                    case _:
                        assert_never(accession)

                if (
                    self.requested.unplaced
                    and self.assembly_report
                    and self.assembly_report.is_unplaced(to_lookup)
                ):
                    matches.append(uniprot.UNPLACED)

                return matches

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
            accession = wgs.parse_wgs_sequence_id(record.id)
            if not accession:
                accession = Accession.build(record.id)

            if matching := self.requested_component(accession):
                LOGGER.info("Found - Accession %s was requested", accession)
                for found in matching:
                    LOGGER.debug("Accession %s matches %s", accession, found)
                    seen.mark(found)

                if not len(record.seq):
                    raise ValueError(f"Record {record.id} is empty")

                yield Found(matching_accession=record.id, record=record)

            else:
                LOGGER.info("Extra - Accession %s is extra", record.id)
                pass

        LOGGER.debug("Requested %s", self.requested)
        LOGGER.debug("Seen %s", seen)

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
                            accession = wgs.parse_wgs_sequence_id(
                                str(sequence_info.accession())
                            )
                            if not accession:
                                accession = sequence_info.accession()
                            if seen.seen_accession(accession):
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
                    if not seen.seen_accession(accession):
                        yield MissingAccession(accession=accession)

                for wgs_sequence in self.requested.wgs_sequences:
                    if not seen.seen_accession(wgs_sequence):
                        yield MissingWgsSequence(sequence_id=wgs_sequence)

                for wgs_prefix in self.requested.wgs_sets:
                    if not seen.seen_accession(wgs_prefix):
                        yield MissingWgsSet(prefix=wgs_prefix)

            case _:
                assert_never(self.requested)
