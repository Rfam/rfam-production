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

import typing as ty

from attrs import define

from rfamseq import uniprot, wgs
from rfamseq.accession import Accession
from rfamseq.utils import assert_never


@define
class Missing:
    """
    This is meant to track what things which were requested but were missed.
    This is basically a collection of sets for Accessions, WgsPrefixes and
    WgsSequeceIds. This has some logic for remove conceptually duplicate entries,
    ie WgsSequenceIds which are members of a WgsPrefix that is already stored.
    """

    accessions: ty.Set[Accession]
    wgs_sets: ty.Set[wgs.WgsPrefix]
    wgs_sequences: ty.Set[wgs.WgsSequenceId]

    @classmethod
    def empty(cls) -> Missing:
        """
        Create an empty Missing object.
        """
        return Missing(accessions=set(), wgs_sets=set(), wgs_sequences=set())

    def add_components(self, comps: uniprot.SelectedComponents):
        """
        Add all entries in the SelectedComponents to this Missing object. If
        the SelectedComponents includes unplaced items this is an error. This
        cannot track which ids are unplaced since that can only be known by
        parsing an assembly report.
        """

        if comps.unplaced:
            raise ValueError("Cannot add unplaced to Missing")
        for accession in comps.accessions:
            self.add(accession)
        for wgs_set in comps.wgs_sets:
            self.add(wgs_set)
        for wgs_sequence in comps.wgs_sequences:
            self.add(wgs_sequence)

    def add(self, accession: ty.Union[wgs.WgsPrefix, wgs.WgsSequenceId, Accession]):
        """
        Add the given Accession, WgsSequenceId or WgsPrefix to this Missing
        object. If given a WgsPrefix this will remove all WgsSequenceIds which
        have a prefix that is already stored. Additionally, giving a
        WgSequenceId which has a known prefix will not do anything.
        """

        match accession:
            case wgs.WgsSequenceId():
                if not any(w.matches(accession.prefix) for w in self.wgs_sets):
                    self.wgs_sequences.add(accession)
            case wgs.WgsPrefix():
                self.wgs_sets.add(accession)
                to_remove = []
                for wgs_sequence in self.wgs_sequences:
                    if any(w.matches(wgs_sequence.prefix) for w in self.wgs_sets):
                        to_remove.append(wgs_sequence)
                self.wgs_sequences.difference_update(to_remove)
            case Accession():
                self.accessions.add(accession)
            case _:
                assert_never(accession)

    def update(self, other: Missing):
        """
        Update this Missing object with all other entries in the other missing.
        This is the same as updating all involved sets, but does ensure that
        WGS sets/sequences are collapsed properly.
        """

        for accession in other.accessions:
            self.add(accession)
        for wgs_set in other.wgs_sets:
            self.add(wgs_set)
        for wgs_sequence in other.wgs_sequences:
            self.add(wgs_sequence)

    def __bool__(self) -> bool:
        """
        Returns True if there are any accessions, wgs_sets or wgs_sequences.
        """

        return bool(self.accessions) or bool(self.wgs_sets) or bool(self.wgs_sequences)