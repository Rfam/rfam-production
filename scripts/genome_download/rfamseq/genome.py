# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

import enum
import typing as ty

from attrs import frozen

from rfamseq import ncbi, wgs
from rfamseq.accession import Accession
from rfamseq.uniprot.proteome import Proteome
from rfamseq.utils import assert_never

Component = ty.Union[Accession, wgs.WgsPrefix, wgs.WgsSequenceId]


@enum.unique
class AssemblyLevel(enum.Enum):
    CONTIG = "contig"
    CHROMOSOME = "chromosome"
    SCAFFOLD = "scaffold"
    COMPLETE_GENOME = "complete-genome"

    @classmethod
    def from_ncbi(cls, level: ncbi.AssemblyLevel) -> AssemblyLevel:
        if hasattr(cls, level.name):
            return getattr(cls, level.name)
        raise ValueError(f"Unknown level {level}")


@enum.unique
class GenomeSource(enum.Enum):
    ENA = "ENA/EMBL"
    REF_SEQ = "REFSEQ"
    ENSEMBL_FUNGI = "ENSEMBLFUNGI"
    ENSEMBL_PROTISTS = "ENSEMBLPROTISTS"
    WORMBASE = "WORMBASE"
    ENSEMBL_PLANTS = "ENSEMBLPLANTS"
    ENSEMBL_METAZOA = "ENSEMBLMETAZOA"
    ENSEMBL = "ENSEMBL"

    def from_ebi(self) -> bool:
        match self:
            case GenomeSource.ENA:
                return True
            case GenomeSource.REF_SEQ:
                return False
            case GenomeSource.ENSEMBL_FUNGI:
                return True
            case GenomeSource.ENSEMBL_PROTISTS:
                return True
            case GenomeSource.WORMBASE:
                return True
            case GenomeSource.ENSEMBL_PLANTS:
                return True
            case GenomeSource.ENSEMBL_METAZOA:
                return True
            case GenomeSource.ENSEMBL:
                return True
            case _:
                assert_never(self)


# @frozen
# class SelectedComponents:
#     unplaced: bool
#     accessions: ty.List[Accession]
#     wgs_sets: ty.List[wgs.WgsPrefix]
#     wgs_sequences: ty.List[wgs.WgsSequenceId]
#
#     @classmethod
#     def build(cls, components: ty.List[Component]) -> SelectedComponents:
#         """
#         Create a new SelectedComponent. This requires the list of components
#         be not empty.
#         """
#         if not components:
#             raise ValueError("Cannot build empty SelectedComponents")
#         unplaced = False
#         accessions = []
#         wgs_sets = list()
#         wgs_sequences = list()
#         for component in components:
#             match component:
#                 case Unplaced():
#                     unplaced = True
#                 case Accession():
#                     accessions.append(component)
#                 case wgs.WgsPrefix():
#                     wgs_sets.append(component)
#                 case wgs.WgsSequenceId():
#                     wgs_sequences.append(component)
#                 case _:
#                     assert_never(component)
#
#         return cls(
#             unplaced=unplaced,
#             accessions=accessions,
#             wgs_sets=wgs_sets,
#             wgs_sequences=wgs_sequences,
#         )
#
#     def matching(self, component: Component) -> ty.List[Component]:
#         match component:
#             case Unplaced():
#                 if self.unplaced:
#                     return [UNPLACED]
#                 return []
#
#             case Accession():
#                 return [a for a in self.accessions if component.matches(a)]
#
#             case wgs.WgsPrefix():
#                 matches: ty.List[Component] = []
#                 for wgs_set in self.wgs_sets:
#                     if wgs_set.matches(component, within_one_version=True):
#                         matches.append(wgs_set)
#                 for wgs_sequence in self.wgs_sequences:
#                     wgs_set = wgs_sequence.prefix
#                     if wgs_set.matches(component, within_one_version=True):
#                         matches.append(wgs_sequence)
#                 return matches
#
#             case wgs.WgsSequenceId():
#                 matches: ty.List[Component] = []
#                 for wgs_set in self.wgs_sets:
#                     if wgs_set.matches(component.prefix, within_one_version=True):
#                         matches.append(wgs_set)
#                 for wgs_sequence in self.wgs_sequences:
#                     if wgs_sequence.matches(component, within_one_version=True):
#                         matches.append(wgs_sequence)
#                 return matches
#
#             case _:
#                 assert_never(component)
#
#     def includes_unplaced(self) -> bool:
#         return self.unplaced
#
#     def includes_wgs(self):
#         return bool(self.wgs_sets) or bool(self.wgs_sequences)


@frozen
class UniprotGenome:
    accession: str
    genome_source: GenomeSource
    proteome: Proteome

    @property
    def upid(self) -> str:
        return self.proteome.id

    @property
    def version(self) -> ty.Optional[str]:
        if "." not in self.accession:
            return None
        return self.accession.split(".", 1)[1]
