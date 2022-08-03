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

import enum
import typing as ty

from attrs import define
from Bio import SeqIO

from rfamseq import ena, fasta_filter, ncbi, uniprot, wget


@enum.unique
class AssemblyLevel(enum.Enum):
    CONTIG = "contig"
    CHROMOSOME = "chromosome"
    SCAFFOLD = "scaffold"
    COMPLETE_GENOME = "complete-genome"


@enum.unique
class MoleculeType(enum.Enum):
    GENOMIC_DNA = "genomic DNA"


@define
class GenSeq:
    rfamseq_acc: str
    chromosome_name: ty.Optional[str]
    chromosome_type: ty.Optional[str]
    version: str


@define
class Genome:
    assembly_acc: ty.Optional[str]
    assembly_version: ty.Optional[str]
    wgs_acc: ty.Optional[str]
    assembly_name: ty.Optional[str]
    assembly_level: ty.Optional[AssemblyLevel]
    study_ref: ty.Optional[str]
    description: str
    total_length: int
    ungapped_length: int
    ciruclar: ty.Optional[bool]
    ncbi_id: int
    scientific_name: str
    common_name: ty.Optional[str]
    kingdom: str

@define
class RfamSeq:
    rfamseq_acc: str
    accession: str
    version: str
    ncbi_id: int
    mol_type: MoleculeType
    length: int
    description: str
    source: str


@define
class Metadata:
    upid: str
    genseq: ty.List[GenSeq]
    rfamseq: ty.List[RfamSeq]


@define
class FromFasta:
    rfamseq_acc: str
    length: int
    accession: str
    version: str
    description: str

    @classmethod
    def from_record(cls, record: SeqIO.SeqRecord) -> FromFasta:
        accession = record.id
        version = "1"
        parts = accession.split(".", 1)
        if len(parts) == 2:
            accession, version = parts
        return cls(
            rfamseq_acc=record.id,
            length=len(record.seq),
            accession=accession,
            version="%05s" % version,
            description=record.description,
        )

def build(version: str, pinfo: uniprot.ProteomeInfo, records: ty.List[FromFasta]) -> Metadata:
    genseq = []
    rfamseq = []
    for info in records:
        genseq.append(GenSeq(
            rfamseq_acc=info.rfamseq_acc,
            chromosome_name=None,
            chromosome_type=None,
            version=version,
        ))

        rfamseq.append(RfamSeq(
            rfamseq_acc=info.rfamseq_acc,
            accession=info.accession,
            version=info.version,
            ncbi_id=int(pinfo.taxid),
            mol_type=MoleculeType.GENOMIC_DNA,
            length=info.length,
            description=info.description,
            source='UNIPROT; ENA',
        ))

    return Metadata(
        upid=pinfo.upi,
        genseq=genseq,
        rfamseq=rfamseq
    )
