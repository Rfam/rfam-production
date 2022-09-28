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

import datetime as dt
import enum
import typing as ty

import attrs
from attrs import frozen
from Bio import SeqIO
from pypika import Query, Table

from rfamseq import ncbi, uniprot


def as_value(raw) -> ty.Union[int, str, bool]:
    if isinstance(raw, enum.Enum):
        return as_value(raw.value)
    if raw is None or isinstance(raw, (int, str, bool)):
        return raw
    raise ValueError(f"Cannot convert {raw} to sql")


def as_insert(column: str, data: ty.Any) -> Query:
    table = Table(column)
    if not isinstance(data, list):
        data = [data]

    fields = attrs.fields(data[0].__class__)
    columns = [f.metadata.get("column_name", f.name) for f in fields]
    query = Query.into(table).columns(columns)
    for datum in data:
        values = [as_value(getattr(datum, f.name)) for f in fields]
        query = query.insert(tuple(values))
    return query


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
class MoleculeType(enum.Enum):
    GENOMIC_DNA = "genomic DNA"


@frozen
class GenSeq:
    upid: str
    rfamseq_acc: str
    chromosome_name: ty.Optional[str]
    chromosome_type: ty.Optional[str]
    version: str

    @classmethod
    def from_fasta(
        cls,
        upid: str,
        version: str,
        seq_info: ty.Optional[ncbi.NcbiSequenceInfo],
        info: FromFasta,
    ) -> GenSeq:
        chromosome_name = None
        chromosome_type = None
        if seq_info:
            chromosome_name = seq_info.name
            chromosome_type = seq_info.molecule_type

        return cls(
            upid=upid,
            rfamseq_acc=info.rfamseq_acc,
            chromosome_name=chromosome_name,
            chromosome_type=chromosome_type,
            version=version,
        )


@frozen
class Genome:
    upid: str
    assembly_acc: ty.Optional[str]
    assembly_version: ty.Optional[str]
    wgs_acc: ty.Optional[str]
    wgs_version: ty.Optional[int]
    assembly_name: ty.Optional[str]
    assembly_level: ty.Optional[AssemblyLevel]
    study_ref: ty.Optional[str]
    description: ty.Optional[str]
    total_length: ty.Optional[int]
    ungapped_length: ty.Optional[int]
    circular: ty.Optional[bool]
    ncbi_id: int
    scientific_name: str
    common_name: ty.Optional[str]
    kingdom: str
    num_rfam_regions: int
    num_families: int
    is_reference: bool
    is_representative: bool
    created: dt.datetime
    updated: dt.datetime

    @classmethod
    def build(
        cls,
        pinfo: uniprot.ProteomeInfo,
        ncbi_info: ty.Optional[ncbi.NcbiAssemblyReport],
        taxonomy: Taxonomy,
        total_length: int,
    ) -> Genome:
        assembly_name = None
        assembly_level = None
        study_ref = None
        wgs_acc = None
        if ncbi_info:
            wgs_acc = ncbi_info.wgs_project
            study_ref = ncbi_info.bio_project
            assembly_name = ncbi_info.assembly_name
            if ncbi_info.assembly_level:
                assembly_level = AssemblyLevel.from_ncbi(ncbi_info.assembly_level)

        wgs_version = None
        if wgs_acc and "." in wgs_acc:
            wgs_version = int(wgs_acc.split(".", 1)[1])

        now = dt.datetime.now()
        return cls(
            upid=pinfo.upi,
            assembly_acc=pinfo.genome_info.accession,
            wgs_acc=wgs_acc,
            wgs_version=wgs_version,
            assembly_name=assembly_name,
            assembly_level=assembly_level,
            assembly_version=pinfo.genome_info.version,
            study_ref=study_ref,
            description=pinfo.description,
            total_length=total_length,
            ungapped_length=None,
            circular=None,
            ncbi_id=pinfo.taxid,
            scientific_name=pinfo.lineage_info.species,
            common_name=pinfo.lineage_info.common_name,
            kingdom=taxonomy.kingdom(),
            num_rfam_regions=0,
            num_families=0,
            is_reference=pinfo.is_reference,
            is_representative=pinfo.is_representative,
            created=now,
            updated=now,
        )


@frozen
class RfamSeq:
    rfamseq_acc: str
    accession: str
    version: str
    ncbi_id: int
    mol_type: MoleculeType
    length: int
    description: str
    previous_acc: str
    source: str

    @classmethod
    def from_fasta(cls, taxid: int, info: FromFasta) -> RfamSeq:
        accession = info.rfamseq_acc
        version = "1"
        parts = accession.split(".", 1)
        if len(parts) == 2:
            accession, version = parts

        return cls(
            rfamseq_acc=info.rfamseq_acc,
            accession=accession,
            version="%06i" % int(version),
            ncbi_id=taxid,
            mol_type=MoleculeType.GENOMIC_DNA,
            length=info.length,
            description=info.description,
            previous_acc="",
            source="UNIPROT; ENA",
        )


@frozen
class Taxonomy:
    ncbi_id: int
    species: str
    tax_string: str
    tree_display_name: str
    align_display_name: str

    @classmethod
    def from_lineage(cls, info: uniprot.LineageInfo) -> Taxonomy:
        species = info.species
        if info.common_name and "virus" not in info.tax_string:
            species = f"{species} ({info.common_name.lower()})"
        tree_name = species.replace(" ", "_")
        return Taxonomy(
            ncbi_id=info.ncbi_id,
            species=species,
            tax_string=info.tax_string,
            tree_display_name=tree_name,
            align_display_name=f"{tree_name}[{info.ncbi_id}]",
        )

    def kingdom(self) -> str:
        return self.tax_string.split("; ", 1)[0]


@frozen
class Metadata:
    upid: str
    genseq: ty.List[GenSeq]
    rfamseq: ty.List[RfamSeq]
    genome: Genome
    taxonomy: Taxonomy

    @classmethod
    def build(
        cls,
        version: str,
        pinfo: uniprot.ProteomeInfo,
        assembly_info: ty.Optional[ncbi.NcbiAssemblyReport],
        records: ty.List[FromFasta],
    ) -> Metadata:

        genseq = []
        rfamseq = []
        total_length = 0
        for info in records:
            seq_info = None
            if assembly_info:
                seq_info = assembly_info.info_for(info.rfamseq_acc)
            rfamseq.append(RfamSeq.from_fasta(int(pinfo.taxid), info))
            genseq.append(GenSeq.from_fasta(pinfo.upi, version, seq_info, info))
            total_length += info.length

        taxonomy = Taxonomy.from_lineage(pinfo.lineage_info)
        return Metadata(
            upid=pinfo.upi,
            genseq=genseq,
            rfamseq=rfamseq,
            genome=Genome.build(pinfo, assembly_info, taxonomy, total_length),
            taxonomy=taxonomy,
        )

    def as_inserts(self) -> ty.List[Query]:
        return [
            as_insert("rfamseq", self.rfamseq),
            as_insert("genseq", self.genseq),
            as_insert("genome", self.genome),
            as_insert("taxonomy", self.taxonomy),
        ]


@frozen
class FromFasta:
    rfamseq_acc: str
    length: int
    description: str

    @classmethod
    def from_record(cls, record: SeqIO.SeqRecord) -> FromFasta:
        return cls(
            rfamseq_acc=record.id,
            length=len(record.seq),
            description=record.description,
        )
