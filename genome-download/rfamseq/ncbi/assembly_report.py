# -*- coding: utf-8 -*-

"""
Copyright [2009-${2022}] EMBL-European Bioinformatics Institute
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

import csv
import enum
import logging
import typing as ty
from functools import lru_cache

import attrs
import cattrs
from attrs import field, frozen
from sqlitedict import SqliteDict

from rfamseq import wget
from rfamseq.accession import Accession
from rfamseq.ncbi import ftp
from rfamseq.ncbi.assembly_summary import AssemblyLevel
from rfamseq.ncbi.utils import maybe

LOGGER = logging.getLogger(__name__)


class UnknownGenomeId(Exception):
    """
    Raised when the genome id cannot be found.
    """


@enum.unique
class SequenceRole(enum.Enum):
    ASSEMBLED_MOLECULE = "assembled-molecule"
    UNPLACED_SCAFFOLD = "unplaced-scaffold"
    UNLOCALIZED_SCAFFOLD = "unlocalized-scaffold"
    NOVEL_PATCH = "novel-patch"
    FIX_PATCH = "fix-patch"
    ALTERNATE_SCAFFOLD = "alt-scaffold"


@enum.unique
class SequenceRelationship(enum.Enum):
    EQUAL = "="
    DIFFERENT = "<>"

    def is_equal(self):
        return self is SequenceRelationship.EQUAL


@frozen
class NcbiSequenceInfo:
    genbank_accession: ty.Optional[Accession] = field(
        metadata={"ncbi_name": "GenBank-Accn"}
    )
    refseq_accession: ty.Optional[Accession] = field(
        metadata={"ncbi_name": "RefSeq-Accn"}
    )
    relationship: SequenceRelationship = field(metadata={"ncbi_name": "Relationship"})
    name: str = field(metadata={"ncbi_name": "Sequence-Name"})
    role: SequenceRole = field(metadata={"ncbi_name": "Sequence-Role"})
    molecule_type: ty.Optional[str] = field(
        metadata={"ncbi_name": "Assigned-Molecule-Location/Type"}
    )
    length: ty.Optional[int] = field(metadata={"ncbi_name": "Sequence-Length"})

    def accession(self, merged=False) -> Accession:
        """
        Get the accession for this sequence. This will prefer the
        genbank_accession over the refseq one if both are set and different. If
        they are equal then an accession with both present will be used. This
        can be forced by setting merged=True.
        """

        primary = None
        if self.genbank_accession:
            primary = self.genbank_accession
        elif self.refseq_accession:
            primary = self.refseq_accession
        else:
            raise ValueError("May not build sequence info without accessions")

        if (
            self.refseq_accession
            and self.refseq_accession != primary
            and (self.relationship.is_equal() or merged)
        ):
            return primary.alias(self.refseq_accession)
        return primary

    def matches(self, accession: Accession) -> bool:
        if self.genbank_accession and self.genbank_accession.matches(accession):
            return True
        if self.relationship.is_equal() and self.refseq_accession:
            return self.refseq_accession.matches(accession)
        return False

    def is_unplaced(self) -> bool:
        return self.role is SequenceRole.UNPLACED_SCAFFOLD


@frozen
class NcbiAssemblyReport:
    taxid: int = field(metadata={"ncbi_name": "Taxid"})
    assembly_name: str = field(metadata={"ncbi_name": "Assembly name"})
    assembly_level: ty.Optional[AssemblyLevel] = field(
        metadata={"ncbi_name": "Assembly level"}
    )
    organism_name: str = field(metadata={"ncbi_name": "Organism name"})
    bio_sample: ty.Optional[str] = field(metadata={"ncbi_name": "BioSample"})
    bio_project: ty.Optional[str] = field(metadata={"ncbi_name": "BioProject"})
    wgs_project: ty.Optional[str] = field(metadata={"ncbi_name": "WGS project"})
    sequence_info: ty.Dict[str, NcbiSequenceInfo]

    def info_for(self, accession: Accession) -> ty.Optional[NcbiSequenceInfo]:
        """
        Find the NcbiSequenceInfo for the given accession.
        """

        versionless = accession.strip_version()
        if info := self.sequence_info.get(str(versionless), None):
            if info.matches(accession):
                return info
        for alias in versionless.aliases:
            if info := self.info_for(alias):
                return info
        return None

    def is_unplaced(self, accession: Accession) -> bool:
        if info := self.info_for(accession):
            return info.is_unplaced()
        return False


@lru_cache
def sequence_header(raw_headers: ty.Tuple[str]) -> ty.List[str]:
    mapping = {}
    for name, field in attrs.fields_dict(NcbiSequenceInfo).items():
        ncbi_name = field.metadata["ncbi_name"]
        mapping[ncbi_name] = name

    headers = []
    for raw in raw_headers:
        if raw in mapping:
            headers.append(mapping[raw])
        else:
            headers.append(raw)
    return headers


def parse_header_line(line: str) -> ty.Optional[ty.Tuple[str, ty.Optional[str]]]:
    for field in attrs.fields(NcbiAssemblyReport):
        if field.name == "sequence_info":
            continue
        prefix = f"# {field.metadata['ncbi_name']}:"
        if not line.startswith(prefix):
            continue
        _, value = line.split(":", 1)
        return (field.name, maybe(value.strip()))
    return None


def parse_header(handle: ty.IO) -> ty.Dict[str, ty.Any]:
    header: ty.Dict[str, ty.Any] = {
        f.name: None for f in attrs.fields(NcbiAssemblyReport)
    }

    seen = set()
    location = handle.tell()
    while line := handle.readline():
        if line.startswith("# Sequence-Name"):
            handle.seek(location)
            return header
        location = handle.tell()
        result = parse_header_line(line)
        if result:
            name, value = result
            if name in seen:
                raise ValueError(f"Multiple values for {name} found")
            header[name] = value
            seen.add(name)

    return header


def parse_sequence_lines(handle: ty.IO) -> ty.List[ty.Dict[str, ty.Optional[str]]]:
    raw_header = handle.readline()
    if not raw_header.startswith("# Sequence-Name"):
        LOGGER.info("Assembly file appears to be empty")
        return []
    header = sequence_header(tuple(raw_header[2:].strip().split("\t")))
    reader = csv.DictReader(handle, delimiter="\t", fieldnames=header)
    sequences = []
    for row in reader:
        converted = {}
        for key, value in row.items():
            converted[key] = maybe(value)
        sequences.append(converted)
    return sequences


def parse_assembly_info(handle: ty.IO) -> ty.Optional[NcbiAssemblyReport]:
    raw: ty.Dict[str, ty.Any] = parse_header(handle)
    sequences = parse_sequence_lines(handle)
    if not sequences:
        LOGGER.info("No sequences found, returning no assembly info")
        return None

    raw["sequence_info"] = {}
    for sequence in sequences:
        info = cattrs.structure(sequence, NcbiSequenceInfo)
        if info.genbank_accession:
            acc = str(info.genbank_accession.strip_version())
            raw["sequence_info"][acc] = sequence
        if info.refseq_accession:
            acc = str(info.refseq_accession.strip_version())
            raw["sequence_info"][acc] = sequence
    assert len(raw["sequence_info"]) >= len(sequences)
    return cattrs.structure(raw, NcbiAssemblyReport)


def ftp_path(info: SqliteDict, accession: str) -> ty.Optional[str]:
    return ftp.ftp_path(info, accession, "assembly_report.txt")


def fetch_assembly_report(info: SqliteDict, accession: str) -> NcbiAssemblyReport:
    LOGGER.info("Getting NCBI assembly information for %s", accession)
    path = ftp_path(info, accession)
    if not path:
        raise UnknownGenomeId(accession)
    LOGGER.info("Using path %s", path)

    try:
        with wget.wget(path) as handle:
            ainfo = parse_assembly_info(handle)
            if not ainfo:
                raise UnknownGenomeId(accession)
            return ainfo
    except Exception as err:
        LOGGER.debug(err)

    raise ValueError("Failed")
