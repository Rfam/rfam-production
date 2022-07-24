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

import re
import csv
import logging
import typing as ty
import enum
from xml.etree import ElementTree as ET

import attrs
import cattrs
import requests
from attrs import define, field
from sqlitedict import SqliteDict
from Bio import SeqIO

from rfamseq import wget, fasta

LOGGER = logging.getLogger(__name__)

NCBI_SEQ_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"

NCBI_SUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={accessions}"


class UnknownGCF(Exception):
    """
    Raised if an Unknown GCF id is given
    """

    pass


class InvalidGenomeId(Exception):
    """
    Raised if given a non GCA/GCF id.
    """

    pass


@enum.unique
class SequenceRole(enum.Enum):
    ASSEMBLED_MOLECULE = "assembled-molecule"
    UNPLACED_SCAFFOLD = "unplaced-scaffold"
    UNLOCALIZED_SCAFFOLD = "unlocalized-scaffold"
    NOVEL_PATCH = "novel-patch"
    FIX_PATCH = "fix-patch"
    ALTERNATE_SCAFFOLD = "alt-scaffold"


@define
class NcbiSequenceInfo:
    genbank_accession: str
    name: str
    role: SequenceRole
    molecule_type: ty.Optional[str]
    length: int


@define
class NcbiSequenceSummary:
    accession: str
    title: str
    length: int


@define
class NcbiAssemblyInfo:
    taxid: int = field(metadata={"ncbi_name": "Taxid"})
    assembly_name: str = field(metadata={"ncbi_name": "Assembly name"})
    organism_name: str = field(metadata={"ncbi_name": "Organism name"})
    bio_sample: ty.Optional[str] = field(metadata={"ncbi_name": "BioSample"})
    bio_project: ty.Optional[str] = field(metadata={"ncbi_name": "BioProject"})
    wgs_project: ty.Optional[str] = field(metadata={"ncbi_name": "WGS project"})
    sequence_info: ty.List[NcbiSequenceInfo]


def maybe(raw: str) -> ty.Optional[str]:
    if raw == "na":
        return None
    return raw


def add_version_if_missing(info: SqliteDict, id: str) -> str:
    if "." in id:
        return id
    possible = {}
    pattern = re.compile(f"^{id}.(\\d+)$")
    for key in info.iterkeys():
        if match := re.match(pattern, key):
            index = int(match.group(1))
            possible[index] = key
    to_use = max(possible.keys())
    return possible[to_use]


def ftp_path(info: SqliteDict, accession: str, suffix: str) -> ty.Optional[str]:
    versioned = add_version_if_missing(info, accession)
    if not versioned or versioned not in info:
        LOGGER.info("Accession %s not found in ncbi db", versioned)
        return None

    path = info[accession]["ftp_path"]
    parts = path.split("/")
    name = parts[-1]
    return f"{path}/{name}_{suffix}"


def genome_ftp_path(info: SqliteDict, accession: str) -> ty.Optional[str]:
    return ftp_path(info, accession, "genomic.fna.gz")


def assembly_info_path(info: SqliteDict, accession: str) -> ty.Optional[str]:
    return ftp_path(info, accession, "assembly_report.txt")


def parse_header_line(line: str) -> ty.Optional[ty.Tuple[str, str]]:
    for field in attrs.fields(NcbiAssemblyInfo):
        if field.name == 'sequence_info':
            continue
        prefix = f"# {field.metadata['ncbi_name']}:"
        if not line.startswith(prefix):
            continue
        _, value = line.split(":", 1)
        return (field.name, value.strip())
    return None


def parse_header(handle: ty.IO) -> ty.Dict[str, str]:
    header = {}
    location = handle.tell()
    while line := handle.readline():
        if line.startswith("# Sequence-Name"):
            handle.seek(location)
            return header
        location = handle.tell()
        result = parse_header_line(line)
        if result:
            name, value = result
            if name in header:
                raise ValueError(f"Multiple values for {name} found")
            header[name] = value
    return header


def parse_sequence_lines(handle: ty.IO) -> ty.List[NcbiSequenceInfo]:
    raw_header = handle.readline()
    assert raw_header.startswith("# Sequence-Name"), "Not at header"
    raw_header = raw_header[2:]
    header = raw_header.split("\t")
    reader = csv.DictReader(handle, delimiter="\t", fieldnames=header)
    sequences = []
    for row in reader:
        print(row)
        sequences.append(
            NcbiSequenceInfo(
                genbank_accession=row["GenBank-Accn"],
                name=row["Sequence-Name"],
                role=cattrs.structure(row["Sequence-Role"], SequenceRole),
                molecule_type=maybe(row["Assigned-Molecule-Location/Type"]),
                length=int(row["Sequence-Length"]),
            )
        )
    return sequences


def parse_assembly_info(handle: ty.IO) -> NcbiAssemblyInfo:
    header = parse_header(handle)
    sequences = parse_sequence_lines(handle)
    return NcbiAssemblyInfo(
        taxid=int(header["taxid"]),
        assembly_name=header["assembly_name"],
        organism_name=header["organism_name"],
        bio_sample=header.get("bio_sample", None),
        bio_project=header.get("bio_project", None),
        wgs_project=header.get("wgs_project", None),
        sequence_info=sequences,
    )


def assembly_info(info: SqliteDict, accession: str) -> NcbiAssemblyInfo:
    path = assembly_info_path(info, accession)
    if not path:
        raise ValueError(f"Failed to build path to {accession}")
    with wget.wget(path) as handle:
        return parse_assembly_info(handle)


def parse_doc_sum(entry: ET.Element) -> NcbiSequenceInfo:
    pass


def esummary(accessions: ty.List[str]):
    accs = ",".join(accessions)
    response = requests.get(NCBI_SUMMARY_URL.format(accessions=accs))
    response.raise_for_status()
    root = ET.fromstring(response.content)
    data = {}
    for entry in root:
        if entry.tag == "DocSum":
            accession, info = parse_doc_sum(entry)
            data[accession] = info
        else:
            raise ValueError("Unknown type of element %s" % entry.tag)
    return data


def efetch_fasta(accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    LOGGER.info("Trying efetch for %s", accession)
    url = NCBI_SEQ_URL.format(accession=accession)
    with wget.wget(url) as handle:
        yield from fasta.parse(handle)


def ftp_fasta(info: SqliteDict, accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    LOGGER.info("Trying FTP access to %s", accession)
    prefix = accession[0:4]
    if prefix not in {"GCA_", "GCF_"}:
        raise InvalidGenomeId(accession)
    url = genome_ftp_path(info, accession)
    if url is None and prefix == "GCF_":
        raise UnknownGCF(accession)
    if url is None:
        raise Exception("Not yet implemented")
    with wget.wget(url) as handle:
        yield from fasta.parse(handle)


def fetch_fasta(info: SqliteDict, accession: str) -> ty.Iterable[SeqIO.SeqRecord]:
    if accession.startswith("GCA_") or accession.startswith("GCF_"):
        yield from ftp_fasta(info, accession)
    else:
        yield from efetch_fasta(accession)
