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

import csv
import enum
import logging
import re
import typing as ty
from functools import lru_cache
from io import StringIO

import attrs
import cattrs
import requests
from attrs import field, frozen
from Bio import SeqIO
from ratelimit import limits, sleep_and_retry
from sqlitedict import SqliteDict

from rfamseq import fasta, wget

LOGGER = logging.getLogger(__name__)

NCBI_SEQ_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"

NCBI_SUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={accessions}"

NCBI_WGS_URL = "https://www.ncbi.nlm.nih.gov/Traces/wgs/{accession}/contigs/tsv"


class UnknownGCF(Exception):
    """
    Raised if an Unknown GCF id is given
    """


class InvalidGenomeId(Exception):
    """
    Raised if given a non GCA/GCF id.
    """


class UnknownGenomeId(Exception):
    """
    Raised if the genome id looks valid but NCBI assembly data does not know
    about it.
    """


@enum.unique
class SequenceRole(enum.Enum):
    ASSEMBLED_MOLECULE = "assembled-molecule"
    UNPLACED_SCAFFOLD = "unplaced-scaffold"
    UNLOCALIZED_SCAFFOLD = "unlocalized-scaffold"
    NOVEL_PATCH = "novel-patch"
    FIX_PATCH = "fix-patch"
    ALTERNATE_SCAFFOLD = "alt-scaffold"


@frozen
class NcbiSequenceInfo:
    genbank_accession: str = field(metadata={"ncbi_name": "GenBank-Accn"})
    name: str = field(metadata={"ncbi_name": "Sequence-Name"})
    role: SequenceRole = field(metadata={"ncbi_name": "Sequence-Role"})
    molecule_type: ty.Optional[str] = field(
        metadata={"ncbi_name": "Assigned-Molecule-Location/Type"}
    )
    length: ty.Optional[int] = field(metadata={"ncbi_name": "Sequence-Length"})


@frozen
class NcbiSequenceSummary:
    accession: str
    title: str
    length: int


@frozen
class NcbiAssemblyInfo:
    taxid: int = field(metadata={"ncbi_name": "Taxid"})
    assembly_name: str = field(metadata={"ncbi_name": "Assembly name"})
    assembly_level: ty.Optional[AssemblyLevel] = field(
        metadata={"ncbi_name": "Assembly level"}
    )
    organism_name: str = field(metadata={"ncbi_name": "Organism name"})
    bio_sample: ty.Optional[str] = field(metadata={"ncbi_name": "BioSample"})
    bio_project: ty.Optional[str] = field(metadata={"ncbi_name": "BioProject"})
    wgs_project: ty.Optional[str] = field(metadata={"ncbi_name": "WGS project"})
    sequence_info: ty.List[NcbiSequenceInfo]

    def info_for(self, accession: str) -> ty.Optional[NcbiSequenceInfo]:
        for info in self.sequence_info:
            if info.genbank_accession == accession:
                return info
        return None


@enum.unique
class RefseqCategory(enum.Enum):
    REFERENCE_GENOME = "reference genome"
    REPRESENTATIVE_GENOME = "representative genome"


@enum.unique
class AssemblyVersionStatus(enum.Enum):
    LATEST = "latest"
    REPLACED = "replaced"
    SUPPRESSED = "suppressed"


@enum.unique
class AssemblyLevel(enum.Enum):
    COMPLETE_GENOME = "Complete Genome"
    CHROMOSOME = "Chromosome"
    SCAFFOLD = "Scaffold"
    CONTIG = "Contig"


@enum.unique
class AssemblyReleaseType(enum.Enum):
    MAJOR = "Major"
    MINOR = "Minor"
    PATCH = "Patch"


@enum.unique
class GenomeAssemblyRep(enum.Enum):
    FULL = "Full"
    PARTIAL = "Partial"


@frozen
class NcbiAssemblySummary:
    assembly_accession: str
    bioproject: str
    biosample: str
    wgs_master: ty.Optional[str]
    refseq_category: ty.Optional[RefseqCategory]
    taxid: int
    species_taxid: int
    organism_name: str
    infraspecific_name: ty.Optional[str]
    isolate: ty.Optional[str]
    version_status: AssemblyVersionStatus
    assembly_level: AssemblyLevel
    release_type: AssemblyReleaseType
    genome_rep: GenomeAssemblyRep
    seq_rel_date: str
    asm_name: str
    submitter: str
    gbrs_paired_asm: ty.Optional[str]
    paired_asm_comp: ty.Optional[str]
    ftp_path: ty.Optional[str]
    excluded_from_refseq: str
    relation_to_type_material: str
    asm_not_live_date: ty.Optional[str]

    @classmethod
    def from_ncbi_row(cls, row: ty.Dict[str, str]) -> NcbiAssemblySummary:
        updated = {k: maybe(v) for k, v in row.items()}
        return cattrs.structure(updated, cls)


def maybe(raw: str) -> ty.Optional[str]:
    if raw == "na" or raw == "":
        return None
    return raw


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

    path = info[versioned].ftp_path
    if path == "na" or path is None:
        LOGGER.info("Accession %s has no path", versioned)
        return None
    parts = path.split("/")
    name = parts[-1]
    return f"{path}/{name}_{suffix}"


def genome_ftp_path(info: SqliteDict, accession: str) -> ty.Optional[str]:
    return ftp_path(info, accession, "genomic.fna.gz")


def assembly_info_path(info: SqliteDict, accession: str) -> ty.Optional[str]:
    return ftp_path(info, accession, "assembly_report.txt")


def parse_header_line(line: str) -> ty.Optional[ty.Tuple[str, ty.Optional[str]]]:
    for field in attrs.fields(NcbiAssemblyInfo):
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
        f.name: None for f in attrs.fields(NcbiAssemblyInfo)
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


def parse_assembly_info(handle: ty.IO) -> ty.Optional[NcbiAssemblyInfo]:
    raw: ty.Dict[str, ty.Any] = parse_header(handle)
    sequences = parse_sequence_lines(handle)
    if not sequences:
        LOGGER.info("No sequences found, returning no assembly info")
        return None

    raw["sequence_info"] = sequences
    return cattrs.structure(raw, NcbiAssemblyInfo)


def assembly_info(info: SqliteDict, accession: str) -> NcbiAssemblyInfo:
    LOGGER.info("Getting NCBI assembly information for %s", accession)
    path = assembly_info_path(info, accession)
    if not path:
        raise UnknownGenomeId(accession)
    with wget.wget(path) as handle:
        ainfo = parse_assembly_info(handle)
        if not ainfo:
            raise UnknownGenomeId(accession)
        return ainfo


@sleep_and_retry
@limits(3, period=1)
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


def resolve_wgs(accession: str) -> ty.Optional[ty.List[str]]:
    LOGGER.info("Trying to resolve WGS set %s", accession)
    url = NCBI_WGS_URL.format(accession=accession)
    LOGGER.debug("Fetching %s", url)
    response = requests.get(url)
    try:
        response.raise_for_status()
    except:
        LOGGER.info("Request to resolve wgs set failed")
        return None

    reader = csv.DictReader(StringIO(response.text), delimiter="\t")
    accessions = [r["accession"] for r in reader]
    if not accessions:
        LOGGER.info("Failed to get load any accessions from response")
        return None
    return accessions
