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

import cattrs
from attrs import frozen

from rfamseq.ncbi.utils import maybe

LOGGER = logging.getLogger(__name__)


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
    asm_submitter: str
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


def cleaned_assembly(handle: ty.IO) -> ty.Iterable[str]:
    for index, line in enumerate(handle):
        if index == 0:
            continue
        if index == 1:
            yield re.sub(r"#\s*", "", line)
        else:
            yield line


def parse_assembly_files(filenames: ty.List[str]) -> ty.Iterable[NcbiAssemblySummary]:
    """
    Iterates over the list of filenames and parses all assemblies from each
    one. Each file must contain at least one valid assembly summary otherwise
    it will fail. This will log any inputs that fail to be parsed but will not
    fail in that case.
    """

    for filename in filenames:
        LOGGER.debug("Parsing assemblies from %s", filename)
        count = 0
        with open(filename, "r") as handle:
            lines = cleaned_assembly(handle)
            for row in csv.DictReader(lines, delimiter="\t"):
                try:
                    yield NcbiAssemblySummary.from_ncbi_row(row)
                    count += 1
                except Exception:
                    LOGGER.error("Could not parse row %s", row)
                    continue

        if not count:
            raise ValueError(f"Failed to parse any assemblies from {filename}")
