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

import enum
import json
import logging
import typing as ty
from pathlib import Path

import requests
from attrs import frozen

from rfamseq.converter import camel_case_converter
from rfamseq.uniprot import taxonomy as tax
from rfamseq.utils import batched

LOGGER = logging.getLogger(__file__)


@frozen
class ProteomeTaxonomy:
    taxon_id: str
    scientific_name: None | str = None
    memonic: ty.Optional[str] = None
    common_name: ty.Optional[str] = None


@frozen
class ProteomeGenomeAssembly:
    source: str
    assembly_id: str | None = None
    genome_assembly_url: str | None = None
    level: str | None = None

    def version(self) -> str | None:
        if not self.assembly_id:
            return None
        return self.assembly_id.split(".", 1)[1]


@enum.unique
class ProteomeType(enum.Enum):
    REPRESENTATIVE = "Representative proteome"
    REFERENCE = "Reference proteome"
    REPRESENTATIVE_REFERENCE = "Reference and representative proteome"
    OTHER = "Other proteome"

    def is_reference(self) -> bool:
        return (
            self == ProteomeType.REFERENCE
            or self == ProteomeType.REPRESENTATIVE_REFERENCE
        )

    def is_representative(self) -> bool:
        return (
            self == ProteomeType.REPRESENTATIVE
            or self == ProteomeType.REPRESENTATIVE_REFERENCE
        )


@enum.unique
class ComponentSource(enum.Enum):
    GENOME_ACCESSION = "GenomeAccession"
    BIOSAMPLE = "Biosample"


@frozen
class CrossReference:
    database: ComponentSource
    id: str


@frozen
class GenomeAnnotation:
    source: str


@frozen
class ProteomeComponent:
    name: str
    description: str
    genome_annotation: GenomeAnnotation
    proteome_cross_references: None | ty.List[CrossReference] = None


@frozen
class Proteome:
    id: str
    taxonomy: ProteomeTaxonomy
    modified: str
    proteome_type: ProteomeType
    genome_assembly: ProteomeGenomeAssembly
    components: ty.List[ProteomeComponent]
    taxon_lineage: None | ty.List[tax.LineageEntry] = None
    description: None | str = None


def parse_response(
    data: ty.Dict[str, ty.Any], ignore: None | ty.Set[str] = None
) -> ty.Iterable[Proteome]:
    """Parse a response from the UniProt API containing Proteome data. This will
    generate an iterable of Proteome entries, skipping any in the ignore set.
    It is an error if the response contains no entries or if all entries are
    skipped.
    """

    seen = False
    converter = camel_case_converter()
    if "results" not in data:
        proteome = converter.structure(data, Proteome)
        if ignore and proteome.id in ignore:
            LOGGER.debug("Skipping ignored proteome %s", proteome)
        else:
            seen = True
            yield proteome
    else:
        for entry in data.get("results", []):
            try:
                proteome = converter.structure(entry, Proteome)
            except Exception as err:
                LOGGER.warn("Failed to parse proteome: %s", entry)
                LOGGER.exception(err)
                continue
            if ignore and proteome.id in ignore:
                LOGGER.debug("Skipping ignored proteome %s", proteome)
                continue
            seen = True
            yield proteome

    if not seen:
        raise ValueError(f"Found no proteomes")


def fetch(ids: ty.Set[str], batch_size=200) -> ty.Iterable[Proteome]:
    """Fetch the proteome information for the given proteomes. This will fetch
    proteomes in batch size chunks at a time. It is an error to give no ids to
    fetch. This will also fail if any request to the uniprot API fails. This
    will validate that all requested ids are fetched and there are no duplicates.
    """

    if not ids:
        raise ValueError("Cannot give empty ids to fetch")

    seen = set()
    session = requests.Session()
    for batch in batched(ids, batch_size):
        terms = [f"(upid:{id})" for id in batch]
        query = " OR ".join(terms)
        response = session.get(
            "https://rest.uniprot.org/proteomes/stream",
            params={"format": "json", "query": query},
        )
        response.raise_for_status()
        data = response.json()
        for proteome in parse_response(data):
            if proteome.id in seen:
                raise ValueError(f"Somehow got duplicate proteome: {proteome.id}")
            seen.add(proteome.id)
            yield proteome

    if seen != ids:
        missed = ids - seen
        raise ValueError(f"Did not fetch all requested proteomes {missed}")


def parse(path: Path, ignore: None | ty.Set[str] = None) -> ty.Iterable[Proteome]:
    """Parse a file containing a response from the uniprot API containing Proteomes.
    This will load the entire file into memory, and will yield an iterable of
    proteomes. This will fail if the file contains no proteome entries. Also,
    this will skip all proteomes with a UPID in ignore, it is an error if all
    proteomes are ignored.
    """

    with path.open("r") as raw:
        data = json.load(raw)
        yield from parse_response(data, ignore=ignore)
