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
import logging
import typing as ty
import xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path

import cattrs
import requests
from attrs import field, frozen

from rfamseq import wgs
from rfamseq.accession import Accession
from rfamseq.utils import assert_never

LOGGER = logging.getLogger(__name__)

NS = {
    "uni": "http://uniprot.org/uniprot",
    "pro": "http://uniprot.org/proteome",
}

T = ty.TypeVar("T")


@frozen
class All:
    all: bool


@frozen
class Unplaced:
    unplaced: bool


UNPLACED = Unplaced(unplaced=True)

ALL_CHROMOSOMES = All(all=True)

Component = ty.Union[Accession, wgs.WgsPrefix, wgs.WgsSequenceId, Unplaced]


def structure_component(v: ty.Any, _) -> Component:
    if isinstance(v, dict):
        if "unplaced" in v:
            return Unplaced(**v)
        if "accession" in v:
            if "aliases" in v:
                v["aliases"] = tuple(v["aliases"])
            return Accession(**v)
        if "wgs_id" in v:
            return wgs.WgsPrefix(**v)
        if "prefix" in v:
            return wgs.WgsSequenceId(**v)
    raise ValueError(f"Cannot unstructure {v} to Component")


cattrs.register_structure_hook(Component, structure_component)


@frozen
class SelectedComponents:
    unplaced: bool
    accessions: ty.List[Accession]
    wgs_sets: ty.List[wgs.WgsPrefix]
    wgs_sequences: ty.List[wgs.WgsSequenceId]

    @classmethod
    def build(cls, components: ty.List[Component]) -> SelectedComponents:
        unplaced = False
        accessions = []
        wgs_sets = list()
        wgs_sequences = list()
        for component in components:
            match component:
                case Unplaced():
                    unplaced = True
                case Accession():
                    accessions.append(component)
                case wgs.WgsPrefix():
                    wgs_sets.append(component)
                case wgs.WgsSequenceId():
                    wgs_sequences.append(component)
                case _:
                    assert_never(component)

        return cls(
            unplaced=unplaced,
            accessions=accessions,
            wgs_sets=wgs_sets,
            wgs_sequences=wgs_sequences,
        )

    def matching_accessions(self, accession: Accession) -> ty.List[Accession]:
        return [a for a in self.accessions if a.matches(accession)]

    def matching_wgs_sets(self, prefix: wgs.WgsPrefix) -> ty.List[wgs.WgsPrefix]:
        return [a for a in self.wgs_sets if a.matches(prefix, within_one_version=True)]

    def matching_wgs_sequences(
        self, seq_id: wgs.WgsSequenceId
    ) -> ty.List[wgs.WgsSequenceId]:
        return [a for a in self.wgs_sequences if a.matches(seq_id)]

    def includes_unplaced(self) -> bool:
        return self.unplaced

    def includes_wgs(self):
        return bool(self.wgs_sets) or bool(self.wgs_sequences)


Components = ty.Union[All, SelectedComponents]


@enum.unique
class GenomeSource(enum.Enum):
    ENA = "ENA/EMBL"
    REF_SEQ = "Refseq"
    ENSEMBL_FUNGI = "EnsemblFungi"
    ENSEMBL_PROTISTS = "EnsemblProtists"
    WORMBASE = "WormBase"
    ENSEMBL_PLANTS = "EnsemblPlants"
    ENSEMBL_METAZOA = "EnsemblMetazoa"
    ENSEMBL = "Ensembl"


@frozen
class GenomeInfo:
    accession: ty.Optional[str]
    description: ty.Optional[str]
    components: Components
    source: ty.Optional[GenomeSource]

    @property
    def version(self) -> ty.Optional[str]:
        if not self.accession or "." not in self.accession:
            return None
        return self.accession.split(".", 1)[1]


@frozen
class LineageInfo:
    ncbi_id: int
    species: str
    common_name: ty.Optional[str]
    tax_string: str


@frozen
class ProteomeInfo:
    upi: str
    taxid: int
    is_reference: bool
    is_representative: bool
    proteome_description: ty.Optional[str]
    genome_info: GenomeInfo
    lineage_info: LineageInfo

    @property
    def description(self) -> ty.Optional[str]:
        return self.proteome_description or self.genome_info.description


def node_text(node: ty.Optional[ET.Element]) -> ty.Optional[str]:
    if node is not None and node.text is not None:
        return node.text
    return None


def all_matching_node_text(query: str, root: ET.Element) -> ty.Optional[ty.List[str]]:
    nodes = root.findall(query, namespaces=NS)
    if len(nodes) == 0:
        return None
    found = []
    for node in nodes:
        accession = node_text(node)
        if not accession:
            return None
        found.append(accession)
    return found


def as_bool(raw: str) -> bool:
    if raw == "true":
        return True
    if raw == "false":
        return False
    raise ValueError("Cannot convert %s to a bool" % raw)


def proteome_value(
    xml: ET.Element, tag: str, convert: ty.Callable[[str], T] = str
) -> T:
    value = xml.find(tag, NS)
    if value is None:
        raise ValueError(f"Invalid xml, no {tag}")
    value = value.text
    if value is None:
        raise ValueError(f"Invalid xml, empty {tag}")
    return convert(value)


def genome_info(upid: str, root: ET.Element) -> GenomeInfo:
    accessions = all_matching_node_text(
        ".//pro:genomeAssembly/pro:genomeAssembly", root
    )
    if accessions and len(accessions) > 1:
        raise ValueError(f"Cannot handle >1 genome accession in {upid}")

    source = None
    sources = all_matching_node_text(
        ".//pro:genomeAssembly/pro:genomeAssemblySource", root
    )
    if sources:
        LOGGER.debug("Found genome sources: %s", sources)
        if len(sources) == 1:
            source = cattrs.structure(sources[0], GenomeSource)

    components: ty.List[Component] = []
    description = None
    saw_genome = False
    for component in root.findall("pro:component", NS):
        name = component.attrib.get("name", "")
        if "unplaced" in name.lower():
            components.append(UNPLACED)
            continue

        comp_accs = all_matching_node_text("pro:genomeAccession", component)
        if not comp_accs:
            raise ValueError(f"Missing component accession of {upid}")
        if name == "Genome":
            saw_genome = True
            descriptions = all_matching_node_text("pro:description", component)
            if descriptions and len(descriptions) == 1:
                description = descriptions[0]

        if len(comp_accs) > 1 and name != "Genome":
            raise ValueError(f"Cannot handle >1 accession unless it is a genome {upid}")

        for raw_acc in comp_accs:
            if wgs_component := wgs.parse_wgs_accession(raw_acc):
                components.append(wgs_component)
            else:
                components.append(Accession.build(raw_acc))

    # Handle the case where there is not overarching genome assembly to work
    # with and just a series of components to fetch.
    if not accessions:
        # If the only component is a GCA_ sequence then use that and all
        # components.
        if len(components) == 1:
            possible_gca = components[0]
            if isinstance(possible_gca, Accession):
                possible_gca = str(possible_gca)
            if isinstance(possible_gca, str) and possible_gca.startswith("GCA_"):
                return GenomeInfo(
                    accession=possible_gca,
                    description=description,
                    components=ALL_CHROMOSOMES,
                    source=source,
                )

            # If there is one component marked as a genome use that and all
            # chromosomes in it.
            if saw_genome:
                if not isinstance(possible_gca, str):
                    raise ValueError(f"Invalid state for seeing genome in {upid}")
                return GenomeInfo(
                    accession=possible_gca,
                    description=description,
                    components=ALL_CHROMOSOMES,
                    source=source,
                )

        # Otherwise only use the given components
        return GenomeInfo(
            accession=None,
            description=description,
            components=SelectedComponents.build(components),
            source=source,
        )

    assert len(accessions) == 1, f"Multiple genomes, impossibility in {upid}"

    # If not given any components then assume we want all sequences
    if not components:
        return GenomeInfo(
            accession=accessions[0],
            description=description,
            components=ALL_CHROMOSOMES,
            source=source,
        )

    # Handle being asked for a single component that is the versionless version
    # of the genome accession. We assume that we want all sequences in the
    # genome.
    if len(components) == 1:
        versionless = accessions[0].split(".", 1)[0]
        if versionless == components[0]:
            return GenomeInfo(
                accession=accessions[0],
                description=description,
                components=ALL_CHROMOSOMES,
                source=source,
            )

    # Use the simple case of just using the requested components of the given
    # genome
    return GenomeInfo(
        accession=accessions[0],
        description=description,
        components=SelectedComponents.build(components),
        source=source,
    )


def find_description(xml: ET.Element) -> ty.Optional[str]:
    try:
        return proteome_value(xml, "pro:description").strip()
    except:
        return None


@lru_cache
def lineage_info(taxid: str) -> LineageInfo:
    response = requests.get(f"https://rest.uniprot.org/taxonomy/{taxid}.json")
    response.raise_for_status()
    data = response.json()
    parents = []
    should_hide = False
    for taxon in reversed(data["lineage"]):
        if taxon["taxonId"] == 2759:  # Hide in euks
            should_hide = True
        if should_hide and taxon.get("hidden", False):
            continue
        if taxon["scientificName"] == "cellular organisms":
            continue
        parents.append(taxon)
    tax_string = "; ".join(l["scientificName"] for l in parents)
    tax_string += "."

    common_name = data.get("commonName", None)
    if common_name == "":
        common_name = None

    return LineageInfo(
        ncbi_id=data["taxonId"],
        species=data["scientificName"],
        common_name=common_name,
        tax_string=tax_string,
    )


def proteome(xml: ET.Element) -> ProteomeInfo:
    upi = proteome_value(xml, "pro:upid")
    taxid = proteome_value(xml, "pro:taxonomy", convert=int)
    return ProteomeInfo(
        upi=upi,
        taxid=taxid,
        is_reference=proteome_value(xml, "pro:isReferenceProteome", convert=as_bool),
        is_representative=proteome_value(
            xml, "pro:isRepresentativeProteome", convert=as_bool
        ),
        proteome_description=find_description(xml),
        genome_info=genome_info(upi, xml),
        lineage_info=lineage_info(taxid),
    )


def proteomes(path: Path, ignore: ty.Set[str]) -> ty.Iterable[ProteomeInfo]:
    xml = ET.parse(path)
    proteomes = xml.getroot()
    seen = False
    for element in proteomes:
        parsed = proteome(element)
        if parsed.upi in ignore:
            LOGGER.info("Skipping upid %s as requested", parsed.upi)
            continue
        seen = True
        yield parsed

    if not seen:
        raise ValueError("Found no proteomes")
