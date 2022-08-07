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

import logging
import typing as ty
import xml.etree.ElementTree as ET
from pathlib import Path

import cattrs
from attrs import define, field

LOGGER = logging.getLogger(__name__)

NS = {
    "uni": "http://uniprot.org/uniprot",
    "pro": "http://uniprot.org/proteome",
}

T = ty.TypeVar("T")

MaybeStrs = ty.Optional[ty.List[str]]


@define
class All:
    all: bool


@define
class Unplaced:
    unplaced: bool


UNPLACED = Unplaced(unplaced=True)

ALL_CHROMOSOMES = All(all=True)

Component = ty.Union[str, Unplaced]


def structure_component(v: ty.Any, _) -> Component:
    if isinstance(v, str):
        return v
    if isinstance(v, dict):
        return Unplaced(**v)
    raise ValueError(f"Cannot unstructure {v} to Component")


cattrs.register_structure_hook(Component, structure_component)


@define
class SelectedComponents:
    accessions: ty.List[Component] = field()

    @accessions.validator
    def _check_accessions(self, attribute, value):
        assert value, "Cannot create empty selected components"
        for v in value:
            assert isinstance(v, Component)

    def __len__(self) -> int:
        return len(self.accessions)

    def __contains__(self, accession: Component) -> bool:
        return accession in self.accessions

    def __iter__(self):
        return iter(self.accessions)


Components = ty.Union[All, SelectedComponents]


@define
class GenomeInfo:
    accession: ty.Optional[str]
    description: ty.Optional[str]
    components: Components

    @property
    def version(self) -> ty.Optional[str]:
        if not self.accession or "." not in self.accession:
            return None
        return self.accession.split(".", 1)[1]


@define
class ProteomeInfo:
    upi: str
    taxid: str
    is_reference: bool
    is_representative: bool
    proteome_description: ty.Optional[str]
    genome_info: GenomeInfo

    @property
    def description(self) -> ty.Optional[str]:
        if self.proteome_description:
            return self.proteome_description
        return genome_info.description


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

        components.extend(comp_accs)

    # Handle the case where there is not overarching genome assembly to work
    # with and just a series of components to fetch.
    if not accessions:

        # If the only component is a GCA_ sequence then use that and all
        # components.
        if len(components) == 1:
            possible_gca = components[0]
            if isinstance(possible_gca, str) and possible_gca.startswith("GCA_"):
                return GenomeInfo(
                    accession=possible_gca,
                    description=description,
                    components=ALL_CHROMOSOMES,
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
                )

        # Otherwise only use the given components
        return GenomeInfo(
            accession=None,
            description=description,
            components=SelectedComponents(components),
        )

    assert len(accessions) == 1, f"Multiple genomes, impossibility in {upid}"

    # If not given any components then assume we want all sequences
    if not components:
        return GenomeInfo(
            accession=accessions[0], description=description, components=ALL_CHROMOSOMES
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
            )

    # Use the simple case of just using the requested components of the given
    # genome
    return GenomeInfo(
        accession=accessions[0],
        description=description,
        components=SelectedComponents(components),
    )


def find_description(xml: ET.Element) -> ty.Optional[str]:
    try:
        return proteome_value(xml, "pro:description")
    except:
        return None


def proteome(xml: ET.Element) -> ProteomeInfo:
    upi = proteome_value(xml, "pro:upid")
    ginfo = genome_info(upi, xml)
    return ProteomeInfo(
        upi=upi,
        taxid=proteome_value(xml, "pro:taxonomy"),
        is_reference=proteome_value(xml, "pro:isReferenceProteome", convert=as_bool),
        is_representative=proteome_value(
            xml, "pro:isRepresentativeProteome", convert=as_bool
        ),
        proteome_description=find_description(xml),
        genome_info=ginfo,
    )


def proteomes(path: Path, ignore: ty.Set[str]) -> ty.Iterable[ProteomeInfo]:
    xml = ET.parse(path)
    proteomes = xml.getroot()
    for element in proteomes:
        try:
            upid = element.find("pro:upid", NS).text
        except:
            upid = None
        if upid in ignore:
            LOGGER.info("Skipping upid as requested")
            continue
        yield proteome(element)