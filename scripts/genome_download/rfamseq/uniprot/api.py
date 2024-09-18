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
import io
import logging
import typing as ty
import xml.etree.ElementTree as ET
from pathlib import Path

import cattrs
import requests
from attrs import frozen

from rfamseq import wgs
from rfamseq.accession import Accession
from rfamseq.utils import assert_never, batched

LOGGER = logging.getLogger(__name__)

API_TEMPLATE = "https://rest.uniprot.org/proteomes/stream"

NS = {
    "uni": "http://uniprot.org/uniprot",
    "pro": "http://uniprot.org/proteome",
}

T = ty.TypeVar("T")


@frozen
class All:
    all: bool


@frozen(hash=True)
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
        """
        Create a new SelectedComponent. This requires the list of components
        be not empty.
        """
        if not components:
            raise ValueError("Cannot build empty SelectedComponents")
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

    def matching(self, component: Component) -> ty.List[Component]:
        match component:
            case Unplaced():
                if self.unplaced:
                    return [UNPLACED]
                return []

            case Accession():
                return [a for a in self.accessions if component.matches(a)]

            case wgs.WgsPrefix():
                matches: ty.List[Component] = []
                for wgs_set in self.wgs_sets:
                    if wgs_set.matches(component, within_one_version=True):
                        matches.append(wgs_set)
                for wgs_sequence in self.wgs_sequences:
                    wgs_set = wgs_sequence.prefix
                    if wgs_set.matches(component, within_one_version=True):
                        matches.append(wgs_sequence)
                return matches

            case wgs.WgsSequenceId():
                matches: ty.List[Component] = []
                for wgs_set in self.wgs_sets:
                    if wgs_set.matches(component.prefix, within_one_version=True):
                        matches.append(wgs_set)
                for wgs_sequence in self.wgs_sequences:
                    if wgs_sequence.matches(component, within_one_version=True):
                        matches.append(wgs_sequence)
                return matches

            case _:
                assert_never(component)

    def includes_unplaced(self) -> bool:
        return self.unplaced

    def includes_wgs(self):
        return bool(self.wgs_sets) or bool(self.wgs_sequences)


Components = ty.Union[All, SelectedComponents]


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


def proteome_components(
    upid: str, root: ET.Element
) -> ty.Tuple[bool, ty.List[Component], ty.Optional[str]]:
    description = None
    saw_genome = False
    components = []
    for component in root.findall("pro:component", NS):
        comp_accs = all_matching_node_text("pro:genomeAccession", component)
        name = component.attrib.get("name", "")
        if "unplaced" in name.lower():
            components.append(UNPLACED)
            if not comp_accs:
                continue
        if name == "Genome":
            saw_genome = True
            descriptions = all_matching_node_text("pro:description", component)
            if descriptions and len(descriptions) == 1:
                description = descriptions[0]

        if not comp_accs:
            raise ValueError(f"Missing component accession of {upid}")

        for raw_acc in comp_accs:
            if wgs_component := wgs.parse_wgs_accession(raw_acc):
                components.append(wgs_component)
            else:
                components.append(Accession.build(raw_acc))
    return (saw_genome, components, description)


def genome_info(upid: str, root: ET.Element) -> GenomeInfo:
    accessions = all_matching_node_text(
        ".//pro:genomeAssembly/pro:genomeAssembly", root
    )
    if accessions and len(accessions) > 1:
        raise ValueError(f"Cannot handle >1 genome accession in {upid}")

    source = None
    if sources := all_matching_node_text(
        ".//pro:genomeAssembly/pro:genomeAssemblySource", root
    ):
        LOGGER.debug("Found genome sources: %s", sources)
        if len(sources) == 1:
            source = cattrs.structure(sources[0].upper(), GenomeSource)

    saw_genome, components, description = proteome_components(upid, root)

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
    if len(components) == 1 and isinstance(components[0], Accession):
        if Accession.build(accessions[0]).matches(components[0]):
            return GenomeInfo(
                accession=accessions[0],
                description=description,
                components=ALL_CHROMOSOMES,
                source=source,
            )

    # If given a GCA/GCF accession, just use all components in it. This
    # simplifies the downloading logic and possible erorr cases quite a bit.
    # This does increase the amount we download and does mean we use 'extra'
    # sequences but I think that is ok. It also doesn't appear, at least by
    # manual checking, make us miss any sequences.
    if accessions[0].startswith("GCA_") or accessions[0].startswith("GCF_"):
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


def proteomes(
    path: Path | io.StringIO, ignore: ty.Set[str]
) -> ty.Iterable[ProteomeInfo]:
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


# https://rest.uniprot.org/proteomes/stream?format=xml&query=%28%28upid%3AUP000001579%29+OR+%28upid%3AUP000000251%29%29
def fetch_proteomes(ids: ty.Iterable[str], batch_size=10) -> ty.Iterable[ET.Element]:
    session = requests.Session()
    for batch in batched(ids, batch_size):
        terms = [f"(upid:{id})" for id in batch]
        query = " OR ".join(terms)
        response = session.get(
            API_TEMPLATE, params={"compressed": "true", "format": "xml", "query": query}
        )
        tree = ET.fromstring(response.text)
        yield from tree.findall("proteome", namespaces=NS)