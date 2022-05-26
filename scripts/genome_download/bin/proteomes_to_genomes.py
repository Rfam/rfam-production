#!/usr/bin/env python3

import csv
import typing as ty
import logging
import xml.etree.ElementTree as ET

import click
import requests

LOGGER = logging.getLogger(__name__)

NS = {
    'uni': "http://uniprot.org/uniprot",
    'pro': "http://uniprot.org/proteomes",
}

def proteome_xml(proteome: str) -> ET.Element:
    url = "http://www.uniprot.org/proteomes/%s.xml" % proteome
    response = requests.get(url)
    response.raise_for_status()

    return ET.fromstring(response.content)


def node_text(node: ty.Optional[ET.Element]) -> ty.Optional[str]:
    if node is not None and node.text is not None:
        return node.text
    return None


def all_matching_nodes(kind: str, query: str, root: ET.Element) -> ty.Optional[ty.List[ty.Tuple[str, str]]]:
    nodes = root.findall(query, namespaces=NS)
    if len(nodes) == 0:
        return None
    found = []
    for node in nodes:
        accession = node_text(node)
        if not accession:
            return None
        found.append((kind, accession))
    return found


def only_matching_node(kind: str, query: str, root: ET.Element):
    nodes = all_matching_nodes(kind, query, root)
    if nodes is None:
        return None
    if len(nodes) != 1:
        return None
    return nodes


def xml_to_ncbi(root: ET.Element) -> ty.Optional[ty.List[ty.Tuple[str, str]]]:
    result = only_matching_node('ncbi', './/pro:genomeAssembly/pro:genomeAssembly', root)
    if not result:
        return None
    if result[0][1][0:3] not in {'GCA', 'GCF'}:
        return None
    return result


def xml_to_accession(root: ET.Element) -> ty.Optional[ty.List[ty.Tuple[str, str]]]:
    return only_matching_node('other', './/pro:component[@name="Genome"]/pro:genomeAccession', root)


def xml_to_all_components(root: ET.Element) -> ty.Optional[ty.List[ty.Tuple[str, str]]]:
    return all_matching_nodes('other', './/pro:component/pro:genomeAccession', root)


def handle_proteome(xml: ET.Element, to_skip: ty.Set[str]) -> ty.Iterable[ty.Tuple[str, str, str]]:
    methods = [
        xml_to_ncbi,
        xml_to_accession,
        xml_to_all_components,
    ]
    upi = xml.find('pro:upid', NS)
    if upi is None:
        raise ValueError("Invalid xml, no upid")
    upi = upi.text
    if upi is None:
        raise ValueError("Invalid xml, empty upi")
    if upi in to_skip:
        return

    for method in methods:
        result = method(xml)
        if not result:
            continue

        for (kind, accession) in result:
            yield (upi, kind, accession)
        break
    else:
        raise ValueError("Failed to get accession for %s" % upi)


@click.command()
@click.option('--ignorable', default='to-skip', type=click.File('r'))
@click.argument('summary')
@click.argument('output', default='-', type=click.File('w'))
def main(summary, output, ignorable=None):
    to_skip = set()
    if ignorable:
        to_skip.update(l.strip() for l in ignorable)

    writer = csv.writer(output)
    xml = ET.parse(summary)
    proteomes = xml.getroot()
    for proteome in proteomes:
        info = handle_proteome(proteome, to_skip)
        writer.writerows(info)

if __name__ == '__main__':
    main()
