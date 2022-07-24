#!/usr/bin/env python3

import json
import typing as ty
import logging
import xml.etree.ElementTree as ET
from pathlib import Path

import click


T = ty.TypeVar('T')

LOGGER = logging.getLogger(__name__)

NS = {
    "uni": "http://uniprot.org/uniprot",
    "pro": "http://uniprot.org/proteomes",
}


class GenomeInfo(ty.TypedDict):
    upi: str
    is_reference: bool
    is_representative: bool
    accession: str
    ids: ty.Union[ty.List[str], str]


def node_text(node: ty.Optional[ET.Element]) -> ty.Optional[str]:
    if node is not None and node.text is not None:
        return node.text
    return None


def all_matching_nodes(
    kind: str, query: str, root: ET.Element
) -> ty.Optional[ty.List[ty.Dict[str, str]]]:
    nodes = root.findall(query, namespaces=NS)
    if len(nodes) == 0:
        return None
    found = []
    for node in nodes:
        accession = node_text(node)
        if not accession:
            return None
        found.append({"kind": kind, "accession": accession})
    return found


def only_matching_node(kind: str, query: str, root: ET.Element):
    nodes = all_matching_nodes(kind, query, root)
    if nodes is None:
        return None
    if len(nodes) != 1:
        return None
    return nodes


def is_gca_only(accession: str, ids: ty.List[str]) -> bool:
    if len(ids) != 1:
        return False
    versionless = accession.split(".", 1)[0]
    return versionless == ids[0]


def xml_to_ncbi(root: ET.Element) -> ty.Optional[ty.Dict[str, ty.Any]]:
    result = only_matching_node(
        "ncbi", ".//pro:genomeAssembly/pro:genomeAssembly", root
    )
    if not result:
        return None
    if len(result) != 1:
        raise ValueError("Multiple GCA for a genome?")
    result = dict(result[0])
    if result["accession"][0:3] not in {"GCA", "GCF"}:
        return None

    components = all_matching_nodes("ena", ".//pro:component/pro:genomeAccession", root)
    if not components:
        result["ids"] = ["*"]
    else:
        accessions = [c["accession"] for c in components]
        if is_gca_only(result["accession"], accessions):
            result["ids"] = ["*"]
        else:
            result["ids"] = accessions

    return result


def xml_to_all_components(root: ET.Element) -> ty.Optional[ty.Dict[str, ty.Any]]:
    results = only_matching_node(
        "ena", './/pro:component[@name="Genome"]/pro:genomeAccession', root
    )
    if results is None:
        results = all_matching_nodes(
            "ena", ".//pro:component/pro:genomeAccession", root
        )
    if results is None:
        return None
    if len(results) == 1:
        if results[0]["accession"].startswith("GCA"):
            return {"kind": "ncbi", "accession": results[0]["accession"], "ids": ["*"]}

    return {"kind": "ena", "accession": None, "ids": [r["accession"] for r in results]}


def is_allowed(xml: ET.Element, to_skip: ty.Set[str]) -> bool:
    upi = xml.find("pro:upid", NS)
    if upi is None:
        raise ValueError("Invalid xml, no upid")
    upi = upi.text
    if upi is None:
        raise ValueError("Invalid xml, empty upi")

    return upi not in to_skip


def as_bool(raw: str) -> bool:
    if raw == 'true':
        return True
    if raw == 'false':
        return False
    raise ValueError("Cannot convert %s to a bool" % raw)


def proteome_value(xml: ET.Element, tag: str, convert: ty.Callable[[str], T]=str) -> T:
    value = xml.find(tag, NS)
    if value is None:
        raise ValueError(f"Invalid xml, no {tag}")
    value = value.text
    if value is None:
        raise ValueError(f"Invalid xml, empty {tag}")
    return convert(value)


def handle_proteome(xml: ET.Element) -> ty.Iterable[GenomeInfo]:
    methods = [
        xml_to_ncbi,
        xml_to_all_components,
    ]

    info = {
        "upi": proteome_value(xml, 'pro:upid'),
        "taxid": proteome_value(xml, "pro:taxonomy"),
        "is_reference": proteome_value(xml, "pro:isReferenceProteome", convert=as_bool),
        "is_representative": proteome_value(xml, "pro:isRepresentativeProteome", convert=as_bool),
    }

    for method in methods:
        result = method(xml)
        if not result:
            LOGGER.debug("Method %s failed on %s", xml_to_ncbi.__name__, info['upi'])
            continue
        LOGGER.debug("Method %s worked on %s",  xml_to_ncbi.__name__, info['upi'])
        result.update(info)
        yield result
        break
    else:
        raise ValueError("Failed to get accession for %s" % info["upi"])


@click.command()
@click.option("--ignorable", default="to-skip", type=click.File("r"))
@click.argument("summary", type=click.Path())
@click.argument("output", default=".", type=click.Path())
def main(summary, output, ignorable=None):
    logging.basicConfig(level=logging.INFO)
    to_skip = set()
    if ignorable:
        to_skip.update(l.strip() for l in ignorable)

    base = Path(output)
    handles = {}
    xml = ET.parse(summary)
    proteomes = xml.getroot()
    for proteome in proteomes:
        if not is_allowed(proteome, to_skip):
            LOGGER.info("Skipping proteome %s", proteome.find("pro:upid", NS).text)
            continue
        results = handle_proteome(proteome)
        for info in results:
            out = handles.get(info["kind"], None)
            if not out:
                out = (base / f"{info['kind']}.jsonl").open("w")
                handles[info["kind"]] = out
            json.dump(info, out)
            out.write("\n")

    for handle in handles.values():
        handle.flush()
        handle.close()


if __name__ == "__main__":
    main()
