# -*- coding: utf-8 -*-

# Copyright [2009-2024] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import typing as ty
import xml.etree.ElementTree as ET
from datetime import date
from enum import Enum, unique
from pathlib import Path

from attr import NOTHING, AttrsInstance, field, fields
from loguru import logger
from rfam_export.context import Context

MY_TYPE_METADATA = "__RFAM_EXPORT_METADATA"


class NoEntries(Exception):
    """Raised when trying to build an XML object and there are no entries in
    it.
    """


def assert_never(x: ty.NoReturn) -> ty.NoReturn:
    """A way to assert both statically and at run-time that something will
    never happen.

    This is used in match statements to say a type cannot occur. When typing
    checking if a branch is not covered, this will cause an error, at run-time
    this will raise an error.

    The following example should raise an error when run and family type
    checking.

    ```
    >>> a = 4
    >>> match a
        case str():
            print("A string")
        case _:
            assert_never(a)
    ...
    AssertionError "Unhandled type: int"
    ```
    """
    assert False, "Unhandled type: {}".format(type(x).__name__)


@unique
class EntryKind(Enum):
    ELEMENT = "element"
    ADDITIONAL = "additional_fields"
    CROSS_REFERENCE = "cross_references"


@unique
class EntryType(Enum):
    CLAN = "Clan"
    FAMILY = "Family"
    GENOME = "Genome"
    MOTIF = "Motif"
    SEQUENCE = "Sequence"


def element(
    name=None,
    default=NOTHING,
    validator=None,
    repr=True,
    eq=True,
    order=None,
    hash=None,
    init=True,
    metadata=None,
    converter=None,
) -> ty.Any:
    metadata = metadata or {}
    metadata[MY_TYPE_METADATA] = {
        "kind": EntryKind.ELEMENT,
        "name": name,
    }
    return field(
        default=default,
        validator=validator,
        repr=repr,
        eq=eq,
        order=order,
        hash=hash,
        init=init,
        metadata=metadata,
        converter=converter,
    )


def additional_field(
    name=None,
    default=NOTHING,
    validator=None,
    repr=True,
    eq=True,
    order=None,
    hash=None,
    init=True,
    metadata=None,
    converter=None,
) -> ty.Any:
    metadata = metadata or {}
    metadata[MY_TYPE_METADATA] = {
        "kind": EntryKind.ADDITIONAL,
        "name": name,
    }
    return field(
        default=default,
        validator=validator,
        repr=repr,
        eq=eq,
        order=order,
        hash=hash,
        init=init,
        metadata=metadata,
        converter=converter,
    )


def cross_reference(
    dbname,
    name=None,
    default=NOTHING,
    validator=None,
    repr=True,
    eq=True,
    order=None,
    hash=None,
    init=True,
    metadata=None,
    converter=None,
) -> ty.Any:
    metadata = metadata or {}
    metadata[MY_TYPE_METADATA] = {
        "kind": EntryKind.CROSS_REFERENCE,
        "dbname": dbname,
        "name": name,
    }
    return field(
        default=default,
        validator=validator,
        repr=repr,
        eq=eq,
        order=order,
        hash=hash,
        init=init,
        metadata=metadata,
        converter=converter,
    )


SimpleValue = ty.Union[str, float, Enum]


def as_values(
    raw: SimpleValue | ty.List[SimpleValue] | ty.Set[SimpleValue],
) -> ty.List[str]:
    if isinstance(raw, str):
        return [raw]
    if isinstance(raw, (int, float)):
        return [str(raw)]
    if isinstance(raw, Enum):
        return as_values(raw.value)
    if isinstance(raw, (list, set)):
        result = []
        for v in raw:
            result.extend(as_values(v))
        return result
    assert_never(raw)


def as_xml(raw: AttrsInstance) -> ET.Element:
    """Convert a `attrs` decorated class into an `Element` suitable for EBI's
    search XML.

    The object must have one field for each element that is to be added to the
    XML. The fields must be created with `cross_reference`, `additional_field`
    or `element` to indicate where in the XML the each field should go.

    For field annotated as `additional_field` this will use `Enum` values or
    each entry in a list.
    """

    entry = ET.Element("entry", attrib={"id": raw.entry_id})
    cross_references = ET.Element("cross_references")
    additional_fields = ET.Element("additional_fields")
    for field in fields(type(raw)):
        kind = field.metadata[MY_TYPE_METADATA]["kind"]
        xml_name = field.metadata[MY_TYPE_METADATA].get("name", None)
        if xml_name is None:
            xml_name = field.name
        if not xml_name:
            raise ValueError(f"Cannot find name for XML attribue from {field}")

        raw_values = getattr(raw, field.name)
        if raw_values is None:
            continue
        values = as_values(raw_values)
        for value in values:
            match kind:
                case EntryKind.ELEMENT:
                    ET.SubElement(entry, xml_name).text = value
                case EntryKind.ADDITIONAL:
                    attr = {"name": xml_name}
                    ET.SubElement(additional_fields, "field", attr).text = value
                case EntryKind.CROSS_REFERENCE:
                    attrs = {
                        "dbkey": value,
                        "dbname": field.metadata[MY_TYPE_METADATA]["dbname"],
                    }
                    ET.SubElement(cross_references, "ref", attrs)
                case _:
                    assert_never(kind)

    dates = ET.SubElement(entry, "dates")
    today = date.today().strftime("%d %b %Y")
    ET.SubElement(dates, "date", {"type": "created", "value": today})
    ET.SubElement(dates, "date", {"type": "updated", "value": today})
    entry.append(additional_fields)
    entry.append(cross_references)
    return entry


def write_file(context: Context, entries: ty.Iterable[ET.Element], handle: ty.TextIO):
    """Generate the XML document with all the given entries in the handle.

    This assumes this is an Rfam database file for the context's `rfam_version`.
    This will use a small amount of memory as it does not build the entire
    object up at once, and instead writes each entry seperatly. This does lead
    a slightly ugly XML document, but it is parsable and correct.

    Each element of `entries` should be an xml Entry element to write.

    It is an error if entries contains no elements to write and will raise
    `NoEntries`.
    """

    logger.debug("Writing xml header")
    handle.write("<database>\n")
    handle.write("<name>Rfam</name>\n")
    handle.write(
        "<description>A database for non-protein coding RNA families</description>\n"
    )
    handle.write(f"<release>{context.rfam_version}</release>\n")
    timestamp = date.today().strftime("%d/%m/%Y")
    handle.write(f"<release_date>{timestamp}</release_date>\n")

    handle.write("<entries>\n")
    count = 0
    for entry in entries:
        handle.write(ET.tostring(entry, encoding="unicode"))
        handle.write("\n")
        count += 1
    handle.write("</entries>\n")
    logger.debug("Wrote {} entries", count)

    if not count:
        raise NoEntries()

    handle.write(f"<entry_count>{count}</entry_count>\n")
    handle.write("</database>")
