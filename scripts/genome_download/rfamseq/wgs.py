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
import re
import typing as ty

import cattrs
from attrs import frozen

from rfamseq import ena, ncbi, wget

LOGGER = logging.getLogger(__name__)


PATTERN = re.compile(r"^([A-Z]{4,6}\d{2}S?)\d+$")

CONTIG_PATTERN = re.compile(r"(^[A-Z]+)(\d+)$")


class InvalidWgsAccession(Exception):
    """
    Raised if the WGS accession is not formatted as expected.
    """


def wgs_endpoint(raw: str) -> ty.Tuple[str, int]:
    if not (match := re.match(PATTERN, raw)):
        raise InvalidWgsAccession(raw)

    prefix = match.group(1)
    id = raw[len(prefix) :]
    minimal = re.sub("^0+", "", id) or "0"
    return (prefix, int(minimal))


def wgs_id(length: int, prefix: str, index: int) -> str:
    padding = length - len(prefix)
    pattern = f"{prefix}%0{padding}i"
    return pattern % index


@frozen
class ContigInfo:
    prefix: str
    start: int
    stop: int

    @classmethod
    def build(cls, raw: str) -> ContigInfo:
        raw_start, raw_stop = raw.split("-", 1)
        if not (match := re.match(CONTIG_PATTERN, raw_start)):
            raise ValueError(f"Cannot parse {raw}")
        prefix = match.group(1)
        start = match.group(2)

        if not (match := re.match(CONTIG_PATTERN, raw_stop)):
            raise ValueError(f"Cannot parse {raw}")

        assert prefix == match.group(1), f"Mismatch prefix {raw}"
        stop = match.group(2)

        return cls(prefix=prefix, start=int(start), stop=int(stop))

    def as_range(self) -> str:
        return f"{self.prefix}{self.start}-{self.prefix}{self.stop}"

    def ids(self) -> ty.Iterable[str]:
        for index in range(self.start, self.stop + 1):
            yield f"{self.prefix}{index}"


@enum.unique
class WgsSequenceKind(enum.Enum):
    SCAFFOLD = "WGS_SCAFLD"
    SEQUENCE = "WGS"
    CONTIG = "CON"


@frozen
class WgsSequence:
    wgs_id: str
    kind: WgsSequenceKind
    start: int
    stop: int
    id_length: int

    @classmethod
    def from_endpoint(cls, kind: WgsSequenceKind, endpoint: str) -> WgsSequence:
        prefix, start = wgs_endpoint(endpoint)
        return cls(
            wgs_id=prefix, kind=kind, start=start, stop=start, id_length=len(endpoint)
        )

    @classmethod
    def from_range(cls, kind: WgsSequenceKind, range: str) -> WgsSequence:
        raw_start, raw_stop = range.split("-", 1)
        prefix1, start = wgs_endpoint(raw_start)
        prefix2, stop = wgs_endpoint(raw_stop)

        assert prefix1 == prefix2, f"Cannot create range across prefixes {range}"
        assert len(raw_start) == len(
            raw_stop
        ), f"Cannot create range with bad lengths {range}"

        return cls(
            wgs_id=prefix1, kind=kind, start=start, stop=stop, id_length=len(raw_start)
        )

    @classmethod
    def build(cls, kind: WgsSequenceKind, range: str) -> WgsSequence:
        if "-" not in range:
            return cls.from_endpoint(kind, range)
        return cls.from_range(kind, range)

    def within_one_version(self, accession: str) -> bool:
        accession = re.sub(r"0+$", "", accession)
        pattern = re.compile(r"([A-Z]+)0+(\d+)S?$")
        if not (acc_match := re.search(pattern, accession)):
            return False
        if not (match := re.search(pattern, self.wgs_id)):
            return False

        if not match.group(1) == acc_match.group(1):
            return False

        current = int(match.group(2))
        suggested = int(acc_match.group(2))

        return abs(current - suggested) == 1

    def record_id(self) -> str:
        return wgs_id(self.id_length, self.wgs_id, 0)

    def as_range(self) -> str:
        start = wgs_id(self.id_length, self.wgs_id, self.start)
        stop = wgs_id(self.id_length, self.wgs_id, self.stop)
        if self.start == self.stop:
            return start
        return f"{start}-{stop}"

    def ids(self) -> ty.Iterable[str]:
        for index in range(self.start, self.stop + 1):
            yield wgs_id(self.id_length, self.wgs_id, index)

    def __contains__(self, accession: str) -> bool:
        for id in self.ids():
            if id == accession:
                return True
        return False


@frozen
class WgsSummary:
    wgs_id: str
    contigs: ty.Optional[ContigInfo]
    sequences: ty.List[WgsSequence]
    ncbi_ids: ty.List[str]

    def ids(self) -> ty.Iterable[str]:
        if self.contigs:
            yield from self.contigs.ids()

        for wgs in self.sequences:
            yield from wgs.ids()

        for id in self.ncbi_ids:
            yield id

    def record_ids(self) -> ty.Set[str]:
        return {wgs.record_id() for wgs in self.sequences}

    def within_one_version(self, accession: str) -> bool:
        for info in self.sequences:
            if isinstance(info, ContigInfo):
                continue
            if info.within_one_version(accession):
                return True
        return False

    def __contains__(self, accession: str) -> bool:
        for id in self.ids():
            if id == accession:
                return True
        return False


def resolve_ena_wgs(
    accession: str,
) -> ty.Tuple[ty.Optional[ContigInfo], ty.List[WgsSequence]]:
    LOGGER.info("Fetching EMBL formatted file for %s", accession)
    url = ena.ENA_EMBL_URL.format(accession=accession)
    contigs = None
    info: ty.List[WgsSequence] = []
    with wget.wget(url) as handle:
        for line in handle:
            prefix = line[0:3]
            try:
                kind = cattrs.structure(prefix, WgsSequenceKind)
            except:
                continue
            try:
                info.append(WgsSequence.build(kind, line[3:].strip()))
            except InvalidWgsAccession as err:
                if kind is WgsSequenceKind.CONTIG:
                    contigs = ContigInfo.build(line[3:].strip())
                else:
                    raise err

    return (contigs, info)


def resolve_wgs(accession: str) -> ty.Optional[WgsSummary]:
    (contigs, ena_info) = resolve_ena_wgs(accession)
    ncbi_ids = ncbi.resolve_wgs(accession)
    if not ncbi_ids and not ena_info and not contigs:
        return None

    wgs_id = ena_info[0].wgs_id
    return WgsSummary(
        wgs_id=wgs_id, contigs=contigs, sequences=ena_info, ncbi_ids=(ncbi_ids or [])
    )
