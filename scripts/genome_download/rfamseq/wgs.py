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

import re
import enum
import typing as ty
import logging

from attrs import frozen
import cattrs

from rfamseq import ena, wget, ncbi

LOGGER = logging.getLogger(__name__)


PATTERN = re.compile(r"^([A-Z]{4,6}\d{2}S?)\d+$")


def wgs_endpoint(raw: str) -> ty.Tuple[str, int]:
    if not (match := re.match(PATTERN, raw)):
        raise ValueError(f"Cannot parse {raw}")

    prefix = match.group(1)
    id = raw[len(prefix):]
    minimal = re.sub('^0+', '', id) or '0'
    return (prefix, int(minimal))


def wgs_id(length: int, prefix: str, index: int) -> str:
    padding = length - len(prefix)
    pattern = f"{prefix}%0{padding}i"
    return pattern % index


@enum.unique
class WgsSequenceKind(enum.Enum):
    SCAFFOLD = "WGS_SCAFLD"
    SEQUENCE = "WGS"


@frozen
class WgsInfo:
    wgs_id: str
    kind: WgsSequenceKind
    start: int
    stop: int
    id_length: int

    @classmethod
    def from_endpoint(cls, kind: WgsSequenceKind, endpoint: str) -> WgsInfo:
        prefix, start = wgs_endpoint(endpoint)
        return cls(wgs_id=prefix, kind=kind, start=start, stop=start,
            id_length=len(endpoint))

    @classmethod
    def from_range(cls, kind: WgsSequenceKind, range: str) -> WgsInfo:
        raw_start, raw_stop = range.split('-', 1)
        prefix1, start = wgs_endpoint(raw_start)
        prefix2, stop = wgs_endpoint(raw_stop)

        assert prefix1 == prefix2, f"Cannot create range across prefixes {range}"
        assert len(raw_start) == len(raw_stop), f"Cannot create range with bad lengths {range}"

        return cls(wgs_id=prefix1, kind=kind, start=start, stop=stop,
                id_length=len(raw_start))

    @classmethod
    def build(cls, kind: WgsSequenceKind, range: str) -> WgsInfo:
        if '-' not in range:
            return cls.from_endpoint(kind, range)
        return cls.from_range(kind, range)

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
    ena_info: ty.List[WgsInfo]
    ncbi_ids: ty.List[str]

    def ids(self) -> ty.Iterable[str]:
        for wgs in self.ena_info:
            yield from wgs.ids()

        for id in self.ncbi_ids:
            yield id

    def record_ids(self) -> ty.Set[str]:
        return {wgs.record_id() for wgs in self.ena_info}

    def __contains__(self, accession: str) -> bool:
        for id in self.ids():
            if id == accession:
                return True
        return False



def resolve_ena_wgs(accession: str) -> ty.List[WgsInfo]:
    LOGGER.info("Fetching EMBL formatted file for %s", accession)
    url = ena.ENA_EMBL_URL.format(accession=accession)
    info = []
    with wget.wget(url) as handle:
        for line in handle:
            prefix = line[0:3]
            try:
                kind = cattrs.structure(prefix, WgsSequenceKind)
            except:
                continue
            info.append(WgsInfo.build(kind, line[3:].strip()))
    return info


def resolve_wgs(accession: str) -> ty.Optional[WgsSummary]:
    ena_info = resolve_ena_wgs(accession)
    ncbi_ids = ncbi.resolve_wgs(accession)
    if not ncbi_ids and not ena_info:
        return None

    return WgsSummary(ena_info=ena_info, ncbi_ids=(ncbi_ids or []))
