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

import collections as coll
import enum
import logging
import re
import typing as ty

import cattrs
from attrs import frozen

from rfamseq import ena, ncbi, wget

LOGGER = logging.getLogger(__name__)


PATTERN = re.compile(r"^([A-Z]{4,6}\d{2}S?)\d+(\.\d+)?$")

CONTIG_PATTERN = re.compile(r"(^[A-Z]+)(\d+)$")


class InvalidWgsAccession(Exception):
    """
    Raised if the WGS accession is not formatted as expected.
    """


def looks_like_wgs_accession(raw: str) -> bool:
    try:
        wgs_endpoint(raw)
        return True
    except InvalidWgsAccession:
        return False


def wgs_endpoint(raw: str) -> ty.Tuple[str, int]:
    if not (match := re.match(PATTERN, raw)):
        raise InvalidWgsAccession(raw)

    prefix = match.group(1)
    version_suffix = match.group(2)

    id = raw[len(prefix) :]
    if version_suffix:
        id = id.replace(version_suffix, "")

    minimal = re.sub("^0+", "", id) or "0"

    # TODO: Check that the versions are the same. The last two numbers after
    #       the letters in a WGS id are versions. These should match, ie
    #       WWJG01000002.1 has version 1. The 01 after WWJG and .1 are the
    #       same.
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
        if "-" not in raw:
            raw_start = raw_stop = raw
        else:
            raw_start, raw_stop = raw.split("-", 1)
        if not (match := re.match(CONTIG_PATTERN, raw_start)):
            raise ValueError(f"Cannot parse {raw}")
        prefix = match.group(1)
        start = match.group(2)

        if not (match := re.match(CONTIG_PATTERN, raw_stop)):
            raise ValueError(f"Cannot parse {raw}")

        assert prefix == match.group(1), f"Mismatch prefix {raw}"
        stop = match.group(2)

        return cls(prefix=prefix, start=start, stop=stop)

    def as_range(self) -> str:
        return f"{self.prefix}{self.start}-{self.prefix}{self.stop}"

    def ids(self) -> ty.Iterable[str]:
        num_leading_zeros = len(self.start) - len(str(int(self.start)))
        for index in range(int(self.start), int(self.stop) + 1):
            yield f"{self.prefix}{str(index).zfill(num_leading_zeros + len(str(index)))}"


@enum.unique
class WgsSequenceKind(enum.Enum):
    """
    This represents the different assembly levels that a WGS sequence can have.

    Note the ordering of these values matters. They are arranged largest to
    smallest which is used later when deciding which are the largest chunks to
    use.
    """

    CONTIG = "CON"
    SCAFFOLD = "WGS_SCAFLD"
    SEQUENCE = "WGS"


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
    contigs: ty.List[ContigInfo]
    sequences: ty.List[WgsSequence]
    ncbi_ids: ty.List[str]

    @property
    def wgs_prefix(self) -> str:
        return self.wgs_id[0:4]

    @property
    def wgs_version(self) -> int:
        return int(self.wgs_id[4:6])

    def ids(self) -> ty.Iterable[str]:
        for contig in self.contigs:
            yield from contig.ids()

        for wgs in self.sequences:
            yield from wgs.ids()

        for id in self.ncbi_ids:
            yield id

    def largest_ids(self) -> ty.Iterable[str]:
        if self.contigs:
            for contig in self.contigs:
                yield from contig.ids()
        elif self.sequences:
            mapped: ty.Dict[WgsSequenceKind, ty.List[WgsSequence]] = coll.defaultdict(
                list
            )
            for sequence in self.sequences:
                mapped[sequence.kind].append(sequence)

            for kind in WgsSequenceKind:
                seen = False
                for wgs_sequence in mapped.get(kind, []):
                    yield from wgs_sequence.ids()
                    seen = True
                if seen:
                    break
            else:
                raise ValueError("Failed to produce and wgs sequences")

        elif self.ncbi_ids:
            yield from self.ncbi_ids
        else:
            raise ValueError("Somehow failed to have any ids in a wgs set")

    def record_ids(self) -> ty.Set[str]:
        return {wgs.record_id() for wgs in self.sequences}

    def id_matches(self, raw: str, within_one_version=False):
        max_diff = int(within_one_version)
        try:
            prefix, version = wgs_endpoint(raw)
            return (
                prefix[0:4] == self.wgs_prefix
                and abs(version - int(self.wgs_version)) == max_diff
            )
        except InvalidWgsAccession:
            return False

    def __contains__(self, accession: str) -> bool:
        for id in self.ids():
            if id == accession:
                return True
        return False


def resolve_ena_wgs(
    accession: str,
) -> ty.Tuple[ty.List[ContigInfo], ty.List[WgsSequence]]:
    LOGGER.info("Fetching EMBL formatted file for %s", accession)
    url = ena.ENA_EMBL_URL.format(accession=accession)
    contigs: ty.List[ContigInfo] = []
    info: ty.List[WgsSequence] = []
    try:
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
                    if kind is not WgsSequenceKind.CONTIG:
                        raise err

                    value = line[3:].strip()
                    parts = value.split(";")
                    for part in parts:
                        contigs.append(ContigInfo.build(part))
    except wget.FetchError as e:
        LOGGER.debug(e)

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
