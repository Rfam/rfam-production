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
from attrs import field, frozen

from rfamseq import ena, ncbi, wget

LOGGER = logging.getLogger(__name__)

PATTERN = re.compile(r"^([A-Z]{4,6}\d{2}S?)(\d+)(\.\d+)?$")

PREFIX_PATTERN = re.compile(r"^([A-Z]{4,6})(\d{2}S?)$")

CONTIG_PATTERN = re.compile(r"(^[A-Z]+)(\d+)$")


class InvalidWgsAccession(Exception):
    """
    Raised if the WGS accession is not formatted as expected.
    """


class InvalidWgsPrefix(Exception):
    """
    Raised if the WGS prefix is not formatted as expected.
    """


@frozen
class WgsPrefix:
    """
    This represents the parsed prefix of WGS ids. WGS ids start with a unique
    identifier for the WGS then a version number, 0 padded to at least 2 digits
    and then an optional S, I assume to mean scaffold. This represents that
    string as more useful object that can be compared against each other.
    """

    wgs_id: str
    wgs_version: int
    scaffold: bool
    length: int

    @classmethod
    def build(cls, raw: str) -> WgsPrefix:
        """
        The method to build WgsPrefix object. This requires that the given raw
        id is a valid WGS prefix id with no extra text.

        >>> WgsPrefix.build("JABDTM02")
        WgsPrefix(wgs_id="JABDTM", wgs_version=2, scaffold=False, length=8)
        >>> WgsPrefix.build("ACTP02")
        WgsPrefix(wgs_id="ACTP", wgs_version=2, scaffold=False, length=6)
        >>> WgsPrefix.build("ALWZ04S")
        WgsPrefix(wgs_id="ALWZ", wgs_version=4, scaffold=True, length=7)
        """
        if not (match := re.match(PREFIX_PATTERN, raw)):
            raise InvalidWgsPrefix(raw)

        scaffold = raw.endswith("S")
        raw_version = re.sub("^0+", "", match.group(2))
        raw_version = re.sub("S$", "", raw_version)
        return cls(
            wgs_id=match.group(1),
            wgs_version=int(raw_version),
            scaffold=scaffold,
            length=len(raw),
        )

    def to_wgs_string(self) -> str:
        """
        Convert the object into a wgs prefix string. Basically the reverse of
        building the object.

        >>> WgsPrefix(wgs_id="JABDTM", wgs_version=2, scaffold=False, length=8).wgs_prefix
        "JABDTM02"
        >>> WgsPrefix(wgs_id="ACTP", wgs_version=2, scaffold=False, length=6).wgs_prefix
        "ACTP02"
        >>> WgsPrefix(wgs_id="ALWZ", wgs_version=2, scaffold=True, length=7).wgs_prefix
        "ALWZ04S"
        """
        version_length = self.length - (len(self.wgs_id) + int(self.scaffold))
        current = f"{self.wgs_id}{self.wgs_version:0{version_length}}"
        if self.scaffold:
            current += "S"
        return current

    def matches(self, other: WgsPrefix, within_one_version=False) -> int:
        """
        Fuzzy compare two WgsPrefixes. This will allow two prefixes to match if
        they only differ by one version, if allowed. This is sometimes useful
        if the old wgs set has disappeared, which can happen.
        """
        if not within_one_version:
            return self == other

        return (
            self.wgs_id == other.wgs_id
            and abs(self.wgs_version - other.wgs_version) <= 1
            and self.scaffold == other.scaffold
        )


@frozen
class WgsSequenceId:
    """
    This represents a single WGS sequence. This is basically a prefix with an
    index to indicate the sequence number. This does not prevent the index from
    being 0, which is commonly used to refer to the set. There is also an
    optional sequence version at the end. These seem to only show up in some
    ENA fasta files and may not be needed according the NCBI, but this handles
    them anyway.
    """

    prefix: WgsPrefix
    sequence_index: int
    sequence_version: ty.Optional[str]
    length: int

    @classmethod
    def build(cls, raw: str) -> WgsSequenceId:
        """
        Create a new WgsSequenceId. This will parse the given raw string to
        produce a WgsSequenceId object. This will fail with InvalidWgsAccession
        if it is not in the expected format.

        >>> WgsSequenceId.build("ALWZ04S3033285")
        WgsSequenceId(prefix=WgsPrefix(wgs_id="ALWZ", wgs_version=4, scaffold=True, length=7), sequence_index=3033285, sequence_version=None, length=14)
        """
        if not (match := re.match(PATTERN, raw)):
            raise InvalidWgsAccession(raw)

        wgs_prefix = match.group(1)
        sequence_version = match.group(3)
        if not sequence_version:
            sequence_version = None
        else:
            sequence_version = sequence_version[1:]

        return cls(
            prefix=WgsPrefix.build(wgs_prefix),
            sequence_index=int(match.group(2)),
            sequence_version=sequence_version,
            length=len(raw),
        )

    def is_wgs_set_reference(self) -> bool:
        """
        Check if this WgsSequenceId is actually a reference to the entire WGS
        set. This happens when the sequence index is 0. Sometimes ids are given
        with all zeros but these generally mean the entire WGS set.

        >>> WgsSequenceId(prefix=WgsPrefix(wgs_id="ALWZ", wgs_version=4, scaffold=True, length=7), sequence_index=3033285, sequence_version=None, length=14).is_wgs_set_reference()
        False
        >>> WgsSequenceId(prefix=WgsPrefix(wgs_id="ALWZ", wgs_version=4, scaffold=True, length=7), sequence_index=0, sequence_version=None, length=14).is_wgs_set_reference()
        True
        """
        return self.sequence_index == 0

    def to_wgs_string(self) -> str:
        """
        Convert this parsed id back into a string.

        >>> WgsSequenceId(prefix=WgsPrefix(wgs_id="ALWZ", wgs_version=4, scaffold=True, length=7), sequence_index=3033285, sequence_version=None, length=14).to_wgs_string()
        "ALWZ04S3033285"
        """
        prefix = self.prefix.to_wgs_string()
        suffix = ""
        if self.sequence_version is not None:
            suffix = f".{self.sequence_version}"
        length = self.length - (len(prefix) + len(suffix))
        return f"{prefix}{self.sequence_index:0{length}}{suffix}"


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
class WgsSequenceRange:
    """
    This represents one or more WgsSequenceIds. This can be created from a
    single sequence id like, "ALWZ04S3033285", or a range of them like:
    "ALWZ04S0000001-ALWZ04S3033285". This requires that all sequences in the
    set have the same prefix. Additionally, the kind of sequence must be
    provided.
    """

    prefix: WgsPrefix
    kind: WgsSequenceKind
    start: WgsSequenceId
    stop: WgsSequenceId

    @classmethod
    def from_endpoint(cls, kind: WgsSequenceKind, endpoint: str) -> WgsSequenceRange:
        """
        Build a sequence set that represents a single sequence.
        """
        parsed = WgsSequenceId.build(endpoint)
        return cls(
            prefix=parsed.prefix,
            kind=kind,
            start=parsed,
            stop=parsed,
        )

    @classmethod
    def from_range(cls, kind: WgsSequenceKind, range: str) -> WgsSequenceRange:
        """
        Build a sequence set that represents a range of sequences.
        """
        raw_start, raw_stop = range.split("-", 1)
        endpoint1 = WgsSequenceId.build(raw_start)
        endpoint2 = WgsSequenceId.build(raw_stop)

        assert (
            endpoint1.prefix == endpoint2.prefix
        ), f"Cannot create range across prefixes {range}"
        assert len(raw_start) == len(
            raw_stop
        ), f"Cannot create range with bad lengths {range}"

        return cls(
            prefix=endpoint1.prefix,
            kind=kind,
            start=endpoint1,
            stop=endpoint2,
        )

    @classmethod
    def build(cls, kind: WgsSequenceKind, range: str) -> WgsSequenceRange:
        if "-" not in range:
            return cls.from_endpoint(kind, range)
        return cls.from_range(kind, range)

    def is_single_range(self) -> bool:
        return self.start == self.stop

    def to_wgs_string(self) -> str:
        if self.start == self.stop:
            return self.start.to_wgs_string()
        return f"{self.start.to_wgs_string()}-{self.stop.to_wgs_string()}"

    def ids(self) -> ty.Iterable[WgsSequenceId]:
        start = self.start.sequence_index
        stop = self.stop.sequence_index + 1
        for index in range(start, stop):
            yield WgsSequenceId(
                prefix=self.start.prefix,
                sequence_index=index,
                sequence_version=None,
                length=self.start.length,
            )

    def includes(
        self, accession: ty.Union[str, WgsSequenceId], within_one_version=False
    ) -> bool:
        if isinstance(accession, str):
            accession = WgsSequenceId.build(accession)

        if not self.prefix.matches(
            accession.prefix, within_one_version=within_one_version
        ):
            return False

        return (
            self.start.sequence_index
            <= accession.sequence_index
            <= self.stop.sequence_index
        )


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


@frozen
class WgsSummary:
    prefix: WgsPrefix
    contigs: ty.List[ContigInfo] = field(hash=False)
    sequences: ty.List[WgsSequenceRange] = field(hash=False)
    ncbi_ids: ty.List[str] = field(hash=False)

    def largest_ids(self) -> ty.Iterable[str]:
        if self.contigs:
            for contig in self.contigs:
                yield from contig.ids()
        elif self.sequences:
            mapped: ty.Dict[
                WgsSequenceKind, ty.List[WgsSequenceRange]
            ] = coll.defaultdict(list)
            for sequence in self.sequences:
                mapped[sequence.kind].append(sequence)

            for kind in WgsSequenceKind:
                seen = False
                for wgs_sequence in mapped.get(kind, []):
                    for seqid in wgs_sequence.ids():
                        yield seqid.to_wgs_string()
                    seen = True
                if seen:
                    break
            else:
                raise ValueError("Failed to produce and wgs sequences")

        elif self.ncbi_ids:
            yield from self.ncbi_ids
        else:
            raise ValueError("Somehow failed to have any ids in a wgs set")

    def includes_sequence(
        self, endpoint: WgsSequenceId, within_one_version=False
    ) -> bool:
        if self.prefix.matches(endpoint.prefix, within_one_version=within_one_version):
            return False

        for sequence in self.sequences:
            if sequence.includes(endpoint, within_one_version=within_one_version):
                return True

        return False


def looks_like_wgs_accession(raw: str) -> bool:
    try:
        WgsSequenceId.build(raw)
        return True
    except InvalidWgsAccession:
        return False


def resolve_ena_wgs(
    accession: str,
) -> ty.Tuple[ty.List[ContigInfo], ty.List[WgsSequenceRange]]:
    LOGGER.info("Fetching EMBL formatted file for %s", accession)
    url = ena.ENA_EMBL_URL.format(accession=accession)
    contigs: ty.List[ContigInfo] = []
    info: ty.List[WgsSequenceRange] = []
    try:
        with wget.wget(url) as handle:
            for line in handle:
                prefix = line[0:3]
                try:
                    kind = cattrs.structure(prefix, WgsSequenceKind)
                except:
                    continue
                try:
                    info.append(WgsSequenceRange.build(kind, line[3:].strip()))
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
    return WgsSummary(
        prefix=ena_info[0].prefix,
        contigs=contigs,
        sequences=ena_info,
        ncbi_ids=(ncbi_ids or []),
    )
