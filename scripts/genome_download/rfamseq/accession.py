# -*- coding: utf-8 -*-

"""
Copyright [2009-${2022}] EMBL-European Bioinformatics Institute
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
import re
import typing as ty

import cattrs
from attrs import field, frozen

WGS_SET_PATTERN = r""
WGS_SEQUENCE_PATTERN = r""
ACCESSION_PATTERN = re.compile(r"[A-Z]+.+")


@enum.unique
class AccessionKind(enum.Enum):
    GENERIC = "Generic Sequence"
    REFSEQ_GENOME = "RefSeq Genome"
    GENBANK_GENOME = "Standard Genome"
    WGS_SET = "WGS Set"
    WGS_SEQUENCE = "WGS Seqeunce"

    @classmethod
    def guess(cls, raw: str) -> None | AccessionKind:
        """Guess the kind of Accession a sequence is. This will check if the
        string matches several patterns to guess what kind of accession it will
        be. This may not be perfect but it should work for most cases. If no
        pattern matches it will return None. It is an error to provide this an
        empty string.

        This can handle generic accession with and without versions.

        >>> AccessionKind.guess("D00050")
        AccessionKind.GENERIC
        >>> AccessionKind.guess("D00050.1")
        AccessionKind.GENERIC
        >>> AccessionKind.guess("CY130004")
        AccessionKind.GENERIC
        >>> AccessionKind.guess("NC_004448.1")
        AccessionKind.GENERIC
        >>> AccessionKind.guess("NC_004448.10")
        AccessionKind.GENERIC

        This can handle RefSeq genomes with and without versions.

        >>> AccessionKind.guess("GCF_000512975.1")
        AccessionKind.REFSEQ_GENOME
        >>> AccessionKind.guess("GCF_000512975")
        AccessionKind.REFSEQ_GENOME

        This handles generic GCA accesssions with and without versions.

        >>> AccessionKind.guess("GCA_000769135.1")
        AccessionKind.STANDARD_GENOME
        >>> AccessionKind.guess("GCA_000769135")
        AccessionKind.GENBANK_GENOME

        This can handle WGS Set ids. In both short and long forms. It can also
        handle the ones with the 'S' suffix.

        >>> AccessionKind.guess("JABDTM02")
        AccessionKind.WGS_SET
        >>> AccessionKind.guess("ACTP02")
        AccessionKind.WGS_SET
        >>> AccessionKind.guess("ALWZ04S")
        AccessionKind.WGS_SET
        >>> AccessionKind.guess("ALWZ04S0000000")
        AccessionKind.WGS_SET

        This can handle WGS Sequence ids with or without a version suffix. As
        well as those of varying lengths.

        >>> AccessionKind.guess("ALWZ04S3033285")
        AccessionKind.WGS_SEQUENCE
        >>> AccessionKind.guess("JABDTM020000010")
        AccessionKind.WGS_SEQUENCE
        >>> AccessionKind.guess("ACTP02000001")
        AccessionKind.WGS_SEQUENCE
        >>> AccessionKind.guess("BAMX01000058.1")
        AccessionKind.WGS_SEQUENCE

        Given random strings it gives None.

        >>> AccessionKind.guess("1")
        None
        >>> AccessionKind.guess("a")
        None
        """

        if not raw:
            raise ValueError("Cannot guess accession kind of empty string")

        if raw.startswith("GCA_"):
            return cls.GENBANK_GENOME
        if raw.startswith("GCF_"):
            return cls.REFSEQ_GENOME
        if re.match(WGS_SET_PATTERN, raw):
            return cls.WGS_SET
        if re.match(WGS_SEQUENCE_PATTERN, raw):
            return cls.WGS_SEQUENCE
        if re.match(ACCESSION_PATTERN, raw):
            return cls.GENERIC
        return None


@frozen(hash=True)
class Accession:
    """
    This is meant to represent an 'accession' that sequences are identified by.
    Accessions are generally a '<string>.<version>'. Depending upon where the
    sequence came from and in what context it may or may not have a version.
    This class is used instead of treating accession as a raw string as
    sometimes we have versioned ones (with the '.<version>') and sometimes we do
    not. However, we always want to compare them ignoring the version. However,
    we still need to track the version for some usages.

    :accession: The accession without the version, if it exists
    :version: The version for this accession
    """

    accession: str
    version: ty.Optional[str]
    aliases: ty.Tuple[Accession] = field(factory=tuple, converter=tuple)

    @classmethod
    def build(cls, raw: str, aliases=None) -> Accession:
        """Build a new Accession from the given raw string."""
        if not raw:
            raise ValueError("Cannot parse an empty accession")
        aliases = aliases or tuple()
        parts = raw.split(".", 1)
        if len(parts) == 2:
            return cls(parts[0], parts[1], aliases=aliases)
        return cls(raw, None, aliases=aliases)

    def alias(self, alias: Accession) -> Accession:
        """Create a new one with which adds the given Accession as an alias."""

        aliases = list(self.aliases)
        aliases.append(alias)
        return Accession(
            accession=self.accession, version=self.version, aliases=tuple(aliases)
        )

    def matches(self, accession: Accession) -> bool:
        """Check if the given accession matches this one, ignoring the version if
        any of both accessions. If a versioned check needs to be done then just
        using `==` between two Accessions will work. Or converting this to a
        string and comparing those.
        """
        return (
            accession.accession == self.accession
            or any(a.matches(accession) for a in self.aliases)
            or any(a.matches(self) for a in accession.aliases)
        )

    def is_newer(self, other: Accession) -> bool:
        """Check if this accession is newer than the given one. This will
        return False if the other accession does not have same base or if the
        other is newer.

        >>> Accession.build("a.2").is_newer(Accession.build("a.1"))
        True
        >>> Accession.build("a.1").is_newer(Accession.build("a.2"))
        False
        >>> Accession.build("a").is_newer(Accession.build("a.1"))
        False
        >>> Accession.build("a").is_newer(Accession.build("a"))
        False
        >>> Accession.build("a.1").is_newer(Accession.build("a.1"))
        False
        >>> Accession.build("b.1").is_newer(Accession.build("a.1"))
        False
        >>> Accession.build("b.1").is_newer(Accession.build("a.2"))
        False
        >>> Accession.build("b.2").is_newer(Accession.build("a.1"))
        False
        """

        if other.accession != self.accession:
            return False
        return self.accession > other.accession

    def strip_version(self) -> Accession:
        """Generate a versionless accession.

        >>> Accession.build("a.1").strip_version()
        Accession.build("a")
        >>> Accession.build("a").strip_version()
        Accession.build("a")
        """

        return Accession(
            accession=self.accession,
            version=None,
            aliases=tuple(a.strip_version() for a in self.aliases),
        )

    def unaliased(self) -> Accession:
        return Accession(
            accession=self.accession,
            version=self.version,
            aliases=tuple(),
        )

    def is_versioned(self) -> bool:
        """Check if this accession is versioned.

        >>> Accession.build("a").is_versioned()
        False
        >>> Accession.build("a.2").is_versioned()
        True
        """
        return self.version is not None

    def increment(self, size=1, default_version=1) -> Accession:
        """Will produce a new Accession with the version incremented. If this
        Accession has no version, it will be assumed to be `default_version`
        and thus incremented to 2. This will remove all aliases from the
        Accession.

        >>> Accession.build("a.1").increment()
        Accession.build("a.2")
        >>> Accession.build("a").increment(size=2)
        Accession.build("a.3")
        >>> Accession.build("a").increment(size=2, default_version=2)
        Accession.build("a.4")
        """

        assert size > 0, "Must give positive size"
        assert default_version > 0, "Must give positive default version"
        version = default_version
        if self.version:
            version = int(self.version)

        version = version + size
        return Accession(
            accession=self.accession,
            version=str(version),
            aliases=tuple(),
        )

    def is_genomic(self):
        """Check if this is a GCA or GCF accession, which represent a genome.
        Other generic sequences may or may not represent a genome, but this
        does not check that.

        >>> Accession.build("GCA_00001.1").is_genomic()
        True
        >>> Accession.build("GCA_00001").is_genomic()
        True
        >>> Accession.build("NM00001.1").is_genomic()
        False
        >>> Accession.build("NM00001").is_genomic()
        False
        """

        return self.accession.startswith("GCA_") or self.accession.startswith("GCF_")

    def __str__(self) -> str:
        if self.version:
            return f"{self.accession}.{self.version}"
        return self.accession


def structure_accession(raw: ty.Any, _) -> Accession:
    if isinstance(raw, str):
        return Accession.build(raw)
    if isinstance(raw, dict):
        return Accession(**raw)
    raise ValueError(f"Cannot structure {raw} to Accession")


cattrs.register_structure_hook(Accession, structure_accession)
