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

import typing as ty

import cattrs
from attrs import field, frozen


@frozen
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
    aliases: ty.Tuple[Accession] = field(factory=tuple)

    @classmethod
    def build(cls, raw: str, aliases=None) -> Accession:
        """
        Build a new Accession from the given raw string.
        """
        aliases = aliases or tuple()
        parts = raw.split(".", 1)
        if len(parts) == 2:
            return cls(parts[0], parts[1], aliases=aliases)
        return cls(raw, None, aliases=aliases)

    def alias(self, alias: Accession) -> Accession:
        """
        Create a new one with which adds the given Accession as an alias.
        """
        aliases = list(self.aliases)
        aliases.append(alias)
        return Accession(
            accession=self.accession, version=self.version, aliases=tuple(aliases)
        )

    def matches(self, accession: Accession) -> bool:
        """
        Check if the given accession matches this one, ignoring the version if
        any of both accessions. If a versioned check needs to be done then just
        using `==` between two Accessions will work. Or converting this to a
        string and comparing those.
        """
        return accession == self.accession or any(
            a.matches(accession) for a in self.aliases
        )

    def versioned(self) -> bool:
        """
        Check if this accession is versioned.
        """
        return self.version is not None

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
