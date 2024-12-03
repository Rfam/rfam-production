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

from __future__ import annotations

import logging

from attrs import frozen

LOGGER = logging.getLogger(__name__)


def none_if(value: str, flag: str) -> None | str:
    if not value or value == flag:
        return None
    return value


@frozen
class MAGInfo:
    accession: str
    gca_accession: None | str
    n50: int
    completeness: float
    contamination: float
    taxid: None | int
    species_level: bool

    @classmethod
    def from_mgnify(cls, data) -> MAGInfo:
        species_level = data["Species_level"].lower() == "true"
        taxid = None
        if data["Taxid"]:
            taxid = int(data["Taxid"])

        return cls(
            accession=data["Genome"],
            gca_accession=none_if(data["GCA_accession"], "N/A"),
            n50=int(data["N50"]),
            completeness=float(data["Completeness"]),
            contamination=float(data["Contamination"]),
            taxid=taxid,
            species_level=species_level,
        )

    @property
    def genome_url(self) -> str:
        return f"https://www.ebi.ac.uk/metagenomics/api/v1/genomes/{self.accession}/downloads/{self.accession}.fna"


@frozen
class Selector:
    """A class to select MAGs based upon some simple quality criteria

    :param require_species bool: Limit the MAGs to species level.
    :param min_completeness float: Require all MAGs to have at least this much
    completness
    :param max_contamination float: Require all MAGs to have less than this
    much contamination.
    """

    require_species: bool
    min_completeness: float
    max_contamination: float

    def accepts(self, mag: MAGInfo) -> bool:
        """Check if a given MAG passes the cutoffs this selector requires."""

        if self.require_species and not mag.species_level:
            LOGGER.debug("Reject %s, not species", mag.accession)
            return False

        if mag.completeness < self.min_completeness:
            LOGGER.debug(
                "Reject %s, incomplete %f < %f",
                mag.accession,
                mag.completeness,
                self.min_completeness,
            )
            return False

        if mag.contamination >= self.max_contamination:
            LOGGER.debug(
                "Reject %s, too containinated, %f > %f",
                mag.accession,
                mag.contamination,
                self.max_contamination,
            )
            return False

        return True
