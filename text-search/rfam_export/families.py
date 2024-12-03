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

import enum
import typing as ty
from pathlib import Path

import cattrs
from attr import frozen
from loguru import logger
from rfam_export import xml
from rfam_export.context import Context


class MissingFamilies(Exception):
    """This is raised if no families are loaded for some reason"""


class UnknownFamily(Exception):
    """This is raise if given a that could not be loaded from the database."""


class MissingTaxonomy(Exception):
    """This is raised if a family has no taxonomic information."""


ALL_FAMILIES = """
SELECT
    rfam_acc accession,
    rfam_id id,
    description,
    type rna_type
FROM family
"""

FAMILY_INFO = """
SELECT
    f.rfam_acc as id,
    f.rfam_id as name,
    f.description,f.author,
    f.number_of_species as num_species,
    f.number_3d_structures as num_3d_structures,
    f.num_seed,
    f.num_full,
    f.type as rna_type,
    f.created,
    f.updated,
    group_concat(distinct concat(dl.db_id,':',dl.db_link)) as dbxrefs,
    group_concat(distinct concat(fl.pmid)) as pmids
FROM family f JOIN database_link dl using (rfam_acc)
JOIN family_literature_reference fl USING (rfam_acc)
WHERE
    f.rfam_acc = '%(accession)s'
    AND dl.db_id in ('GO', 'SO')
GROUP BY f.rfam_acc, f.rfam_id, f.description, f.author, f.number_of_species, f.number_3d_structures, f.num_seed, f.num_full, f.type, f.created, f.updated
"""

FAMILY_TAXONOMY = """
SELECT
    ncbi_id,
    lineage,
FROM taxonomy
JOIN full_regions
"""


@frozen
class Family:
    accession: str
    id: str
    description: str
    rna_type: str


@enum.unique
class PkSupport(enum.Enum):
    SEED_PK_COVARY = "seed with covariation support"
    SEED_PK_NO_COVARY = "seed no covariation support"
    RSCAPE_PK_COVARY = "predicted with covariation support"
    RSCAPE_PK_NO_COVARY = "predicted no covariation support"


@frozen
class FamilyElement:
    accession: str = xml.element(name="name")
    description: str = xml.element()
    ncbi_ids: ty.Set[int] = xml.cross_reference("ncbi_taxonomy_id")
    entry_type: xml.EntryType = xml.additional_field()
    author: ty.Set[str] = xml.additional_field()
    has_3d_structure: bool = xml.additional_field()
    has_pseudoknot: int = xml.additional_field()
    num_3d_structures: int = xml.additional_field()
    num_full: int = xml.additional_field()
    num_seed: int = xml.additional_field()
    num_species: int = xml.additional_field()
    pseudoknot_evidence: ty.Set[PkSupport] = xml.additional_field()
    rna_type: str = xml.additional_field()
    tax_string: ty.Set[str] = xml.additional_field()

    @property
    def entry_id(self) -> str:
        return self.accession


@frozen
class Exporter:
    context: Context
    cursor: pymysql.Cursor
    base_path: Path

    @classmethod
    def build(
        cls, context: Context, cursor: pymysql.Cursor, base_path: Path
    ) -> Exporter:
        return cls(context, cursor, base_path)

    def pk_evidence(self, family: Family) -> ty.List[PkSupport]:
        pass

    def family(self, family: Family) -> FamilyElement:
        self.cursor.execute(FAMILY_INFO, {"accession": family.accession})
        data = self.cursor.fetchone()
        if not data:
            raise UnknownFamily(family)

        pk_evidence = self.pk_evidence(family)

        data.update(
            {
                "pseudoknot_evidence": pk_evidence,
                "ncbi_ids": set(),
                "tax_string": set(),
                "entry_type": xml.EntryType.FAMILY,
                "has_3d_structure": bool(data.get("num_3d_structures", 0)),
                "has_pseudoknot": int(bool(pk_evidence)),
            }
        )
        self.cursor.execute(FAMILY_TAXONOMY, {"accession": family.accession})
        for result in self.cursor:
            data["ncbi_ids"].add(result["ncbi_id"])
            data["tax_string"].add(result["tax_string"])

        if not data["ncbi_ids"] or not data["tax_string"]:
            raise MissingTaxonomy(family)

        return cattrs.structure(data, FamilyElement)

    def export_family(self, family: Family):
        """Write the search XML for the given family to disk. This will create
        a new file in the self.base_path for the family. The file will only
        exist if export is successful.

        :family Family: The family to export.
        """
        path = self.base_path / f"{family.accession}.xml"
        logger.info("Exporting family {} to {}", family, path)
        try:
            entry = self.family(family)
            xml_element = xml.as_xml(entry)
            with path.open("w") as out:
                xml.write_file(self.context, [xml_element], out)
        except Exception as err:
            path.unlink(missing_ok=True)
            raise err


def families(cursor) -> ty.List[Family]:
    cursor.execute(ALL_FAMILIES)
    families = []
    for row in cursor:
        family = cattrs.structure(row, Family)
        families.append(family)
    if not families:
        raise MissingFamilies()
    return families


def family_mapping(cursor) -> ty.Iterable[ty.Tuple[str, Family]]:
    for family in families(cursor):
        yield (family.accession, family)
