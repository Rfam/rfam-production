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

import cattrs
from attr import frozen


class MissingFamilies(Exception):
    """This is raised if no families are loaded for some reason"""


ALL_FAMILIES = """
SELECT
    rfam_acc accession,
    rfam_id id,
    description,
    type rna_type
FROM family
"""


@frozen
class Family:
    accession: str
    id: str
    description: str
    rna_type: str


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
