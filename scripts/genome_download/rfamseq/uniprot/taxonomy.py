# -*- coding: utf-8 -*-

"""
Copyright [2009-${2024}] EMBL-European Bioinformatics Institute
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

import typing as ty
from functools import lru_cache

import requests
from attrs import frozen


@frozen
class LineageInfo:
    ncbi_id: int
    species: str
    common_name: ty.Optional[str]
    tax_string: str


@frozen
class LineageEntry:
    taxon_id: int
    scientific_name: str
    rank: str
    hidden: bool = False


@lru_cache
def lineage_info(taxid: str) -> LineageInfo:
    response = requests.get(f"https://rest.uniprot.org/taxonomy/{taxid}.json")
    response.raise_for_status()
    data = response.json()
    parents = []
    should_hide = False
    for taxon in reversed(data["lineage"]):
        if taxon["taxonId"] == 2759:  # Hide in euks
            should_hide = True
        if should_hide and taxon.get("hidden", False):
            continue
        if taxon["scientificName"] == "cellular organisms":
            continue
        parents.append(taxon)
    tax_string = "; ".join(l["scientificName"] for l in parents)
    tax_string += "."

    common_name = data.get("commonName", None)
    if common_name == "":
        common_name = None

    return LineageInfo(
        ncbi_id=data["taxonId"],
        species=data["scientificName"],
        common_name=common_name,
        tax_string=tax_string,
    )


# def build_info(entries: ty.List[LineageEntry]) -> LineageInfo:
#     pass
