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

import re
import typing as ty
import xml.etree.ElementTree as ET
from enum import Enum, unique
from pathlib import Path

import cattrs
from attr import frozen
from loguru import logger

from rfam_export import xml
from rfam_export.context import Context
from rfam_export.genomes import Genome

FULL_REGION_FIELDS = """
SELECT
    fr.rfamseq_acc,
    CONCAT(fr.rfamseq_acc, '/', fr.seq_start, ':', fr.seq_end) accession,
    CONCAT(tx.species, ' ', f.rfam_id) description,
    f.rfam_acc,
    fr.seq_start,
    fr.seq_end,
    rs.description as rfamseq_acc_description,
    fr.cm_start,
    fr.cm_end,
    fr.evalue_score,
    fr.bit_score,
    fr.type as alignment_type,
    fr.truncated,
    fr.rfam_acc,
    f.rfam_id,
    f.type as rna_type,
    rs.ncbi_id as ncbi_id,
    tx.species as scientific_name,
    tx.tax_string
FROM full_region fr, family f, unique_genseq gs, rfamseq rs, taxonomy tx
WHERE
    fr.rfamseq_acc = gs.rfamseq_acc
    AND rs.ncbi_id = tx.ncbi_id
    AND gs.rfamseq_acc = rs.rfamseq_acc
    AND fr.rfam_acc = f.rfam_acc
    AND fr.is_significant = 1
    AND fr.type = 'full'
    AND gs.upid = %(genome_accession)s
    AND gs.version = %(rfamseq_version)s
order by rfamseq_acc, seq_start, seq_end
"""

FULL_REGION_SEEDS = """
SELECT
    fr.rfamseq_acc,
    CONCAT(fr.rfamseq_acc, '/', fr.seq_start, ':', fr.seq_end) accession,
    CONCAT(tx.species, ' ', f.rfam_id) description,
    fr.seq_start,
    fr.seq_end,
    fr.cm_start,
    fr.cm_end,
    fr.evalue_score,
    fr.bit_score,
    fr.type as alignment_type,
    fr.truncated,
    fr.rfam_acc,
    f.rfam_id,
    f.type as rna_type,
    rs.description as rfamseq_acc_description,
    rs.ncbi_id as ncbi_id,
    tx.species as scientific_name,
    tx.tax_string
FROM full_region fr, family f, unique_genseq gs, rfamseq rs, taxonomy tx
WHERE
    fr.rfamseq_acc=gs.rfamseq_acc
    AND rs.ncbi_id=tx.ncbi_id
    AND gs.rfamseq_acc=rs.rfamseq_acc
    AND fr.rfam_acc=f.rfam_acc
    AND fr.is_significant=1
    AND fr.type = 'full'
    AND gs.upid = %(genome_accession)s
    AND gs.version = %(rfamseq_version)s
order by rfamseq_acc, seq_start, seq_end
"""

RNACENTRAL_MAPPING = """
SELECT
    rm.rfamseq_acc,
    rm.seq_start,
    rm.seq_end,
    rm.rnacentral_id
from rnacentral_matches rm
join full_region fr on
    rm.rfamseq_acc = fr.rfamseq_acc
    AND fr.rfamseq_acc = rm.rfamseq_acc
    AND fr.seq_start = rm.seq_start
    AND fr.seq_end = rm.seq_end
    AND fr.is_significant = 1
join unique_genseq gs on gs.rfamseq_acc = fr.rfamseq_acc
WHERE
    rm.rnacentral_id IS NOT NULL
    AND gs.upid = %(genome_accession)s
    AND gs.version = %(rfamseq_version)s
"""


class NoHits(Exception):
    """Raised when we do not find any hits for a given genome"""


@unique
class AlignmentType(Enum):
    SEED = "seed"
    FULL = "full"


@frozen
class Hit:
    accession: str = xml.element(name="name")
    description: str = xml.element()
    entry_type: xml.EntryType = xml.additional_field()
    rfamseq_acc: str = xml.additional_field()
    rfamseq_acc_description: str = xml.additional_field()
    seq_start: int = xml.additional_field()
    seq_end: int = xml.additional_field()
    cm_start: int = xml.additional_field()
    cm_end: int = xml.additional_field()
    evalue_score: float = xml.additional_field()
    bit_score: float = xml.additional_field()
    alignment_type: AlignmentType = xml.additional_field()
    truncated: int = xml.additional_field()
    tax_string: str = xml.additional_field()
    rna_type: ty.List[str] = xml.additional_field()
    scientific_name: str = xml.additional_field()
    ncbi_id: int = xml.cross_reference("ncbi_taxonomy_id")
    ena_id: None | str = xml.cross_reference("ENA")
    rfam_acc: str = xml.cross_reference("RFAM")
    rnacentral_id: str = xml.cross_reference("RNACENTRAL")
    proteome_id: None | str = xml.cross_reference("Uniprot")

    @property
    def entry_id(self) -> str:
        return self.accession.replace(".", "_").replace("/", "_")


def rnacentral_mapping(
    cursor, context: Context, accession: str
) -> ty.Iterator[ty.Tuple[str, ty.Dict[str, str]]]:
    cursor.execute(
        RNACENTRAL_MAPPING,
        {"genome_accession": accession, "rfamseq_version": context.rfamseq_version},
    )
    for row in cursor:
        name = "%s/%s:%s" % (row["rfamseq_acc"], row["seq_start"], row["seq_end"])
        yield (name, (row["rnacentral_id"]))


def fetch_entries(
    cursor, query: str, context: Context, genome: Genome
) -> ty.Iterable[Hit]:
    """Extract all FULL hits from the database for the given `Genome`. This
    will execute a query, iterate the results and produce an `Iterable` of
    `FullHit` objects.
    """

    cursor.execute(
        query,
        {
            "genome_accession": genome.accession,
            "rfamseq_version": context.rfamseq_version,
        },
    )
    with (
        context.cached("rnacentral_mapping") as rnacentral_mapping,
        context.cached("families") as families,
    ):

        last = None
        for row in cursor:
            family = families[row["rfam_acc"]]

            rnacentral_id = row["rfamseq_acc"]
            if row["accession"] in rnacentral_mapping:
                urs = rnacentral_mapping[row["accession"]]
                taxid = row["ncbi_id"]
                rnacentral_id = f"{urs}_{taxid}"

            proteome_id = None
            if genome.accession.startswith("UP"):
                proteome_id = genome.accession

            ena_id = None
            if "." in row["rfamseq_acc"]:
                ena_id = row["rfamseq_acc"].split(".", 1)[0]
            elif row["rfamseq_acc"].startswith("URS"):
                ena_id = row["rfamseq_acc"]

            chromosome = genome.chromosome(row["rfamseq_acc"])
            update = {
                "entry_type": xml.EntryType.SEQUENCE,
                "rna_type": re.sub(r";$", "", family.rna_type).split("; "),
                "rfam_id": family.id,
                "rfamseq_acc_description": chromosome.description,
                "rnacentral_id": rnacentral_id,
                "ena_id": ena_id,
                "proteome_id": proteome_id,
            }

            row.update(update)
            current = cattrs.structure(row, Hit)
            if current == last:
                logger.trace("Found duplicate entry: {}", current.accession)
                continue
            yield current
            last = current

        if last is None:
            raise NoHits(genome)


def entries(cursor, context: Context, genome: Genome) -> ty.Iterable[ET.Element]:
    context.cache(
        "rnacentral_mapping", rnacentral_mapping(cursor, context, genome.accession)
    )
    for hit in fetch_entries(cursor, FULL_REGION_FIELDS, context, genome):
        yield xml.as_xml(hit)
    for hit in fetch_entries(cursor, FULL_REGION_SEEDS, context, genome):
        yield xml.as_xml(hit)


def export_region(cursor, context: Context, genome: Genome, output: ty.TextIO):
    elements = entries(cursor, context, genome)
    xml.write_file(context, elements, output)
