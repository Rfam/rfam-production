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


import itertools as it
import operator as op
import typing as ty

import cattrs
from attrs import frozen

from rfam_export.context import Context

GENOME = """
SELECT
    upid 'accession',
    genome.ncbi_id,
    scientific_name,
    taxonomy.tax_string
FROM genome
JOIN taxonomy USING (ncbi_id)
WHERE
    upid = %(genome_accession)s
    AND EXISTS(SELECT 1 FROM unique_genseq WHERE unique_genseq.version = %(rfamseq_version)s AND unique_genseq.upid = genome.upid)
"""


CHROMOSOMES = """
SELECT
	rfamseq_acc,
    rfamseq.description
FROM genome
JOIN taxonomy USING (ncbi_id)
JOIN unique_genseq USING (upid)
JOIN rfamseq USING (rfamseq_acc)
WHERE
    upid = %(genome_accession)s
    AND unique_genseq.version = %(rfamseq_version)s
"""


ALL_GENOMES = """
SELECT
    upid 'accession',
    genome.ncbi_id,
    scientific_name,
    taxonomy.tax_string,
	rfamseq_acc,
    rfamseq.description
FROM genome
JOIN taxonomy USING (ncbi_id)
JOIN unique_genseq USING (upid)
JOIN rfamseq USING (rfamseq_acc)
WHERE
    unique_genseq.version = %(rfamseq_version)s
order by upid
"""


ALL_GENOME_ACCESSIONS = """
SELECT
    upid 'accession'
FROM genome
WHERE
    exists(select 1 from unique_genseq where unique_genseq.upid = genome.upid and unique_genseq.version = %(rfamseq_version)s)
order by upid
"""


@frozen
class Chromosome:
    rfamseq_acc: str
    description: str


@frozen
class Genome:
    accession: str
    ncbi_id: int
    scientific_name: str
    tax_string: str
    chromosomes: ty.Dict[str, Chromosome]

    def chromosome(self, accession: str) -> Chromosome:
        return self.chromosomes[accession]


def genome_info(cursor, context: Context, accession: str) -> Genome:
    params = {"genome_accession": accession, "rfamseq_version": context.rfamseq_version}
    cursor.execute(GENOME, params)
    genome = cursor.fetchone()
    cursor.execute(CHROMOSOMES, params)
    chromosomes = {}
    for chromosome in cursor.fetchall():
        chromosomes[chromosome["rfamseq_acc"]] = chromosome
    genome["chromosomes"] = chromosomes
    return cattrs.structure(genome, Genome)


def all_genome_accessions(cursor, context: Context) -> ty.Iterable[str]:
    cursor.execute(ALL_GENOME_ACCESSIONS, {"rfamseq_version": context.rfamseq_version})
    for row in cursor:
        yield row["accession"]


def all_genomes(cursor, context: Context) -> ty.Iterable[Genome]:
    cursor.execute(ALL_GENOMES, {"rfamseq_version": context.rfamseq_version})
    grouped = it.groupby(cursor, key=op.itemgetter("upid"))
    for upi, rows in grouped:
        chromosomes = {}
        ncbi_id = None
        scientific_name = None
        tax_string = None
        for row in rows:
            ncbi_id = row["ncbi_id"]
            scientific_name = row["scientific_name"]
            tax_string = row["tax_string"]
            chromosome = Chromosome(
                rfamseq_acc=row["rfamseq_acc"],
                description=row["description"],
            )
            chromosomes[chromosome.rfamseq_acc] = chromosome
        if scientific_name is None:
            raise ValueError("Never found a scientific name")
        if ncbi_id is None:
            raise ValueError("Never found an ncbi_id")
        if tax_string is None:
            raise ValueError("Never found a tax_string")
        yield Genome(
            accession=upi,
            ncbi_id=int(ncbi_id),
            scientific_name=scientific_name,
            tax_string=tax_string,
            chromosomes=chromosomes,
        )
