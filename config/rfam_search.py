"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

"""
Description: A library of scripts and constants to support rfam_xml_dumper

TO DO: - Set release version and release date automatically
"""

import datetime

# -----------------------------------------------------------------------------

# MAIN XML Fields
DB_NAME = "Rfam"  # DB name
DB_DESC = "A database for non-protein coding RNA families"  # DB description
DB_RELEASE = "15.1"  # release version
# DB_REL_DATE = datetime.date.today()  # datetime.date.today()

# DELIMITERS
RNA_TYPE_DEL = ';'
AUTH_DEL = ','

# ENTRY TYPES
MOTIF = 'M'
CLAN = 'C'
FAMILY = 'F'
GENOME = 'G'
MATCH = 'R'


# 9606 - human
# 10090 - mouse
# 7955 - zebrafish
# 3702 - Arabidopsis thaliana
# 6239 - Caenorhabditis elegans
# 7227 - Drosophila melanogaster
# 559292 - Saccharomyces cerevisiae S288c
# 4896 - Schizosaccharomyces pombe
# 511145 - Escherichia coli str. K-12 substr. MG1655

# MODEL ORGANISMS - popular species
POPULAR_SPECIES = (
    '9606', '10090', '7955', '3702', '6239', '7227', '559292', '4896', '511145')


# RFAM SEARCH QUERIES
REL_FIELDS = "SELECT rfam_release, rfam_release_date FROM version"

# -------------------------------RFAM ACCESSIONS--------------------------

# FETCHING

CLAN_ACC = "SELECT clan_acc FROM clan"

MOTIF_ACC = "SELECT motif_acc FROM motif"

FAM_ACC = "SELECT rfam_acc FROM family"

GENOME_ACC = "SELECT upid FROM genome WHERE num_families > 0"


# ---------------------------------RFAM FIELDS----------------------------

# Select all family related fields joining tables to fetch pubmed,go and so ids
FAM_FIELDS ="""SELECT f.rfam_acc as id, f.rfam_id as name, f.description,f.author, 
f.number_of_species as num_species,
f.number_3d_structures as num_3d_structures, f.num_seed,
f.num_full, f.type as rna_type, f.created, f.updated,
group_concat(distinct concat(dl.db_id,':',dl.db_link)) as dbxrefs,
group_concat(distinct concat(fl.pmid)) as pmids
FROM family f JOIN database_link dl using (rfam_acc)
JOIN family_literature_reference fl USING (rfam_acc)
WHERE f.rfam_acc = '%s'
AND (dl.db_id like 'GO' OR dl.db_id like 'SO')
GROUP BY f.rfam_acc, f.rfam_id, f.description, f.author, f.number_of_species,
f.number_3d_structures, f.num_seed, f.num_full, f.type, f.created, f.updated"""

# Fetching clan fields from the db
CLAN_FIELDS = """
              SELECT c.clan_acc as id, c.id as name, c.description, c.created,
              c.updated, c.author, count(*) as num_families
              FROM clan c, clan_membership cm
              WHERE c.clan_acc=cm.clan_acc
              AND c.clan_acc='%s'
              """

# Fetching clan fields from the db
MOTIF_FIELDS = """
               SELECT m.motif_acc as id, m.motif_id as name, m.description,
               m.created, m.updated, m.author
               FROM motif m
               WHERE m.motif_acc='%s'
               """

# Fetching genome fields from the db
GENOME_FIELDS = """
                SELECT g.upid as id, g.scientific_name as name, g.assembly_acc, g.description,
                g.total_length, tx.tax_string, g.ncbi_id, g.num_rfam_regions,g.num_families,
                g.common_name, g.assembly_name, g.assembly_level, g.created, g.updated
                FROM genome g, taxonomy tx
                WHERE g.ncbi_id=tx.ncbi_id
                AND g.upid='%s'
                """

FULL_REGION_FIELDS = """
    SELECT
    fr.rfamseq_acc, fr.seq_start, fr.seq_end, fr.cm_start, fr.cm_end, fr.evalue_score,
    fr.bit_score, fr.type as alignment_type, fr.truncated, fr.rfam_acc,
    f.rfam_id, f.type as rna_type, rs.description as rfamseq_acc_description,
    rs.ncbi_id as ncbi_id, tx.species as scientific_name, tx.tax_string
    FROM full_region fr, family f, genseq gs, rfamseq rs, taxonomy tx
    WHERE fr.rfamseq_acc=gs.rfamseq_acc
    AND rs.ncbi_id=tx.ncbi_id
    AND gs.rfamseq_acc=rs.rfamseq_acc
    AND fr.rfam_acc=f.rfam_acc
    AND fr.is_significant=1
    AND fr.type = 'full'
    AND gs.upid = '%s'
    AND gs.version = '15.0'
"""

FULL_REGION_SEEDS = """
    SELECT
    fr.rfamseq_acc, fr.seq_start, fr.seq_end, fr.cm_start, fr.cm_end, fr.evalue_score,
    fr.bit_score, fr.type as alignment_type, fr.truncated, fr.rfam_acc,
    f.rfam_id, f.type as rna_type, rs.description as rfamseq_acc_description, 
    tx.species as scientific_name, rs.ncbi_id as ncbi_id, tx.tax_string
    FROM full_region fr, family f, genseq gs, rfamseq rs, taxonomy tx
    WHERE fr.rfamseq_acc=rs.rfamseq_acc
    AND gs.rfamseq_acc=rs.rfamseq_acc
    AND rs.ncbi_id=tx.ncbi_id
    AND fr.rfam_acc=f.rfam_acc
    AND fr.is_significant=1
    AND fr.type='seed'
    AND gs.upid = '%s'
    AND gs.version='15.0'
"""

# -----------------------------CROSS REFERENCES---------------------------

# FAMILIES
# Fetch pdb ids related to a family accession
PDB_IDs_QUERY = """
                SELECT distinct pdb_id
                FROM pdb_full_region
                WHERE rfam_acc='%s'
                AND is_significant=1
                """

PSEUDOKNOTS_QUERY = """
    SELECT distinct pseudoknot_id
    FROM pseudoknots
    WHERE rfam_acc='%s'
    """

# Fetch ncbi ids related to a family accession
NCBI_IDs_QUERY = """
                 SELECT tx.ncbi_id, tx.tax_string
                 FROM taxonomy tx, full_region fr, rfamseq rs
                 WHERE tx.ncbi_id=rs.ncbi_id
                 AND fr.rfamseq_acc=rs.rfamseq_acc
                 AND fr.is_significant=1
                 AND fr.rfam_acc='%s'
                 GROUP BY tx.ncbi_id
                 """

# Fetch upids related to a family accession
FAMILY_UPIDS = """
               SELECT distinct upid
               FROM genseq gs, full_region fr
               WHERE gs.rfamseq_acc=fr.rfamseq_acc
               AND fr.rfam_acc='%s'
               AND fr.is_significant = 1
               AND gs.version='15.0'
               """


# Fetch clan based on rfam_acc
FAM_CLAN = """
           SELECT clan_acc
           FROM clan_membership
           WHERE rfam_acc='%s'
           """

# CLANS
# Fetch clan members
CLAN_FAMS = """
            SELECT rfam_acc
            from clan_membership
            WHERE clan_acc='%s'
            """

# MOTIFS
MOTIF_FAMS = """
             SELECT distinct rfam_acc
             FROM motif_matches
             WHERE motif_acc='%s'
             """

# GENOMES
GENOME_FAMS = """
              SELECT distinct rfam_acc
              FROM full_region fr, genseq gs
              WHERE fr.rfamseq_acc=gs.rfamseq_acc
              AND fr.is_significant = 1
              AND gs.upid='%s'
              AND gs.version='15.0'
              """


# -------------------------ADDITIONAL FIELDS------------------------------

# count # families associated with clan
NUM_FAMS_CLAN = """
                SELECT count(*) FROM clan_membership
                WHERE clan_acc='%s'
                """

# count # families associated with motif
NUM_FAMS_MOTIF = """
                 SELECT count(*) FROM motif_family_stats
                 WHERE motif_acc='%s'
                 """

COUNT_FULL_REGION = """
                    SELECT count(*)
                    FROM full_region fr, genseq gs
                    WHERE fr.rfamseq_acc = gs.rfamseq_acc
                    AND fr.is_significant = 1
                    AND fr.type='full'
                    AND gs.upid = '%s'
                    AND gs.version='15.0'
                    """

AU_ORCIDS = """
            SELECT orcid
            FROM author au, family_author fa
            WHERE au.author_id=fa.author_id
            AND rfam_acc='%s'
            AND orcid <> ''
            """

SEED_PK_WITH_COV = """
                SELECT count(*)
                FROM pseudoknot
                WHERE source='seed'
                AND covariation=1
                AND rfam_acc='%s'
                """

SEED_PK_NO_COV = """
    SELECT count(*)
    FROM pseudoknot
    WHERE source='seed'
    AND covariation=0
    AND rfam_acc='%s'
    """

RSCAPE_PK_WITH_COV = """
    SELECT count(*)
    FROM pseudoknot
    WHERE source='rscape'
    AND covariation=1
    AND rfam_acc='%s'
    """

RSCAPE_PK_NO_COV = """
    SELECT count(*)
    FROM pseudoknot
    WHERE source='rscape'
    AND covariation=0
    AND rfam_acc='%s'
    """
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    pass
