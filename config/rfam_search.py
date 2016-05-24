'''
Created on 13 May 2016

@author: ikalvari

Description: A library of scripts and constants to support rfam_xml_dumper

TO DO: - Set release version and release date automatically

'''
import datetime
# -----------------------------------------------------------------------------

# MAIN XML Fields
DB_NAME = "Rfam"  # DB name
DB_DESC = "A database for non-protein coding RNA families"  # DB description
DB_RELEASE = "12.1"  # release version
DB_REL_DATE = "26/04/2016"  # datetime.date.today()

# DELIMITERS
RNA_TYPE_DEL = ';'
AUTH_DEL = ','

# ENTRY TYPES
MOTIF = 'M'
CLAN = 'C'
FAMILY = 'F'

# RFAM SEARCH QUERIES
REL_FIELDS = ("SELECT rfam_release, rfam_release_date FROM version")

# -------------------------------RFAM ACCESSIONS--------------------------

# FETCHING

CLAN_ACC = ("SELECT clan_acc FROM clan")

MOTIF_ACC = ("SELECT motif_acc FROM motif")

FAM_ACC = ("SELECT rfam_acc FROM family")


# ---------------------------------RFAM FIELDS----------------------------

# Select all family related fields joining tables to fetch pubmed,go and so ids
FAM_FIELDS = ("SELECT f.rfam_acc as id, f.rfam_id as name, f.description,"
              "f.author, f.number_of_species as num_species,"
              "f.number_3d_structures as num_3d_structures, f.num_seed,"
              "f.num_full, f.type as rna_type, f.created, f.updated,"
              "group_concat(distinct concat(dl.db_id,\':\',dl.db_link)) as dbxrefs,"
              "group_concat(distinct concat(fl.pmid)) as pmids\n"
              "FROM family f JOIN full_region fr USING (rfam_acc)\n"
              "JOIN database_link dl using (rfam_acc)\n"
              "JOIN family_literature_reference fl USING (rfam_acc)\n"
              "WHERE f.rfam_acc = \'%s\' AND fr.is_significant=1\n"
              "AND (dl.db_id like \'GO\' OR dl.db_id like \'SO\')\n"
              "GROUP BY f.rfam_acc, f.rfam_id, f.description, f.author, f.number_of_species,"
              "f.number_3d_structures, f.num_seed, f.num_full, f.type, f.created, f.updated")

# Fetching clan fields from the db
CLAN_FIELDS = ("SELECT c.clan_acc as id, c.id as name, c.description, c.created,"
               "c.updated, c.author, count(*) as num_families\n"
               "FROM clan c, clan_membership cm\n"
               "WHERE c.clan_acc=cm.clan_acc\n"
               "AND c.clan_acc=\'%s\'")

# Fetching clan fields from the db

MOTIF_FIELDS = ("SELECT m.motif_acc as id, m.motif_id as name, m.description,"
                "m.created, m.updated, m.author\n"
                "FROM motif m\n"
                "WHERE m.motif_acc=\'%s\'")

# -----------------------------CROSS REFERENCES---------------------------

# FAMILIES
# Fetch pdb ids related to a family accession
PDB_IDs_QUERY = ("SELECT distinct pdb_id FROM pdb_full_region\n"
                 "WHERE rfam_acc=\'%s\'")

# Fetch ncbi ids related to a family accession
NCBI_IDs_QUERY = ("SELECT distinct tx.ncbi_id\n"
                  "FROM taxonomy tx, full_region fr, rfamseq rs\n"
                  "WHERE tx.ncbi_id=rs.ncbi_id\n"
                  "AND fr.rfamseq_acc=rs.rfamseq_acc\n"
                  "AND fr.is_significant=1\n"
                  "AND fr.rfam_acc=\'%s\'")

# Fetch clan based on rfam_acc
FAM_CLAN = ("SELECT clan_acc from clan_membership\n"
            "WHERE rfam_acc=\'%s\'")

# CLANS
# Fetch clan members
CLAN_FAMS = ("SELECT rfam_acc from clan_membership\n"
             "WHERE clan_acc=\'%s\'")

# MOTIFS
MOTIF_FAMS = ("SELECT distinct rfam_acc from motif_matches\n"
              "where motif_acc=\'%s\'")

# -------------------------ADDITIONAL FIELDS------------------------------

# count # families associated with clan
NUM_FAMS_CLAN = ("SELECT count(*) FROM clan_membership\n"
                 "WHERE clan_acc=\'%s\'")

# count # families associated with motif
NUM_FAMS_MOTIF = ("SELECT count(*) FROM motif_family_stats\n"
                  "WHERE motif_acc=\'%s\'")

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    pass
