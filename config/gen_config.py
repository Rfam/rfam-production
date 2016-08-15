'''
Created on 4 Jul 2016

@author: ikalvari

Description: Genome's Project - Configuration File
Purpose: Contains all constants used in the genome's project. Uniprot and ENA
URLS can be accessed though the corresponding REST API webpages.

'''

from config import rfam_local as cfl

# -----------------------------------CONSTANTS---------------------------------

# ENA file formats
FORMATS = {"xml": ".xml", "fasta": ".fa"}

# Taxonomy Domains
DOMAINS = ['eukaryota', 'archaea', 'viruses', 'bacteria']

# Rfam's LSF GPFS location
RFAM_GPFS_LOC = cfl.RFAM_GPFS_LOC
LOC_PATH = cfl.LOC_PATH

# LSF JOB constraints & params
MEM = 5000
TMP_MEM = 1000
USER_EMAIL = cfl.USER_EMAIL

# LSF GROUP constraints
LSF_GROUPS_CMD = cfl.LSF_GROUPS_CMD
LSF_GEN_GROUP = cfl.LSF_GEN_GROUP
JOB_LIMIT = 5

# ENA GCA report file label
GCA_REP_LBL = "Sequence_Report"  # GCA ASSEMPLY REPORT FILE XML LABEL
GCA_REG_LBL = "Regions"  # GCA ASSEMPLY REGIONS FILE XML LABEL

# -----------------------------------URLs--------------------------------------

# This returns the list of reference proteomes from Uniprot
REF_PROT_LIST_URL = "http://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes&format=list"

# Retrieve a list of ref. Proteoms including GCA accessions and lineage
REF_PROT_REST_URL = ("http://www.uniprot.org/proteomes/?fil=reference:yes&"
                     "force=no&format=tab&columns=id,assembly,lineage")

# Retrieve a sorted list of ref. Proteoms including GCA accessions,
# species name, id and lineage
'''
REF_PROT_REST_URL = ("http://www.uniprot.org/proteomes/?sort=&desc=&"
                     "compress=no&query=&fil=reference:yes&force=no&"
                     "format=tab&columns=id,name,organism-id,assembly,lineage")
'''

# non sorted list
'''
REF_PROT_REST_URL = ("http://www.uniprot.org/proteomes/?fil=reference:yes&"
                     "force=no&format=tab&columns=id,name,organism-id,"
                     "assembly,lineage")
'''

# Retrieve the proteome rdf file
PROTEOME_URL = "http://www.uniprot.org/proteomes/%s.rdf"

# Retrieve the genome's xml file
ENA_XML_URL = "http://www.ebi.ac.uk/ena/data/view/%s&display=xml"

# ENA url for file download
ENA_DATA_URL = "http://www.ebi.ac.uk/ena/data/view/%s&display=%s&download=gzip"

# ENA url for assembly data retrieval via taxon id
ENA_TAX_URL = "http://www.ebi.ac.uk/ena/data/warehouse/search?query=\"tax_eq(%s)\"&result=assembly&display=xml"

# path to genome downloader executable
GEN_DWLD_EXEC = cfl.GEN_DWLD_EXEC

# -----------------------------------MODELS------------------------------------

GENOME_MODEL = "RfamLive.Genome"  # table names RfamLive.Genome
GENSEQ_MODEL = "RfamLive.Genseq"  # table names RfamLive.Genseq
