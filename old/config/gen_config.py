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

from config import rfam_local as cfl

# -----------------------------------CONSTANTS---------------------------------

# ENA file formats
FORMATS = {"xml": ".xml", "fasta": ".fa"}

# Taxonomy Domains
DOMAINS = ["eukaryota", "archaea", "viruses", "bacteria"]

# LSF JOB constraints & params
MEM = 5000
TMP_MEM = 1000
USER_EMAIL = cfl.USER_EMAIL
MAX_ALLOWED_FILES = 2000
CPUS = 5

# LSF GROUP constraints
LSF_GROUPS_CMD = cfl.LSF_GROUPS_CMD
LSF_GEN_GROUP = cfl.LSF_GEN_GROUP
JOB_LIMIT = 5
LSF_RFAM_BIN = cfl.BIN_LOC

# ENA GCA report file label
GCA_REP_LBL = "Sequence_Report"  # GCA ASSEMPLY REPORT FILE XML LABEL
GCA_REG_LBL = "Regions"  # GCA ASSEMPLY REGIONS FILE XML LABEL

# Genome Search variables
SPLIT_SIZE = 5427083
SRCH_MEM = 36000
SCAN_MEM = 36000
RFAMSEQ_SIZE = 451031.997884  # size of rfamseq13 in Mb
CM_NO = 2588  # number of cms in Rfam.cm file
CPU_NO = 5
SRCH_GROUP = "/rfam_search"

# -----------------------------------URLs--------------------------------------

# This returns the list of reference proteomes from Uniprot
REF_PROT_LIST_URL = (
    "http://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes&format=list"
)

# Retrieve a list of ref. Proteoms including GCA accessions and lineage
REF_PROT_REST_URL = (
    "http://www.uniprot.org/proteomes/?fil=reference:yes&"
    "force=no&format=tab&columns=id,assembly,lineage"
)

# Retrieve a sorted list of ref. Proteoms including GCA accessions,
# species name, id and lineage
"""
REF_PROT_REST_URL = ("http://www.uniprot.org/proteomes/?sort=&desc=&"
                     "compress=no&query=&fil=reference:yes&force=no&"
                     "format=tab&columns=id,name,organism-id,assembly,lineage")
"""

# non sorted list
"""
REF_PROT_REST_URL = ("http://www.uniprot.org/proteomes/?fil=reference:yes&"
                     "force=no&format=tab&columns=id,name,organism-id,"
                     "assembly,lineage")
"""

# Retrieve the proteome rdf file
PROTEOME_URL = "http://www.uniprot.org/proteomes/%s.rdf"

# Proteome xml file
PROTEOME_XML_URL = "http://www.uniprot.org/proteomes/%s.xml"

NCBI_SEQ_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=fasta&retmode=text"'

# Retrieve the genome's xml file
ENA_XML_URL = "http://www.ebi.ac.uk/ena/data/view/%s&display=xml"

# ENA url for file download
ENA_DATA_URL_GZIP = "http://www.ebi.ac.uk/ena/data/view/%s&display=%s&download=gzip"
ENA_DATA_URL = "http://www.ebi.ac.uk/ena/data/view/%s&display=%s"
# ENA url for assembly data retrieval via taxon id
ENA_TAX_URL = 'http://www.ebi.ac.uk/ena/data/warehouse/search?query="tax_eq(%s)"&result=assembly&display=xml'

# path to genome downloader executable
# GEN_DWLD_EXEC = cfl.GEN_DWLD_EXEC

# path to public WGS sets
# ENA_FTP_WGS_PUB = cfl.ENA_FTP_WGS_PUB

# path to supressed WGS sets
# ENA_FTP_WGS_SUP = cfl.ENA_FTP_WGS_SUP

# path to ena assembly sequence report files
# ENA_GCA_SEQ_REPORT = cfl.ENA_GCA_SEQ_REPORT

# -----------------------------------MODELS------------------------------------

GENOME_MODEL = "RfamLive.Genome"  # table names RfamLive.Genome
GENSEQ_MODEL = "RfamLive.Genseq"  # table names RfamLive.Genseq
