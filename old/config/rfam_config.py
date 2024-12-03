"""
Created on 2 Mar 2016

@author: ikalvari

Description: This is a configuration file with all the paths on the cluster

Edited May 2022 By ecooke - codon migration

"""

from config import rfam_local as cfl

# ---------------------------------PATHS---------------------------------------
ESL_PATH = cfl.ESL_SFETCH
ESL_SEQSTAT = cfl.ESL_SEQSTAT
ESL_SSPLIT = cfl.ESL_SSPLIT
CMSEARCH = cfl.CMSEARCH
CMSCAN = cfl.CMSCAN
CMFILE = cfl.CMFILE
SCANNER = cfl.SCANNER
DWL_SCRIPT = cfl.DWL_SCRIPT
FA_GEN = cfl.FA_GEN
RFAMSEQ_PATH = cfl.RFAMSEQ_PATH
FAM_VIEW_PL = cfl.FAM_VIEW_PL
TMP_PATH = cfl.TMP_PATH
ESL_FSEQ_PATH = cfl.ESL_FSEQ_PATH
ENA_URL = cfl.ENA_URL
TAX_NODES_DUMP = cfl.TAX_NODES_DUMP
TAX_NAMES_DUMP = cfl.TAX_NAMES_DUMP
BIN_LOC = cfl.BIN_LOC
RFAM_SEED_SEQ_14_1 = cfl.RFAM_SEED_SEQ_14_1
PDB_FILES = cfl.PDB_FILES

# ---------------------------------GEN CONFIG----------------------------------
ENA_FTP = cfl.ENA_FTP
SLACK_WEBHOOK = cfl.SLACK_WEBHOOK
USER_EMAIL = cfl.USER_EMAIL
RFAM_EMAIL = cfl.RFAM_EMAIL
BROWSER_HUB_DESC_URL = cfl.BROWSER_HUB_DESC_URL
BULK_REPORT_URL = cfl.BULK_REPORT_URL
APICURON_TOKEN = cfl.APICURON_TOKEN

# ------------------------------DATABASES--------------------------------------
# Dictionaries
RFAMLIVEPUB = cfl.RFAMLIVEPUB
RFAMLIVE = cfl.RFAMLIVE
RFAMLOCAL = cfl.RFAMLOCAL
RFAMLIVELOC = cfl.RFAMLIVELOC
RFAMREL = cfl.RFAMREL
PG = cfl.PG
FB1 = cfl.FB1

# ----------------------------Django settings----------------------------------
# DATABASES
RFAMDEV = cfl.RFAMDEV
RFAMLOC = cfl.RFAMLOC
RFAMLIVE_DJANGO = cfl.RFAMLIVE_DJANGO

# SETTINGS
SECRET_KEY = cfl.SECRET_KEY

# ---------------------------------SEQDBs--------------------------------------
RFAM_SEED_DB = cfl.RFAM_SEED_SEQ_14_1

# -------------------------------LSF------------------------------------
# rfamprod privileges required
FA_EXPORT_GROUP = cfl.FA_EXPORT_GROUP
RFAM_VIEW_GROUP = cfl.RFAM_VIEW_GROUP
LSF_GROUPS_CMD = cfl.LSF_GROUPS_CMD
LSF_GEN_GROUP = cfl.LSF_GEN_GROUP
