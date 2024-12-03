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

# ----------------------------PATHS--------------------------------
ESL_PATH = ""
ESL_SEQSTAT = ""
ESL_SSPLIT = ""
CMSEARCH = ""
CMSCAN = ""
CMFILE = ""
SCANNER = ""
DWL_SCRIPT = ""
FA_GEN = ""
RFAMSEQ_PATH = ""
FAM_VIEW_PL = ""
TMP_PATH = "/tmp"
ESL_FSEQ_PATH = ""
ENA_URL = "http://www.ebi.ac.uk/ena/data/view/%s&display=fasta&range=%s-%s"
TAX_NODES_DUMP = ""
TAX_NAMES_DUMP = ""
BIN_LOC = ""
RFAM_SEED_SEQ_14_1 = ""
PDB_FILES = ""

# ---------------------------------GEN CONFIG----------------------------------
ENA_FTP = ""
SLACK_WEBHOOK = ""
USER_EMAIL = ""
RFAM_EMAIL = ""
BROWSER_HUB_DESC_URL = ""
BULK_REPORT_URL = ""
APICURON_TOKEN = ""

# ------------------------------DATABASES--------------------------------------
RFAMLIVEPUB = {
    "user": "",
    "pwd": "",
    "host": "",
    "db": "",
    "port": "",
}

RFAMLIVE = {
    "user": "",
    "pwd": "",
    "host": "",
    "db": "",
    "port": "",
}

RFAMLIVELOC = {
    "user": "",
    "pwd": "",
    "host": "",
    "db": "",
    "port": "",
}

# ----------------------------Django settings----------------------------------
# DATABASES
RFAMLIVE_DJANGO = {
    "USER": RFAMLIVE["user"],
    "PASSWORD": RFAMLIVE["pwd"],
    "HOST": RFAMLIVE["host"],
    "NAME": RFAMLIVE["db"],
    "PORT": RFAMLIVE["port"],
    "ENGINE": "django.db.backends.mysql",
}
RFAMDEV = {
    "ENGINE": "django.db.backends.mysql",
    "NAME": "",
    "HOST": "",
    "PORT": "",
    "USER": "",
    "PASSWORD": "",
}

RFAMLOC = {
    "ENGINE": "django.db.backends.mysql",
    "NAME": "",
    "HOST": "",
    "PORT": "",
    "USER": "",
    "PASSWORD": "",
}

# SETTINGS
SECRET_KEY = "change secret key in production"

# ---------------------------------SEQDBs--------------------------------------
RFAM_SEED_DB = ""

# -------------------------------LSF------------------------------------
# rfamprod privileges required
FA_EXPORT_GROUP = "/rfam_fa"
RFAM_VIEW_GROUP = "/rfam_view"
LSF_GROUPS_CMD = "bgadd -L %s /rfam_gen/%s"
LSF_GEN_GROUP = "/rfam_gen"

# -----------------------------------------------------------------------------

if __name__ == "__main__":
    pass
