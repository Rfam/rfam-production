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

# ---------------------------------GEN_CONFIG----------------------------------

RFAM_GPFS_LOC = ''
LOC_PATH = ''
GEN_DWLD_EXEC = ''

LSF_GROUPS_CMD = 'bgadd -L %s /rfam_gen/%s'
LSF_GEN_GROUP = '/rfam_gen'

USER_EMAIL = ''


# ------------------------------DATABASES--------------------------------------
# Databases
RFAMLIVEPUB = {
    'user': '',
    'pwd': '',
    'host': '',
    'db': '',
    'port': '',
}

RFAMLIVE = {
    'user': '',
    'pwd': '',
    'host': '',
    'db': '',
    'port': '',
}


RFAM12 = {
    'user': '',
    'pwd': '',
    'host': '',
    'db': '',
    'port': '',
}

RFAMLIVELOC = {
    'user': '',
    'pwd': '',
    'host': '',
    'db': '',
    'port': '',
}

# ----------------------------Django settings----------------------------------

# DATABASES
RFAMDEV = {
    'ENGINE': 'django.db.backends.mysql',
    'NAME': '',
    'HOST': '',
    'PORT': '',
    'USER': '',
    'PASSWORD': '',
}

RFAMLOC = {
    'ENGINE': 'django.db.backends.mysql',
    'NAME': '',
    'HOST': '',
    'PORT': '',
    'USER': '',
    'PASSWORD': '',
}

# SETTINGS
SECRET_KEY = 'change secret key in production'

# ----------------------------RFAM CONFIG PATHS--------------------------------

ESL_PATH = ''
FA_GEN = ''
RFAMSEQ_PATH = ''
FAM_VIEW_PL = ''
TMP_PATH = '/tmp'

ESL_FSEQ_PATH = ''
FSR_PATH = ''
FSR_LOCAL = ''
ENA_URL = 'http://www.ebi.ac.uk/ena/data/view/%s&display=fasta&range=%s-%s'

# Maybe delete these
TAX_NODES_DUMP = ''
TAX_NAMES_DUMP = ''
RFAM_NCBI_IDS = ''
VALID_NCBI_IDS = ''
NCBI_RANKS = ''

# -------------------------------LSF GROUPS------------------------------------
# rfamprod privileges required
FA_EXPORT_GROUP = '/rfam_fa'
RFAM_VIEW_GROUP = '/rfam_view'

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    pass
