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

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import subprocess
from utils import db_utils as db

# ----------------------------------------------------------------------------


def group_cms(cm_dir, no_of_cms=6, dest_dir=None):
    """
    Fetch all family accessions from the database and sort in DESC order according to
    seed size. Split the cms into multiple files defined by no_of_cms making sure that
    we split the large families across the multiple cms

    cm_dir: A directory with all Rfam single covariance model files
    no_of_cms: The number of sub covariance model files to produce
    dest_dir: The path to the destination directory. Will use cm_dir if this parameter
    is None
    """

    rfam_accs = db.fetch_rfam_accs_sorted(order='DESC')

    rfam_cms = [x for x in os.listdir(cm_dir) if x.endswith('.cm')]

    if dest_dir is None:
        dest_dir = cm_dir

    idx = 1
    for cm_file in rfam_cms:
        cm_file_loc = os.path.join(cm_dir, cm_file)

        cmd = "cat %s >> %s" % (cm_file_loc, os.path.join(dest_dir, 'CM'+str(idx)))
        subprocess.call(cmd, shell=True)

        idx += 1

        if idx > no_of_cms:
            idx = 1

# ----------------------------------------------------------------------------

if __name__ == '__main__':

    cm_dir = sys.argv[1]

    group_cms(cm_dir, no_of_cms=6, dest_dir=None)








