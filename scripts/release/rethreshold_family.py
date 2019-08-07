"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

# ----------------------------------------------------------------------------------

import os
import sys
import subprocess
import argparse

# ------------------------------------- GLOBALS ------------------------------------

LSF_GROUP = "/family_srch"
MEMORY = 8000
CPU = 8

# ----------------------------------------------------------------------------------


def checkout_family(rfam_acc):
    """
    Checks out a family from Rfam based on a valid Rfam accession.

    rfam_acc: A valid Rfam accession

    return: None
    """

    cmd = "rfco.pl %s" % rfam_acc

    subprocess.call(cmd, shell=True)

    # add some checks here

# ----------------------------------------------------------------------------------


def submit_new_rfsearch_job(family_dir):
    """
    Submits a new lsf job that runs rfsearch to update SCORES for a new release.
    If no threshold is set with rfsearch.pl, it uses existing thresholds by default.

    family_dir: The physical location of the family directory

    return: None
    """
    # use the pre-process command to change directory to family_dir

    lsf_err_file = os.path.join(family_dir, "auto_rfsearch.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfsearch.out")

    cmd = ("bsub -M %s -R \"rusage[mem=%s]\" -o %s -e %s -n %s -g %s -R \"span[hosts=1]\" "
          "cd %s && rfsearch.pl -cnompi")

    subprocess.call(cmd % (MEMORY, MEMORY, lsf_out_file, lsf_err_file,
                         CPU, LSF_GROUP, family_dir), shell=True)

# ----------------------------------------------------------------------------------

if __name__ == '__main__':

    pass