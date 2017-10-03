#!/usr/bin/python

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
Description: Support code designed to ease Rfam jiffies' execution.

Notes: 1. Call this script as rfamprod
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import subprocess


# -----------------------------------------------------------------------------


def call_jiffy(jiffy, fam_file, outdir=None):
    """
    This function was designed to call the perl script defined by the jiffy
    parameter. Used in Rfam 12.1 with jiffies writeAnnotatedCM.pl,
    writeAnnotatedSeed.pl and writeAnnotatedTree.pl to generate CM, SEED
    and seed_Tree files for the FTP server

    jiffy:  The path to the perl script to call
    fam_file:   A list of all rfam_accessions
    outdir: Destination directory where the files will be generated
    """

    # TO DO
    # convert this to load family accessions from the database

    if outdir is not None:
        os.chdir(os.path.abspath(outdir))

    fp = open(os.path.abspath(fam_file), 'r')

    for rfam_acc in fp:
        cmd = "%s %s" % (jiffy, rfam_acc)
        subprocess.call(cmd, shell=True)

    fp.close()


# -----------------------------------------------------------------------------


def usage():
    # TO DO
    pass


# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # need to check the number of parameters provided
    jiffy = sys.argv[1]
    fam_file = sys.argv[2]
    outdir = sys.argv[3]

    call_jiffy(jiffy, fam_file, outdir)
