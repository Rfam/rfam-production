#!/usr/bin/python

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

import argparse
import os
import subprocess
import sys

JIFFIES = {"SEED": "writeAnnotatedSeed.pl",
           "CM": "writeAnnotatedCM.pl",
           "TREE": "writeAnnotatedTree.pl"}

# -----------------------------------------------------------------------------


def export_ftp_file(type, rfam_acc, dest_dir=None):
    """

    :param type:
    :param rfam_acc:
    :param dest_dir:
    :return:
    """

    # TO DO
    # convert this to load family accessions from the database

    if dest_dir is not None:
        os.chdir(os.path.abspath(dest_dir))
    else:
        sys.exit("Provide a valid destination directory")

    try:
        # select jiffy
        jiffy = JIFFIES[type]

        # call jiffy
        cmd = "%s %s" % (jiffy, rfam_acc)
        subprocess.call(cmd, shell=True)

    except:
        raise OSError

# -----------------------------------------------------------------------------


def parse_arguments():
    """
    """

    parser = argparse.ArgumentParser()

    mutually_exclusive_type = parser.add_mutually_exclusive_group()
    mutually_exclusive_type.add_argument("--seed", help="Export SEED files from SVN",
                                         action="store_true", default=False)
    mutually_exclusive_type.add_argument("--cm", help="Export CM files from SVN",
                                         action="store_true", default=False)
    mutually_exclusive_type.add_argument("--tree", help="Export SEED files from SVN",
                                         action="store_true", default=False)

    mutually_exclusive_input = parser.add_mutually_exclusive_group()
    mutually_exclusive_input.add_argument("--acc", help="Rfam family accession",
                                          action="store", type=str)
    mutually_exclusive_input.add_argument("-f", help="List of Rfam family accessions (.txt)",
                                          action="store", type=str)

    parser.add_argument("--dest-dir", help="Destination directory to generate files to")

    return parser

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    accessions = []
    if os.path.isfile(args.f):
        fp = open(args.f, 'r')
        accessions = [acc.strip() for acc in fp]
        fp.close()
    elif args.acc[0:2] == 'RF':
        accessions.append(args.acc)

    jiffy_type = None

    if args.seed:
        jiffy_type = "SEED"
    elif args.cm:
        jiffy_type = "CM"
    elif args.tree:
        jiffy_type = "TREE"

    for accession in accessions:
        export_ftp_file(jiffy_type, accession, args.dest_dir)
