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

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import argparse
from utils import db_utils as db

# -----------------------------------------------------------------------------

def generate_clanin_file(dest_dir=None):
    """
    Creates a clanin file to be used for clan competition during cmscan

    dest_dir: The path to destination directory. Using currect if no
    directory provided

    returns: void
    """

    # create destination directory or use current if not provided
    if dest_dir is None:
        dest_dir = os.getcwd()

    else:
        if not os.path.exists(dest_dir):
            os.mkdir(dest_dir)

    # fetch clan members from the database in the form of a dictionary
    clan_members = db.fetch_clanin_data()

    fp = open(os.path.join(dest_dir, 'Rfam.clanin'), 'w')

    for clan in clan_members.keys():
        fp.write(clan+'\t')
        fp.write('\t'.join(clan_members[clan])+'\n')

    fp.close()

# -----------------------------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing using Python's argparse

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser(prog="clanin_file_generator.py")
    parser.add_argument("--dest-dir", help="Destination directory to store output to")

    return parser

# -----------------------------------------------------------------------------


if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    generate_clanin_file(args.dest_dir)
