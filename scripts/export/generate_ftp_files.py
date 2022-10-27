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
import shutil
import subprocess
import sys
import tempfile

from utils.db_utils import fetch_all_rfam_accs

JIFFIES = {"SEED": "writeAnnotatedSeed.pl",
           "CM": "writeAnnotatedCM.pl",
           "TREE": "writeAnnotatedTree.pl"}


# -----------------------------------------------------------------------------


def export_ftp_file(file_type, rfam_acc, dest_dir=None):
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
        jiffy = JIFFIES[file_type]

        # call jiffy
        cmd = "%s %s" % (jiffy, rfam_acc)
        subprocess.call(cmd, shell=True)

    except:
        raise Exception('Command {} failed'.format(cmd))


# -----------------------------------------------------------------------------


def rename_files(source_dir, file_type, rfam_acc):
    """

    :param source_dir:
    :param file_type:
    :return:
    """

    if not os.path.exists(source_dir):
        raise IOError

    if file_type == "SEED":
        seed_file_loc = os.path.join(source_dir, rfam_acc)

        if not os.path.exists(seed_file_loc):
            sys.exit("File does not exist %s" % seed_file_loc)

        new_name = os.path.join(source_dir, rfam_acc + '.seed')
        os.rename(seed_file_loc, new_name)

        if not os.path.exists(new_name):
            sys.exit("%s SEED cound not be renamed" % rfam_acc)
    elif file_type == "TREE":
        taxtree = os.path.join(source_dir, rfam_acc + '.taxtree')
        seed_tree = os.path.join(source_dir, rfam_acc + '.seed_tree')
        shutil.copy(taxtree, seed_tree)


# -----------------------------------------------------------------------------


def create_seed_archive(destination):
    """
    Create a combined Rfam.seed file and compress it.
    """
    cwd = os.getcwd()
    os.chdir(destination)
    cmd = "rm -f Rfam.seed && cat *.seed > Rfam.seed && gzip -c Rfam.seed > Rfam.seed.gz"
    status = os.system(cmd.format(destination))
    if status:
        raise Exception('There was a problem generating Rfam.seed.gz in {}'.format(destination))
    os.chdir(cwd)


# -----------------------------------------------------------------------------


def create_combined_cm_file(destination):
    """
    Create a combined Rfam.cm file.
    """
    cwd = os.getcwd()
    os.chdir(destination)
    cmd = "rm -f Rfam.cm && cat *.CM > Rfam.cm"
    status = os.system(cmd.format(destination))
    if status:
        raise Exception('There was a problem generating Rfam.cm in {}'.format(destination))
    os.chdir(cwd)


# -----------------------------------------------------------------------------


def create_tree_archive(destination):
    """
    Create a combined Rfam.seed_tree file.
    """
    cwd = os.getcwd()
    os.chdir(destination)
    cmd = ("rm -f Rfam.seed_tree.tar.gz && "
           "rm -Rf Rfam.seed_tree && "
           "mkdir Rfam.seed_tree && "
           "mv RF0*.seed_tree Rfam.seed_tree && "
           "tar -cf Rfam.seed_tree.tar.gz Rfam.seed_tree")
    status = os.system(cmd.format(destination))
    if status:
        raise Exception('There was a problem generating Rfam.seed_tree.gz in {}'.format(destination))
    os.chdir(cwd)


# -----------------------------------------------------------------------------


def validate_seed_archive(destination, rfam_accs):
    """
    Check that Rfam.seed contains the correct number of entries.
    """
    cwd = os.getcwd()
    os.chdir(destination)
    family_count = 0
    with open(os.path.join(destination, 'Rfam.seed')) as f_seed:
        for line in f_seed:
            if line.startswith('# STOCKHOLM 1.0'):
                family_count += 1
    os.chdir(cwd)
    try:
        assert (family_count == len(rfam_accs))
        print('OK: Found {} families in Rfam.seed'.format(family_count))
    except AssertionError:
        raise Exception('Error: Rfam.seed contains {} families instead of {}'.format(family_count, len(rfam_accs)))


# -----------------------------------------------------------------------------


def get_all_rfam_accessions():
    """
    Fetch a list of all Rfam families from the SVN repository.
    """
    rfam_accessions = []
    svn_url = 'https://xfamsvn.ebi.ac.uk/svn/data_repos/trunk/Families/'
    svn_list = tempfile.NamedTemporaryFile()
    cmd = "svn list {} > {}".format(svn_url, svn_list.name)
    os.system(cmd)
    with open(svn_list.name, 'r') as f_svn_list:
        for line in f_svn_list:
            if line.startswith('RF'):
                rfam_accessions.append(line.strip().replace('/', ''))
    print('Found {} accessions on SVN'.format(len(rfam_accessions)))
    return rfam_accessions


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
    if args.acc == 'all':
        accessions = fetch_all_rfam_accs()
    elif os.path.isfile(args.f):
        fp = open(args.f, 'r')
        accessions = [acc.strip() for acc in fp]
        fp.close()
    elif args.acc[0:2] == 'RF':
        accessions.append(args.acc)

    if args.dest_dir and not os.path.exists(args.dest_dir):
        os.system('mkdir -p {}'.format(args.dest_dir))
        print('Created folder {}'.format(args.dest_dir))

    jiffy_type = None

    if args.seed:
        jiffy_type = "SEED"
    elif args.cm:
        jiffy_type = "CM"
    elif args.tree:
        jiffy_type = "TREE"

    for accession in accessions:
        print(accession)
        export_ftp_file(jiffy_type, accession, args.dest_dir)
        rename_files(args.dest_dir, jiffy_type, accession)

    if args.seed and args.acc == 'all':
        create_seed_archive(args.dest_dir)
        validate_seed_archive(args.dest_dir, accessions)
    elif args.cm and args.acc == 'all':
        create_combined_cm_file(args.dest_dir)
    elif args.tree and args.acc == 'all':
        create_tree_archive(args.dest_dir)
