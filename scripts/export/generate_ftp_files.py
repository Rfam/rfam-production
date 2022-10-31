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

import argparse
import os
import shutil
import subprocess
import sys

from utils.db_utils import fetch_all_rfam_accs

JIFFIES = {"SEED": "writeAnnotatedSeed.pl",
           "CM": "writeAnnotatedCM.pl",
           "TREE": "writeAnnotatedTree.pl"}


def export_ftp_file(file_type, rfam_acc, dest_dir=None):
    """
    Call the perl jiffy script (cm, seed, or tree) to export the ftp file
    :param file_type: cm, seed, or tree
    :param rfam_acc: Rfam ID
    :param dest_dir: where to write ftp file to
    """

    if dest_dir is not None:
        os.chdir(os.path.abspath(dest_dir))
    else:
        sys.exit("Provide a valid destination directory")

    try:
        jiffy = JIFFIES[file_type]
        cmd = "%s %s" % (jiffy, rfam_acc)
        subprocess.call(cmd, shell=True)
    except Exception as e:
        raise Exception('Command {cmd} failed: {err}'.format(cmd=cmd, err=e))


def rename_files(source_dir, file_type, rfam_acc):
    """
    Rename the ftp files

    :param source_dir: source of FTP file
    :param file_type: cm, seed, or tree
    :param rfam_acc: Rfam ID
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


def create_combined_cm_file(destination, accs):
    """
    Create a combined Rfam.cm file.
    Write a .txt file containing a list of all Rfam accessions in the CM file.
    """
    cwd = os.getcwd()
    os.chdir(destination)
    with open("list.txt") as accs_list:
        accs_list.write(accs)
    cmd = "rm -f Rfam.cm && cat *.CM > Rfam.cm"
    status = os.system(cmd.format(destination))
    if status:
        raise Exception('There was a problem generating Rfam.cm in {}'.format(destination))
    os.chdir(cwd)


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


def parse_arguments():
    """
    """
    arg_parser = argparse.ArgumentParser()

    mutually_exclusive_type = arg_parser.add_mutually_exclusive_group()
    mutually_exclusive_type.add_argument("--seed", help="Export SEED files from SVN",
                                         action="store_true", default=False)
    mutually_exclusive_type.add_argument("--cm", help="Export CM files from SVN",
                                         action="store_true", default=False)
    mutually_exclusive_type.add_argument("--tree", help="Export SEED files from SVN",
                                         action="store_true", default=False)

    mutually_exclusive_input = arg_parser.add_mutually_exclusive_group()
    mutually_exclusive_input.add_argument("--acc", help="Rfam family accession",
                                          action="store", type=str)
    mutually_exclusive_input.add_argument("-f", help="List of Rfam family accessions (.txt)",
                                          action="store", type=str)

    arg_parser.add_argument("--dest-dir", help="Destination directory to generate files to")

    return arg_parser


if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    accessions = []
    if args.acc == 'all':
        accessions = fetch_all_rfam_accs()
    elif args.f:
        if os.path.isfile(args.f):
            with open(args.f, 'r') as fp:
                accessions = [acc.strip() for acc in fp]
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
        create_combined_cm_file(args.dest_dir, accessions)
    elif args.tree and args.acc == 'all':
        create_tree_archive(args.dest_dir)
