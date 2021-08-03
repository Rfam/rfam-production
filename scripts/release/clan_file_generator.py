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

import os
import sys
import argparse
from utils import db_utils as rfamdb

# loop over all files and sort using unix sort command to group duplicate seq_acc
# together
# for file in ./*; do sort -k2 -t $'\t' $file > ${file:2:7}s.txt; done
# need to add order by clause in the query to omit linux sort

# -----------------------------------------------------------------------------


def parse_clan_file(clan_list):
    """
    Parses a list of Rfam clan accessions

    clan_list: A plain .txt file containing a list of Rfam Clan Accessions

    return: A list of clan accessions
    """

    fp = open(clan_list, 'r')

    clan_accessions = [x.strip() for x in fp]

    fp.close()

    return clan_accessions

# -----------------------------------------------------------------------------


def clan_file_generator(output_dir, clan_comp_type='FULL', clan_acc=None):
    """
    Generates clan files for clan competition

    output_dir: The path to the output directory. It will be created if
    it does not exist
    clan_comp_type: This can be 'FULL' for clan competition on on full_region_table
    or PDB for clan competition on pdb_full_region

    returns: void
    """

    if os.path.exists(output_dir) is False:
        os.mkdir(output_dir)

    # fetch all clans
    clans = []
    if clan_acc is None:
        clans = rfamdb.fetch_clan_accessions()

    elif os.path.isfile(clan_acc):
        # fetches clan accessions from file
        clans = parse_clan_file(clan_acc)

    elif clan_acc[0:2] == "CL":
        clans.append(clan_acc)

    else:
        sys.exit("Wrong clan input provided!")

    # some sanity checks
    if clan_comp_type.upper() != 'FULL' and clan_comp_type.upper() != 'PDB':
        sys.exit("\nPlease provide correct type for clan competition!")

    # generate clan files
    for clan in clans:
        regions = []
        if clan_comp_type.upper() == 'FULL':
            regions = rfamdb.fetch_clan_full_region_records(clan)

        elif clan_comp_type.upper() == 'PDB':
            regions = rfamdb.fetch_clan_pdb_full_region_records(clan)

        # open a new clan file to dump the sequences
        if len(regions) > 0:
            print("clan: ", clan)
            clan_fp = open(os.path.join(output_dir, clan) + ".txt", 'w')

        # if no available regions, move to next clan
        else:
            continue

        flag = 1
        # loop over rows
        for region in regions:
            # append fields
            if clan_comp_type == 'FULL':
                region = region[:len(region) - 1]
            line = ''
            for field in region:
                if flag == 1:
                    line = str(field)
                    flag = 0
                else:
                    line = line + '\t' + str(field)

            line += '\n'

            clan_fp.write(line)
            flag = 1

        clan_fp.close()

        flag = 1

# -----------------------------------------------------------------------------


def parse_arguments():
    """
    Performs some basic argument parsing

    return: parser object
    """

    parser = argparse.ArgumentParser(description='Generates required clan competition input files')

    mutualy_exclusive = parser.add_mutually_exclusive_group()
    parser.add_argument("--dest-dir", help="Destination directory where to generate the files",
                        type=str, action="store")
    parser.add_argument("--cc-type", help="Specifies clan competition type FULL/PDB", type=str,
                        action="store")
    mutualy_exclusive.add_argument("--clan-acc", help="Clan accession", type=str)
    mutualy_exclusive.add_argument("--all", help="All clans in Rfam",
                                   action="store_true", default=False)
    mutualy_exclusive.add_argument("-f",
                                   help="A file with a list of clan accessions to compete",
                                   action="store")

    return parser

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # arg1: path to output directory
    # arg2: clan competition type 'FULL' for full_region 'PDB' for pdb_full_region
    # arg3: clan accession if we only want to generate a file for a single clan

    parser = parse_arguments()
    args = parser.parse_args()

    if args.clan_acc:
        clan_file_generator(args.dest_dir, args.cc_type, args.clan_acc)

    elif args.all:
            clan_file_generator(args.dest_dir, args.cc_type, None)

    else:
        clan_file_generator(args.dest_dir, args.cc_type, args.f)




