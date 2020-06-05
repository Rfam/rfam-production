import os
import argparse
import random
import datetime

# --------------------------------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--tblout', help="infernal's tblout file", action='store')
    parser.add_argument('--dest-dir', help="destination directory to store output to", action="store")

    return parser

# --------------------------------------------------------------------------------


def convert_tblout_2_pdb_full_region(tblout_file, dest_dir=None):
    """
    Converts Infernal's tblout file to pdb_full_region txt dump
    for direct import to the rfam_live database

    tblout_file: Infernal's tblout file

    return: Void
    """

    if dest_dir is None:
        dest_dir = os.getcwd()

    # Rfam website hex colours to randomly choose from
    hex_colours = ["1fc01f", "c00f0f", "bdc000", "c008ae", "00bac0", "8484c0",
                   "93c090", "c0af92", "8e2511", "f29242", "8585e6", "ff87fa",
                   "008700", "454545", "0003c0", "ebeb30", "ff87a4", "0064f4"]

    rand_hex_colour = random.choice(hex_colours)

    fp_in = open(tblout_file, 'r')

    # create a new pdb_full_region output file and mark the date it has been generated
    fp_out = open(os.path.join(dest_dir, "pdb_full_region_" + str(datetime.date.today()) + ".txt"), 'w')

    new_line = ''

    for line in fp_in:
        # skips comment lines
        if line[0] != '#':

            columns = [x for x in line.strip().split(' ') if x != '']

            pdb_info = columns[2].partition('_')

            new_line = '\t'.join((columns[1], pdb_info[0], pdb_info[2],
                                  columns[7], columns[8], columns[14],
                                  columns[15], columns[5], columns[6],
                                  rand_hex_colour, '1'))

            fp_out.write(new_line+'\n')

    fp_in.close()
    fp_out.close()

# --------------------------------------------------------------------------------


if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    convert_tblout_2_pdb_full_region(args.tblout, dest_dir=args.dest_dir)




