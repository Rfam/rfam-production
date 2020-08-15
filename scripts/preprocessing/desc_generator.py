#!/usr/bin/python3

import os
import argparse

# --------------------------------------------------------------------


def desc_template_generator(desc_file, mirna_name, family_id, second_author=None, dest_dir=None):
    """

    desc_file:
    mirna_name:
    family_id:
    dest_dir:
    return:
    """

    """
    ID   ShortName
    DE   Family description
    AU   Who RU
    SE   Where did the seed come from
    GA   25.00
    TC   30.10
    NC   24.50
    BM   cmbuild -F CM SEED
    CB   cmcalibrate --mpi CM
    SM   cmsearch --cpu 4 --verbose --nohmmonly -T 30.00 -Z 742849.287494 CM SEQDB
    """

    essential_desc_lines = {"GA": "", "TC": "", "NC": "", "BM": "", "CB": "", "SM": ""}

    if dest_dir is None:
        dest_dir = os.path.split(desc_file)[0]

    if desc_file is not None:
        fp_in = open(desc_file, 'r')
        for line in fp_in:
            if line[0:2] in essential_desc_lines:
                line = [x for x in line.strip().split(' ') if x != '']
                new_line = ' '.join(line[1:])
                essential_desc_lines[line[0]] = new_line
        fp_in.close()

    os.rename(desc_file, os.path.join(dest_dir, "DESC_old"))

    author = "Griffiths-Jones SR; 0000-0001-6043-807X"

    if second_author is not None or second_author.find("Griffiths-Jones SR") == -1:
        author = author + '; ' + second_author

    desc_template = """ID   %s
DE   %s microRNA precursor family
AU   %s
SE   Griffiths-Jones SR
SS   Predicted; PFOLD
GA   %s
TC   %s
NC   %s
TP   Gene; miRNA;
BM   %s
CB   %s
SM   %s
DR   MIPF; %s;
DR   URL; http://www.mirbase.org;
DR   SO; 0001244; pre_miRNA;
CC   This family represents the microRNA (miRNA) precursor %s
CC   imported from miRBase.
WK   MicroRNA
"""

    fp_out = open(os.path.join(dest_dir, "DESC"), 'w')

    fp_out.write(desc_template % (mirna_name, mirna_name, author, essential_desc_lines["GA"],
                                  essential_desc_lines["TC"], essential_desc_lines["NC"],
                                  essential_desc_lines["BM"], essential_desc_lines["CB"],
                                  essential_desc_lines["SM"], family_id, mirna_name))

    fp_out.close()

    if os.path.exists(os.path.join(dest_dir, "DESC")):
        return True

    return False


# --------------------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing using python's argparse

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser("Generates a DESC template for a new family")
    parser.add_argument("--input", help="miRBase directory with rfsearch results",
                      action="store", default=None)
    parser.add_argument("--outdir", help="Path to the output directory", action="store",
                      default=None)
    parser.add_argument("--ga-author", help="Name and orcid of gathering threshold author (e.g. Edwards BA; ORCID)",
                        action="store", default=None)

    return parser

# --------------------------------------------------------------------

if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    desc_file = os.path.join(args.input, "DESC")

    mirna_labels = [x for x in os.path.basename(args.input).split("_") if x != '']
    mirna_family_id = mirna_labels[0]
    mirna_name = mirna_labels[1]

    desc_template_generator(desc_file, mirna_name, mirna_family_id, args.ga_author,
                            dest_dir=args.outdir)
