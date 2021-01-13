import os
import argparse

# ----------------------------------------------------------------


def parse_outlist_file(outlist_file):
    """

    :param outlist_file:
    :return:
    """
    outlist_info = {}

    fp = open(outlist_file, 'r')

    for line in fp:
        # if not a comment line
        if line[0] != '#':
            line = line.strip().split()
            unique_accession = "_".join([line[3], line[5], line[6]])
            if unique_accession not in outlist_info:
                outlist_info[unique_accession] = {"evalue": float(line[1]),
                                                  "bit_score": float(line[0]),
                                                  "accession": line[3],
                                                  "start": int(line[5]),
                                                  "end": int(line[6])}

    fp.close()

    return outlist_info


# ----------------------------------------------------------------


def extract_tax_ids_from_species(species_file):
    """

    :param species_file:
    :return:
    """

    tax_ids = {}

    fp = open(species_file, 'r')

    for line in fp:
        # if not a comment line
        if line[0] != '#':
            line = line.strip().split()

            if line[3] not in tax_ids:
                tax_ids[line[3]] = int(line[5])

    fp.close()

    return tax_ids

# ----------------------------------------------------------------


def parse_arguments():
    """

    :return:
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--family-dir", help="Path to family directory")

    return parser

# ----------------------------------------------------------------

if __name__ == "__main__":

    parser = parse_arguments()
    args = parser.parse_args()

    species_file = os.path.join(args.family_dir, "species")
    outlist_file = os.path.join(args.family_dir, "outlist")

    outlist_info = parse_outlist_file(outlist_file)
    taxids = extract_tax_ids_from_species(species_file)

    

