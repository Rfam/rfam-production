#!/usr/lib/env python

import argparse


# ----------------------------------------------------------

def fasta_headers_to_urs_accessions(fasta_header_file):
    """
    Extracts the URS accessions from a file containing
    all zwd fasta header lines. Fasta headers can be
    extracted using grep '>' new_zwd.fasta > fasta_header_file

    fasta_header_file: A .txt file containing all header lines
    from a new ZWD fasta file

    returns: A dictionary with all URS accessions as keys
    """

    urs_accs = {}
    fp = open(fasta_header_file, 'r')

    for line in fp:
        # split header line on spaces and trip first '>' character
        urs_acc = line.strip().split(' ')[0][1:]
        if urs_acc not in urs_accs:
            urs_accs[urs_acc] = ""

    fp.close()

    return urs_accs


# ----------------------------------------------------------

def load_rfam_urs_accessions_from_file(urs_acc_list):
    """
    Loads all existing Rfam URS accessions in a python
    dictionary

    urs_acc_list: A .txt file with all URS accession already
    in Rfam

    return: A python dictionary with all URS accessions as
    keys.
    """

    rfam_urs_accs = {}

    fp = open(urs_acc_list, 'r')

    for line in fp:
        accession = line.strip()
        if accession not in rfam_urs_accs:
            rfam_urs_accs[accession] = ""

    fp.close()

    return rfam_urs_accs


# ----------------------------------------------------------

def get_novel_zwd_accessions(rfam_urs_accs, new_zwd_accs):
    """
    Isolates  all new ZWD URS accessions to imported to Rfam
    with the use of python sets

    rfam_urs_accs: A python dictionary with all Rfam URS accessions
    new_zwd_accs: A python dictionary with all new ZWD URS candidates

    return: A list of novel ZWD URS accessions for import to Rfam
    """

    novel_zwd_accessions = list(set(new_zwd_accs.keys()).difference(set(rfam_urs_accs.keys())))

    return novel_zwd_accessions


# ----------------------------------------------------------

def get_novel_zwd_accessions_from_dict(rfam_urs_accs, new_zwd_accs):
    """
    Isolates  all new ZWD URS accessions to imported to Rfam
    with the use of python dictionaries

    rfam_urs_accs: A python dictionary with all Rfam URS accessions
    new_zwd_accs: A python dictionary with all new ZWD URS candidates

    return: A list of novel ZWD URS accessions for import to Rfam
    """

    novel_zwd_accessions = []

    for accession in new_zwd_accs.keys():
        if accession not in rfam_urs_accs:
            novel_zwd_accessions.append(accession)

    return novel_zwd_accessions


# ----------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing using python's argparse library
    """

    parser = argparse.ArgumentParser(description="Checks for novel ZWD accessions")

    parser.add_argument("--zwd-headers",
                        help="A header file generated directly from ZWD fasta", action="store")
    parser.add_argument("--rfam-urs-list",
                        help="A file listing all ZWD URS accessions in Rfam", action="store")
    # parser.add_argument("--dest-dir",
    # help="Destination directory where output will be stored", action="store")

    parser.add_argument("--type", help="Type of output [URS | TAXIDS]", action="store")

    return parser


# ----------------------------------------------------------

if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    new_zwd_accs = fasta_headers_to_urs_accessions(args.zwd_headers)
    rfam_accs = load_rfam_urs_accessions_from_file(args.rfam_urs_list)

    novel_accessions = get_novel_zwd_accessions(rfam_accs, new_zwd_accs)

    if len(novel_accessions) > 0:

        if args.type.upper() == 'URS':
            for accession in sorted(novel_accessions):
                print (accession)

        # print taxids for taxonomy table update
        elif args.type.upper() == 'TAXIDS':
            taxids = []
            for accession in sorted(novel_accessions):
                taxids.append(accession.partition('_')[2])
            # fetch unique taxids
            unique_taxids = list(set(taxids))
            for taxid in unique_taxids:
                    print (taxid)

    else:
        exit("\nThere are no new ZWD sequences to import!\n")
