#!/usr/bin/python3

import os
import argparse
import requests


# --------------------------------------------------------------------


def desc_template_generator(desc_file, mirna_name, family_id, wiki_links=None, second_author=None, go_terms=None,
                            dest_dir=None):
    """

    desc_file:
    mirna_name:
    family_id:
    wiki_links:
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

    wk_constant = "MicroRNA"

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

    if second_author is not None:
        if second_author.find("Griffiths-Jones SR") == -1:
            author = author + '; ' + second_author

    desc_template = """ID   %s
DE   %s microRNA precursor family
AU   %s
SE   Griffiths-Jones SR
SS   Predicted; RNAalifold
GA   %s
TC   %s
NC   %s
TP   Gene; miRNA;
BM   %s
CB   %s
SM   %s
DR   MIPF; %s;
DR   URL; http://www.mirbase.org;
DR   SO; 0001244; pre_miRNA;"""

    desc_bottom_half = """
CC   This family represents the microRNA (miRNA) precursor %s
CC   imported from miRBase.
WK   %s
"""

    fp_out = open(os.path.join(dest_dir, "DESC"), 'w')

    wiki_link = wk_constant

    if wiki_links and family_id in wiki_links:
        wiki_link = wiki_links[family_id]

    # write DESC upper half
    fp_out.write(desc_template % (mirna_name, mirna_name, author, essential_desc_lines["GA"],
                                  essential_desc_lines["TC"], essential_desc_lines["NC"],
                                  essential_desc_lines["BM"], essential_desc_lines["CB"],
                                  essential_desc_lines["SM"], family_id))

    # write GO terms if any
    if go_terms is not None:
        i = 0

        for go_term_string in go_terms.keys():
            if i < len(go_terms.keys()) - 1:
                if i == 0:
                    fp_out.write("\nDR   %s\n" % go_term_string)
                else:
                    fp_out.write("DR   %s\n" % go_term_string)
            else:
                fp_out.write("DR   %s" % go_term_string)
            i += 1

    # write DESC bottom half
    fp_out.write(desc_bottom_half % (mirna_name, wiki_link))

    fp_out.close()

    if os.path.exists(os.path.join(dest_dir, "DESC")):
        return True

    return False


# --------------------------------------------------------------------


def extract_sequence_accessions_from_seed(seed_file):
    """
    Parses a seed MSA and extracts all sequence
    accessions in the form of a dictionary

    seed_file: An Rfam seed alignment

    return: A dictionary of seed accessions
    """

    accessions = {}
    fp = open(seed_file, 'r')

    for line in fp:
        line = line.strip()
        if len(line) > 1 and line[0] != '#' and line != '':
            line = line.split(' ')
            accession = line[0].partition('/')[0]
            if accession != '':
                accessions[accession] = ""

    fp.close()

    return accessions


# --------------------------------------------------------------------


def fetch_go_terms_from_rnacentral(rnacentral_id):
    """

    rnacentral_id:
    return:
    """

    url = "https://rnacentral.org/api/v1/rna/%s/go-annotations/%s"

    id_elements = rnacentral_id.partition("_")
    upi = id_elements[0]
    taxid = id_elements[2]

    response = requests.get(url % (upi, taxid))

    go_terms = {}

    if response.status_code == 200:
        data = response.json()
        if len(data) == 0:
            return None
        else:
            for item in data:
                if item["go_term_id"] not in go_terms:
                    go_terms[item["go_term_id"]] = item["go_term_name"]

    return go_terms


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
    parser.add_argument("--wiki-links", help="File with family/wiki link mappings", action="store", default=None)
    parser.add_argument("-f", help="a list of accessions to generate DESC files for", action="store", default=None)
    parser.add_argument("--label", help="label to be used for the generation of the DESC file", action='store',
                        default=None)
    parser.add_argument("--no-ref", help="Skips references in DESC file", action="store_true", default=False)

    return parser.parse_args()


# --------------------------------------------------------------------

if __name__ == '__main__':

    args = parse_arguments()
    desc_file = os.path.join(args.input, "DESC")

    if args.label is None:
        label = os.path.basename(args.input)
    else:
        label = args.label

    mirna_labels = [x for x in label.split("_") if x != '']
    mirna_family_id = mirna_labels[0]
    mirna_name = mirna_labels[1]

    wiki_links = None
    if args.wiki_links is not None:
        wiki_links = {}
        fp = open(args.wiki_links, 'r')
        for line in fp:
            line = line.strip().split('\t')
            if line[0] not in wiki_links:
                wiki_links[line[0]] = line[2]
        fp.close()

    seed_file = os.path.join(args.input, "SEED")

    seed_accessions = extract_sequence_accessions_from_seed(seed_file)

    family_go_terms = {}

    for urs_accession in seed_accessions:
        urs_go = fetch_go_terms_from_rnacentral(urs_accession)
        if urs_go is not None:
            for go_id in urs_go.keys():
                go_term = '; '.join([go_id.replace(':', '; '), urs_go[go_id] + ';'])
                if go_term not in family_go_terms:
                    family_go_terms[go_term] = ""

    desc_template_generator(desc_file, mirna_name, mirna_family_id, wiki_links, args.ga_author, family_go_terms,
                            dest_dir=args.outdir)
