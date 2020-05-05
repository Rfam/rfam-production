import os
import xml.etree.ElementTree as ET
import requests
import datetime
import argparse

from utils import db_utils as db

# ------------------------------------------------------------------------------


def fetch_taxid_from_ncbi(accession):
    """
    Uses a NCBI sequence accession and to extract its
    corresponding tax id from NCBI

    accession: A valid NCBI genome accession

    returns: Tax id if available, None otherwise
    """

    cmd = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=%s"

    # make a call to the API
    request = requests.get(cmd % accession)
    # check that everything went alright
    if request.status_code == 200:
        # fetch request text
        text = request.text
        # convert text to xml
        root = ET.fromstring(text)

        if root is not None:
            docsum = root.find("DocSum")

            if docsum is not None:
                items = docsum.findall("Item")

                if items is not None:
                    for item in items:
                        if item.get("Name") == "TaxId":
                            return item.text

    return None

# ------------------------------------------------------------------------------

def fetch_genome_metadata_from_ncbi(accession):
    """
    Fetches genome related metadata from NCBI using a
    valid sequence accession.

    accession: A valid NCBI genome accession

    return: A python dictionary with genome related metadata
    """

    cmd = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=%s"

    genome_metadata = {}

    # make a call to the API
    request = requests.get(cmd % accession)
    # check that everything went alright
    if request.status_code == 200:
        # fetch request text
        text = request.text
        # convert text to xml
        root = ET.fromstring(text)

        if root is not None:
            docsum = root.find("DocSum")

            if docsum is not None:
                items = docsum.findall("Item")

                if items is not None:
                    for item in items:
                        if item.get("Name") == "Title":
                            genome_metadata["description"] = item.text
                        if item.get("Name") == "TaxId":
                            genome_metadata["ncbi_id"] = item.text
                        if item.get("Name") == "Length":
                            genome_metadata["length"] = item.text

        return genome_metadata

    return None

# ------------------------------------------------------------------------------


def padding_rfam_genome_id_with_zeros(number):
    """
    Adds a prefix of zeros to the number provided to
    generate a new Rfam genome id

    number: A number to add zeros to

    return: A 9 character long string of a number with
    a prefix of zeros
    """

    zero_vector = "000000000"

    text_num = str(number)

    num_prefix_zeros = len(zero_vector) - len(text_num)

    padded_number = zero_vector[0:num_prefix_zeros] + text_num

    return padded_number

# ------------------------------------------------------------------------------


def generate_new_rfam_genome_accession(previous_accession = None):
    """
    Generates a new unique Rfam genome accession for a new RfamLive genome
    table entry

    previous_accession: The previous accession if generating new Rfam
    genome accessions for multiple entries

    return: A new genome accession
    """

    # fetch the maximum RG accession from RfamLive or use the previous
    if previous_accession is None:
        max_rfam_genome_id = db.fetch_max_RG_accession_from_genome()
    else:
        max_rfam_genome_id = previous_accession

    acc_int_section = int(max_rfam_genome_id[2:])

    new_number = acc_int_section + 1

    padded_number = padding_rfam_genome_id_with_zeros(new_number)

    new_rg_accession = "RG" + padded_number

    return new_rg_accession

# ------------------------------------------------------------------------------

def generate_genome_table_entry(accession, previous_rg_acc = None):
    """
    Generates a RfamLive genome table entry based on

    accession: A valid NCBI genome accession
    previous_rg_acc: To be used when generating multiple entries

    return: None
    """

    assembly_level = ["contig", "chromosome", "scaffold", "complete-genome"]
    ncbi_genome_metadata = fetch_genome_metadata_from_ncbi(accession)
    rfam_tax_data = db.fetch_taxonomy_fields()

    level = ''
    # try and guess assembly level
    for level in assembly_level:
        if ncbi_genome_metadata["description"].find(level) != -1:
            break

    kingdom = rfam_tax_data["tax_string"].split(";")[0].strip()

    new_genome_accession = generate_new_rfam_genome_accession(previous_rg_acc)

    genome_table_fields = (new_genome_accession, None, None, None, None, None,
                           level, None, ncbi_genome_metadata["description"],
                           ncbi_genome_metadata["length"],
                           ncbi_genome_metadata["length"], None,
                           ncbi_genome_metadata["ncbi_id"],
                           rfam_tax_data["species"],
                           None,
                           kingdom,
                           0, 0, 0, 0)

    return genome_table_fields

# ------------------------------------------------------------------------------


def parse_arguments():
    """
    Basic function to parse arguments

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser(description="Tool converting NCBI data to RfamDB")

    parser.add_argument("--input",
                        help="This can be a valid NCBI accession or a file with an accession list",
                        type=str)
    mutually_exclusive_args = parser.add_mutually_exclusive_group()
    mutually_exclusive_args.add_argument("--taxid-list",
                        help="Generates a taxid list based on the input accession provided",
                        action='store_true')
    mutually_exclusive_args.add_argument("--genome", help="Generates genome table metadata",
                                         action='store_true')

    return parser

# ------------------------------------------------------------------------------


if __name__ == "__main__":

    parser = parse_arguments()
    args = parser.parse_args()

    # generate a new list of tax ids
    if args.taxid_list:
        if not os.path.isfile(args.input):
            print(fetch_taxid_from_ncbi(args.input))
        else:
            fp = open(args.input, 'r')
            for line in fp:
                line = line.strip()
                tax_id = fetch_taxid_from_ncbi(line)
                print ("%s\t%s" % (line, tax_id))

    elif args.genome:
        # case in which the input is a NCBI accession
        if not os.path.isfile(args.input):
            new_entry = generate_genome_table_entry(args.input, previous_rg_acc=None)
        # this is a list of accessions
        else:
            fp = open(args.input, 'r')
            accessions = [x.strip() for x in fp]
            fp.close()

            previous_rg_acc = None
            entry_list = []
            for accession in accessions:
                new_entry = generate_genome_table_entry(accession, previous_rg_acc=previous_rg_acc)
                entry_list.append(new_entry)
                previous_rg_acc = new_entry[0] # fetch current RG entry because it's not
