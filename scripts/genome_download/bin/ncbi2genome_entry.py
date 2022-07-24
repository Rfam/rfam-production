#!/usr/bin/env python3

import xml.etree.ElementTree as ET
import typing as ty
import datetime


import requests
import click

# from utils import db_utils as db


def fetch_genome_metadata_from_ncbi(accession: str) -> ty.Dict[str, str]:
    """
    Fetches genome related metadata from NCBI using a
    valid sequence accession.

    accession: A valid NCBI genome accession

    return: A python dictionary with genome related metadata
    """

    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={accession}"
    response = requests.get(url)
    response.raise_for_status()
    text = response.text

    root = ET.fromstring(text)
    if root is None:
        raise ValueError("Did not fetch any metadata for %s" % accession)

    docsum = root.find("DocSum")
    if docsum is None:
        raise ValueError("Invalid xml, did not find DocSum")

    items = docsum.findall("Item")
    if not items:
        raise ValueError("Invalid xml, did not find any Item entries")

    genome_metadata = {}
    for item in items:
        if item.get("Name") == "Title":
            genome_metadata["description"] = item.text
        if item.get("Name") == "TaxId":
            genome_metadata["ncbi_id"] = item.text
        if item.get("Name") == "Length":
            genome_metadata["length"] = item.text

    return genome_metadata


def generate_new_rfam_genome_accession(previous_accession=None) -> str:
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

    return "RG%09i" % new_number


def generate_genome_table_entry(accession: str, previous_rg_acc=None) -> ty.List[str]:
    """
    Generates a RfamLive genome table entry based on. This requires that
    taxonomy table is populated before the genome table.

    accession: A valid NCBI genome accession
    previous_rg_acc: To be used when generating multiple entries

    return: None
    """

    assembly_level = ["contig", "chromosome", "scaffold", "complete-genome"]
    ncbi_genome_metadata = fetch_genome_metadata_from_ncbi(accession)
    rfam_tax_data = db.fetch_taxonomy_fields(ncbi_genome_metadata["ncbi_id"])

    level = ""
    for level in assembly_level:
        if ncbi_genome_metadata["description"].find(level) != -1:
            break
    else:
        raise ValueError("Failed to find the assembly level for %s" % accession)

    kingdom = rfam_tax_data["tax_string"].split(";")[0].strip().lower()

    if kingdom != "viruses":
        kingdom = "viruses"

    new_genome_accession = generate_new_rfam_genome_accession(previous_rg_acc)
    version = int(accession.partition(".")[2])

    created_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    return [
        new_genome_accession,
        accession,
        str(version),
        r"\N",
        r"\N",
        r"\N",
        level,
        r"\N",
        str(ncbi_genome_metadata["description"]),
        str(ncbi_genome_metadata["length"]),
        str(ncbi_genome_metadata["length"]),
        r"\N",
        str(ncbi_genome_metadata["ncbi_id"]),
        str(rfam_tax_data["species"]),
        r"\N",
        str(kingdom),
        "0",
        "0",
        "0",
        "0",
        str(created_date),
    ]


# def parse_arguments():
#     """
#     Basic function to parse arguments
#     return: Argparse parser object
#     """
#     parser = argparse.ArgumentParser(description="Tool converting NCBI data to RfamDB")
#     parser.add_argument("--input",
#                         help="This can be a valid NCBI accession or a file with an accession list",
#                         type=str)
#     mutually_exclusive_args = parser.add_mutually_exclusive_group()
#     mutually_exclusive_args.add_argument("--taxid-list",
#                         help="Generates a taxid list based on the input accession provided",
#                         action='store_true')
#     mutually_exclusive_args.add_argument("--genome", help="Generates genome table metadata",
#                                          action='store_true')
#     return parser


@click.command()
@click.argument()
def main():
    pass


# if __name__ == "__main__":
#     main()

#     parser = parse_arguments()
#     args = parser.parse_args()

#     # generate a new list of tax ids
#     if args.taxid_list:
#         if not os.path.isfile(args.input):
#             print(fetch_taxid_from_ncbi(args.input))
#         else:
#             fp = open(args.input, 'r')
#             for line in fp:
#                 line = line.strip()
#                 tax_id = fetch_taxid_from_ncbi(line)
#                 print ("%s\t%s" % (line, tax_id))

#     elif args.genome:
#         # case in which the input is a NCBI accession
#         if not os.path.isfile(args.input):
#             new_entry = generate_genome_table_entry(args.input, previous_rg_acc=None)

#         # this works on a list of accessions
#         else:
#             fp = open(args.input, 'r')
#             accessions = [x.strip() for x in fp]
#             fp.close()

#             previous_rg_acc = None
#             entry_list = []

#             for accession in accessions:
#                 new_entry = generate_genome_table_entry(accession, previous_rg_acc=previous_rg_acc)
#                 entry_list.append(new_entry)
#                 print '\t'.join(list(new_entry))
#                 # when populating genome table with multiple entries we need to keep track of
#                 # the current/previous RG accessions to generate unique UPIDs for new genomes
#                 previous_rg_acc = new_entry[0]
