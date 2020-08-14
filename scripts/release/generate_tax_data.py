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

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import httplib
import xml.etree.ElementTree as ET
import requests
import argparse

# -----------------------------------------------------------------------------

def get_taxonomy_entries_from_ncbi(tax_ids):
    """
    Parses NCBI's taxonomy xml and generates entries for  taxonomy table
    in rfam_live

    tax_ids: Expects a file with a list of tax ids or a list

    return: void
    """

    # we can move this one to config
    ncbi_cgi_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
    ncbi_tax_query = "db=taxonomy&id=%s&retmode=xml"

    ncbi_tax_url = ncbi_cgi_url + ncbi_tax_query

    tax_id_list = None
    tax_entry_list = []

    # check data and get tax id list
    if isinstance(tax_ids, str):
        if os.path.isfile(tax_ids):
            tax_ids_fp = open(os.path.abspath(tax_ids), 'r')
            tax_id_list = [x.strip() for x in tax_ids_fp]
            tax_ids_fp.close()

        else:
            # maybe create help and point to that one
            sys.exit("Wrong input. Please provide a list of tax ids or a file.")

    elif isinstance(tax_ids, list):
        tax_id_list = tax_ids

    else:
        sys.exit("Wrong input. Please provide a list of tax ids or a file.")

    for tax_id in tax_id_list:

        response = requests.get(ncbi_tax_url % str(tax_id))

        if response.status_code == httplib.OK:

            xml_root = ET.fromstring(response.content)
            taxon_node = xml_root.find('Taxon')

            # add tax id to the list
            tax_entry_list.append(str(tax_id))

            # find scientific name
            scientific_name = taxon_node.find("ScientificName").text

            # species will be scientific name unless common name is found
            species = scientific_name

            # look for common name
            other_names = taxon_node.find("OtherNames")

            if other_names is not None:
                common_name = other_names.find("GenbankCommonName")
                if common_name is not None:
                    common_name = common_name.text
                    species = scientific_name + " (%s)" % common_name
            # if no common name available, species will be scientific name
            else:
                species = scientific_name

            tax_entry_list.append(species)

            # find superkingdom
            superkingdom = ''
            lineage_nodes = taxon_node.find("LineageEx").findall("Taxon")
            for tax_node in lineage_nodes:
                # search for superkingdom and break loop
                if tax_node.find("Rank").text.find("superkingdom") != -1:
                    superkingdom = tax_node.find("ScientificName").text
                    break

            # get tax string
            lineage_str = taxon_node.find("Lineage").text
            lineage_components = lineage_str.split("; ")
            lineage_cleaned = None

            # discard unwanted taxons
            for component in lineage_components:
                if component.find(superkingdom) != -1:
                    sup_index = lineage_components.index(component)
                    lineage_cleaned = lineage_components[sup_index:]
                    break

            tax_string = "; ".join(lineage_cleaned)
            tax_entry_list.append(tax_string)

            # get tree display name
            tree_display_name = species.replace(' ', '_')
            tax_entry_list.append(tree_display_name)

            # get align display name
            align_display_name = tree_display_name + "[%s]" % str(tax_id)
            tax_entry_list.append(align_display_name)

            print ('\t'.join(tax_entry_list))

            tax_entry_list = []

# -----------------------------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing with Python's argparse

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser("Script to generate taxonomy data for genome import")

    parser.add_argument("-f", help="A file containing a list of valid taxids", action="store")

    return parser

# -----------------------------------------------------------------------------


if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    get_taxonomy_entries_from_ncbi(args.f)

