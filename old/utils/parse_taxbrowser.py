#!/usr/bin/python

"""
Based on:
http://evosite3d.blogspot.co.uk/2013/06/browsing-ncbi-taxonomy-with-python.html

The data is available at:
ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

Note: This module has been modified in order for it to work with rfam_xml_dumper

Modifications: - name_object passed as param to Node's method get_linage
               - name_dict passed as param to read_ncbi_taxonomy_nodes
"""

# ----------------------------------------------------------------------------

import os
import sys

# ----------------------------------------------------------------------------


class Node:
    """
    NCBI Taxonomy node.
    """

    def __init__(self):
        self.tax_id = 0  # Number of the tax id.
        self.parent = 0  # Number of the parent of this node
        self.children = []  # List of the children of this node
        self.tip = 0  # Tip=1 if it's a terminal node, 0 if not.
        # Name of the node: taxa if it's a terminal node, numero if not.
        self.name = ""
        # 1 if name is suppressed in GenBank entry lineage
        self.hidden = 0

    def get_lineage(self, name_object):
        """
        Trace get_lineage from root to leaf
        """
        ancestors = []  # Initialise the list of all nodes from root to leaf.
        tax_id = self.tax_id  # Define leaf

        while 1:
            if name_object.has_key(tax_id):
                ancestors.append(tax_id)
                tax_id = name_object[tax_id].parent
            else:
                break
            if tax_id == "1":
                # If it is the root, we reached the end.
                # Add it to the list and break the loop
                ancestors.append(tax_id)
                break
        return ancestors  # Return the list


# ----------------------------------------------------------------------------


def read_ncbi_names_dmp(names_dmp):
    """
    Load  NCBI names file ("names.dmp")

    names_dmp: NCBIs names.dmp file
    """
    name_dict = {}  # Initialise dictionary with TAX_ID:NAME
    name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID

    name_file = open(names_dmp, "r")
    while 1:
        line = name_file.readline()
        if line == "":
            break
        line = line.rstrip()
        line = line.replace("\t", "")
        tab = line.split("|")
        if tab[3] == "scientific name":
            tax_id, name = tab[0], tab[1]  # Assign tax_id and name ...
            name_dict[tax_id] = name  # ... and load them
            name_dict_reverse[name] = tax_id  # ... into dictionaries
    name_file.close()

    return name_dict, name_dict_reverse


# ----------------------------------------------------------------------------


def read_ncbi_taxonomy_nodes(name_dict, nodes_dmp):
    """
    Load NCBI taxonomy nodes file ("nodes.dmp")

    name_dict: NCBI names dict returned by read_ncbi_names_dmp
    nodes_dmp: NCBIs nodes.dmp file
    """
    name_object = {}
    taxonomy_file = open(nodes_dmp, "r")

    while 1:
        line = taxonomy_file.readline()
        if line == "":
            break
        line = line.replace("\t", "")
        tab = line.split("|")

        tax_id = str(tab[0])
        tax_id_parent = str(tab[1])
        division = str(tab[4])
        hidden = tab[10]

        # Define name of the taxid
        name = "unknown"

        if tax_id in name_dict:
            name = name_dict[tax_id]

        if tax_id not in name_object.keys():
            name_object[tax_id] = Node()
        name_object[tax_id].tax_id = tax_id  # Assign tax_id
        name_object[tax_id].parent = tax_id_parent  # Assign tax_id parent
        name_object[tax_id].name = name  # Assign name
        name_object[tax_id].hidden = hidden

        if tax_id_parent in name_object:
            # If parent is is already in the object
            children = name_object[tax_id].children
            # ... we found its children.
            children.append(tax_id)
            # ... so add them to the parent
            name_object[tax_id].children = children
    taxonomy_file.close()

    return name_object


# -----------------------------------------------------------------------------
