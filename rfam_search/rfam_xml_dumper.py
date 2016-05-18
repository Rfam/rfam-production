'''
Created on 13 May 2016

@author: ikalvari

Description: This module exports Rfam data

TO DO: Need to add functions to parse and update an Xml4dbDumper file and 
       update entries (family, clan, motif)
'''

# ----------------------------------------------------------------------------

import os
import sys
import xml.etree.ElementTree as ET
import datetime
from xml.dom import minidom
from config import rfam_search as rs_conf
from utils import RfamDB  # will load those from rfam-public

# ----------------------------------------------------------------------------

# this should be able to dump a single family, clan, motif or the entire db..


def xml4db_dumper(outdir):
    '''
        exports query results into EB-eye's XML4dbDUMP format
        outdir: Destination directory

        Maybe here provide the fields as a txt file and the dump in txt format
        and according to that dump the xml file
    '''

    # EB_eye_search fixed tags
    db_xml = ET.Element("database")
    ET.SubElement(db_xml, "name").text = rs_conf.DB_NAME
    ET.SubElement(db_xml, "description").text = rs_conf.DB_DESC
    ET.SubElement(db_xml, "release").text = rs_conf.DB_RELEASE
    ET.SubElement(db_xml, "release_date").text = rs_conf.DB_REL_DATE

    entries = ET.SubElement(db_xml, "entries")

    # call family xml builder to add a new family to the xml tree
    family_xml_builder(entries, rfam_acc='RF00177')

    tree = ET.ElementTree(db_xml)

    fp_out = open(os.path.join(outdir, "xml4dbsample.xml"), 'w')

    db_str = ET.tostring(db_xml, 'utf-8')
    db_str_reformated = minidom.parseString(db_str)

    fp_out.write(db_str_reformated.toprettyxml(indent="\t"))

    fp_out.close()

# ----------------------------------------------------------------------------


def family_xml_builder(entries, rfam_acc=None):
    '''
        Expands the Xml4dbDumper object by adding a new family entry.

        entries: The xml entries node to be expanded
        rfam_acc: A specific Rfam family accession
    '''

    cnx = RfamDB.connect()
    cursor = cnx.cursor(dictionary=True)

    # fetch family specific ncbi_ids
    ncbi_ids = fetch_value_list(rfam_acc, rs_conf.NCBI_IDs_QUERY)

    # fetch family specific ncbi_ids
    pdb_ids = fetch_value_list(rfam_acc, rs_conf.PDB_IDs_QUERY)

    # fetch family fields
    cursor.execute(rs_conf.FAM_QUERY % rfam_acc)
    fam_fields = cursor.fetchall()[0]

    # add a new family entry to the xml tree
    entry = ET.SubElement(entries, "entry", id=rfam_acc)

    # entry name
    ET.SubElement(entry, "name").text = str(fam_fields["rfam_id"])

    # entry description
    ET.SubElement(entry, "description").text = str(fam_fields["description"])

    # entry dates - common to motifs and clans
    dates = ET.SubElement(entry, "dates")

    created = fam_fields["created"].date().strftime("%d %b %Y")
    updated = fam_fields["updated"].date().strftime("%d %b %Y")

    ET.SubElement(dates, "date", value=created, type="created")
    ET.SubElement(dates, "date", value=updated, type="updated")

    # loop to add cross references

    # expand xml tree with additional fields
    add_additional_fields(entry, fam_fields, len(pdb_ids), entry_type='Family')

    # CANN function to fetch_pdb_ids

    cursor.close()
    cnx.disconnect()

# ----------------------------------------------------------------------------


def expand_xml_tree(xml_tree_node, value_list):
    '''
        Expands an xml tree from a point onwards by adding new fields.

        TO BE IMPLEMENTED
    '''

    pass
# ----------------------------------------------------------------------------


def add_cross_references(entry):
    # motif cross references

    pass
# ----------------------------------------------------------------------------

# perhaps fields can be a dictionary and move fields to type specific builder??


def add_additional_fields(entry, fields, num_3d_structures, entry_type):
    '''
        This function expands the entry xml field with the additional fields

        entry: This is the xml.etree.ElementTree at the point of entry
        fields: A list of additional fields to expand the entry with

    '''
    add_fields = ET.SubElement(entry, "additional_fields")

    # adding entry type
    ET.SubElement(add_fields, "field", name="entry_type").text = entry_type

    # adding authors
    author_list = get_value_list(fields["author"], rs_conf.AUTH_DEL)

    for author in author_list:
        ET.SubElement(add_fields, "field", name="author").text = author

    if entry_type == "Family":

        # number of species
        ET.SubElement(add_fields, "field", name="num_species").text = str(fields[
            "number_of_species"])
        # number of 3D structures
        ET.SubElement(
            add_fields, "field", name="num_3d_structures").text = str(num_3d_structures)
        # num seed
        ET.SubElement(add_fields, "field", name="num_seed").text = str(fields[
            "num_seed"])
        # num full
        ET.SubElement(add_fields, "field", name="num_full").text = str(fields[
            "num_full"])

        # rna types
        rna_types = get_value_list(fields["type"], rs_conf.RNA_TYPE_DEL)

        for rna_type in rna_types:
            ET.SubElement(add_fields, "field", name="rna_type").text = rna_type

        # has 3d structure
        if num_3d_structures > 0:
            ET.SubElement(
                add_fields, "field", name="has_3d_structure").text = "Yes"
        else:
            ET.SubElement(
                add_fields, "field", name="has_3d_structure").text = "No"

    # perhaps move this to clan and motif xml builder
    else:
        num_families = None

        if entry_type == "Motif":
            num_families = fetch_value(
                rs_conf.NUM_FAMS_MOTIF, fields["motif_acc"])

        elif entry_type == "Clan":
            num_families = fetch_value(
                rs_conf.NUM_FAMS_CLAN, fields["clan_acc"])

        ET.SubElement(add_fields, "field", name="num_families").text = str(
            num_families)

# ----------------------------------------------------------------------------


def motif_xml_builder():
    '''
        TO BE IMPLEMENTED
    '''

    pass

# ----------------------------------------------------------------------------


def clan_xml_builder():
    '''
        TO BE IMPLEMENTED
    '''

    pass

# ----------------------------------------------------------------------------


def get_value_list(val_str, delimeter=','):  # done
    '''
        val_str: A string of family specific values. This string is a
        concatenation of multiple values related to a single family
        delimeter: The delimeter that will be used to split the values' string
    '''

    val_str = val_str.strip()
    values = val_str.split(delimeter)

    value_list = map(lambda x: x.strip(), values)
    value_list = filter(lambda x: x != '', value_list)

    return value_list


# ----------------------------------------------------------------------------

# maybe move this to DB utils


def fetch_value_list(rfam_acc, query):
    '''
        Retrieves and returns a list of all rfam_acc related values, returned
        by executing the query. Values in list are converted to string format.

        rfam_acc: A family specific accession
        query: A string with the MySQL query to be executed
    '''

    cnx = RfamDB.connect()

    cursor = cnx.cursor(raw=True)

    cursor.execute(query % rfam_acc)

    pdb_structures = cursor.fetchall()

    cursor.close()
    cnx.disconnect()

    return map(lambda x: str(x[0]), pdb_structures)

# ----------------------------------------------------------------------------


def fetch_value(query, accession):
    '''
        Retrieves and returns a value from the database depending to the query
        executed.

        query: The query to be executed in the form of string.
        accession: Rfam specific accession (family, clan, motif)
                   to execute the query on.

        *Query should return a single value
    '''

    cnx = RfamDB.connect()

    cursor = cnx.cursor(raw=True)

    cursor.execute(query % accession)

    value = cursor.fetchall()

    cursor.close()
    cnx.disconnect()

    return value[0][0]

# ----------------------------------------------------------------------------

if __name__ == '__main__':
    outdir = "/Users/ikalvari/Desktop/EB_eye_search"
    xml4db_dumper(outdir)
