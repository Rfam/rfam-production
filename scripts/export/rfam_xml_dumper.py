"""
Copyright [2009-2016] EMBL-European Bioinformatics Institute
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

"""
Description: This module exports Rfam data

TO DO:
       - Optimizations (motif_xml_dumper, family_xml_dumper, clan_xml_dumper)
       - Set release version and date automatically
"""

# ----------------------------------------------------------------------------

import logging
import timeit
import sys
import os
import datetime
import argparse
import xml.etree.ElementTree as ET
from sets import Set
from xml.dom import minidom
from config import rfam_search as rs
from config import rfam_config as rfc
from utils import RfamDB
from utils.parse_taxbrowser import *


# ----------------------------------------------------------------------------


def xml4db_dumper(name_dict, name_object, entry_type, entry_acc, hfields, outdir):
    """
    Exports query results into EB-eye's XML4dbDUMP format

    name_dict:  A dictionary with all ncbi names per tax id
    name_object: NCBI tax browser node dictionary
    entry_type: Single char signifying the type of the entry accession
                ('M': Motif, 'F': Family, 'C': Clan)
    entry_acc:  An Rfam related accession (Clan, Motif, Family)
    outdir: Destination directory
    """

    entry_type = entry_type[0].capitalize()

    # EB_eye_search fixed tags
    db_xml = ET.Element("database")
    ET.SubElement(db_xml, "name").text = rs.DB_NAME
    ET.SubElement(db_xml, "description").text = rs.DB_DESC

    # need to fetch from db
    ET.SubElement(db_xml, "release").text = rs.DB_RELEASE

    rel_date = datetime.date.today()
    rel_date = rel_date.strftime("%d/%m/%Y")

    ET.SubElement(db_xml, "release_date").text = rel_date

    # need to add a parameter for this...
    ET.SubElement(db_xml, "entry_count").text = '1'

    entries = ET.SubElement(db_xml, "entries")

    # call family xml builder to add a new family to the xml tree
    if (entry_type == rs.FAMILY):
        family_xml_builder(
            name_dict, name_object, entries, rfam_acc=entry_acc, hfields=hfields)

    elif (entry_type == rs.CLAN):
        clan_xml_builder(entries, clan_acc=entry_acc)

    elif (entry_type == rs.MOTIF):
        motif_xml_builder(entries, motif_acc=entry_acc)

    # export xml tree - writes xml tree into a file
    fp_out = open(os.path.join(outdir, entry_acc + ".xml"), 'w')

    db_str = ET.tostring(db_xml, "utf-8")
    db_str_reformated = minidom.parseString(db_str)

    fp_out.write(db_str_reformated.toprettyxml(indent='\t'))

    fp_out.close()

    # ----------------------------------------------------------------------------


def family_xml_builder(name_dict, name_object, entries, rfam_acc=None, hfields=True):
    """
    Expands the Xml4dbDumper object by adding a new family entry

    name_dict:  A dictionary with all ncbi names per tax id
    name_object: NCBI tax browser node dictionary
    entries:    The xml entries node to be expanded
    rfam_acc:   A specific Rfam family accession
    hfields:    A bool value indicating whether to build hierarchical fields
                for each species. True by default
    """

    entry_type = "Family"

    cross_refs = {}

    # fetch family fields
    fam_fields = fetch_entry_fields(rfam_acc, rs.FAMILY)

    # fetch family specific ncbi_ids
    ncbi_ids = fetch_value_list(rfam_acc, rs.NCBI_IDs_QUERY)
    valid_ncbi_ids = get_valid_family_tax_ids(name_object, ncbi_ids)

    # fetch family specific ncbi_ids
    pdb_ids = fetch_value_list(rfam_acc, rs.PDB_IDs_QUERY)

    # get pubmed ids
    pmids = get_value_list(fam_fields["pmids"], ',')

    # get dbxrefs (GO, SO)
    dbxrefs = get_value_list(fam_fields["dbxrefs"], ',')

    # get associated clan
    clan = fetch_value(rs.FAM_CLAN, rfam_acc)

    # need a function here to split dbxrefs in a pretty way
    go_ids = filter(lambda x: x.find("GO") != -1, dbxrefs)
    so_ids = filter(lambda x: x.find("SO") != -1, dbxrefs)

    # update cross references dictionary
    cross_refs["ncbi_taxonomy_id"] = valid_ncbi_ids
    cross_refs["PDB"] = pdb_ids
    cross_refs["PUBMED"] = pmids
    cross_refs["GO"] = go_ids
    cross_refs["SO"] = so_ids

    if clan is not None:
        cross_refs["RFAM"] = [clan]

    # add a new family entry to the xml tree
    entry = ET.SubElement(entries, "entry", id=rfam_acc)

    # entry name
    ET.SubElement(entry, "name").text = str(fam_fields["name"])

    # entry description
    ET.SubElement(entry, "description").text = str(fam_fields["description"])

    # entry dates - common to motifs and clans
    dates = ET.SubElement(entry, "dates")

    created = fam_fields["created"].date().strftime("%d %b %Y")
    updated = fam_fields["updated"].date().strftime("%d %b %Y")

    ET.SubElement(dates, "date", value=created, type="created")
    ET.SubElement(dates, "date", value=updated, type="updated")

    # loop to add cross references
    cross_refs = build_cross_references(entry, cross_refs)

    # expand xml tree with additional fields
    add_fields = build_additional_fields(
        entry, fam_fields, len(pdb_ids), valid_ncbi_ids, entry_type=entry_type)

    if hfields is True:
        species_tax_trees = get_family_tax_tree(
            name_object, name_dict, valid_ncbi_ids)
        if len(species_tax_trees.keys()) > 1:
            add_hierarchical_fields(add_fields, species_tax_trees, name_dict)


# ----------------------------------------------------------------------------


def clan_xml_builder(entries, clan_acc=None):
    """
    Expands the Xml4dbDumper object by adding a new clan entry

    entries:    The xml entries node to be expanded
    clan_acc:   An Rfam associated clan accession
    """

    entry_type = "Clan"

    cross_ref_dict = {}

    # fetch clan fields
    clan_fields = fetch_entry_fields(clan_acc, rs.CLAN)

    # add a new clan entry to the xml tree
    entry = ET.SubElement(entries, "entry", id=clan_acc)

    ET.SubElement(entry, "name").text = clan_fields["name"]
    ET.SubElement(entry, "description").text = clan_fields["description"]

    # entry dates - common to motifs and clans
    dates = ET.SubElement(entry, "dates")

    created = clan_fields["created"].date().strftime("%d %b %Y")
    updated = clan_fields["updated"].date().strftime("%d %b %Y")

    ET.SubElement(dates, "date", value=created, type="created")
    ET.SubElement(dates, "date", value=updated, type="updated")

    # clan cross references
    clan_fams = fetch_value_list(clan_acc, rs.CLAN_FAMS)
    cross_ref_dict["RFAM"] = clan_fams
    build_cross_references(entry, cross_ref_dict)

    # clan additional fields
    build_additional_fields(
        entry, clan_fields, clan_fields["num_families"], None, entry_type=entry_type)


# ----------------------------------------------------------------------------


def motif_xml_builder(entries, motif_acc=None):
    """
    Expands the Xml4dbDump with a Motif entry

    entries:    Entries node on xml tree
    motif_acc:  An Rfam associated motif accession
    """

    entry_type = "Motif"

    cross_ref_dict = {}

    # fetch clan fields
    motif_fields = fetch_entry_fields(motif_acc, rs.MOTIF)

    # add a new clan entry to the xml tree
    entry = ET.SubElement(entries, "entry", id=motif_acc)

    ET.SubElement(entry, "name").text = motif_fields["name"]
    ET.SubElement(entry, "description").text = motif_fields["description"]

    # entry dates - common to motifs and clans
    dates = ET.SubElement(entry, "dates")

    created = motif_fields["created"].date().strftime("%d %b %Y")
    updated = motif_fields["updated"].date().strftime("%d %b %Y")

    ET.SubElement(dates, "date", value=created, type="created")
    ET.SubElement(dates, "date", value=updated, type="updated")

    # adding cross references
    motif_fams = fetch_value_list(motif_acc, rs.MOTIF_FAMS)
    cross_ref_dict["RFAM"] = motif_fams
    build_cross_references(entry, cross_ref_dict)

    build_additional_fields(
        entry, motif_fields, 0, None, entry_type=entry_type)


# ----------------------------------------------------------------------------


def build_cross_references(entry, cross_ref_dict):
    """
    Expands the entry xml tree by adding the entry's cross references.
    Returns the cross references xml tree

    entry:  The entry node of the xml tree object (xml.etree.ElementTree)
    cross_ref_dict: A dictionary with the entity's cross references in the
                    form of ({db_name:[db_key1,db_key2,..],}) where db_name is
                    a string and values a list of db ids
    """

    # cross_ref_dict will be different for the different types,
    # but the dictionary has to be in the same format

    cross_refs = ET.SubElement(entry, "cross_references")

    for db_name in cross_ref_dict.keys():

        # get db_keys
        db_keys = cross_ref_dict[db_name]
        if len(db_keys) > 0:
            for value in db_keys:
                ET.SubElement(
                    cross_refs, "ref", dbkey=str(value), dbname=db_name)

    return cross_refs


# ----------------------------------------------------------------------------


def add_hierarchical_fields(xml_tree_node, tax_tree_dict, name_dict):
    """
    Expands the cross references xml tree by adding hierarchical references
    for the ncbi ids in valid_ncbi_ids.

    xml_tree_node:  An existing xml tree node to expand with hierarchical
                    fields
    tax_tree_dict:  Species taxonomy tree dictionary as generated by
                    get_family_tax_tree
    name_dict:  NCBI's name dictionary as returned by read_ncbi_names_dmp
    """

    # add a new hierarchical ref for every tax_id in the family
    for tax_id in tax_tree_dict.keys():
        hfields = ET.SubElement(xml_tree_node, "hierarchical_field",
                                name="taxonomy_lineage")

        # fetch lineage
        lineage = tax_tree_dict[tax_id]
        tax_tree = lineage[::-1]

        for tax_tree_node in tax_tree:
            # create the root node
            if tax_tree_node == '1':
                # need to create one root node in order for the xml dump to
                # validate
                ET.SubElement(hfields, "root", label="root").text = '1'
            else:
                # skip root while creating child nodes
                if tax_tree_node != '1':
                    ET.SubElement(
                        hfields, "child", label=name_dict[tax_tree_node]).text = tax_tree_node


# ----------------------------------------------------------------------------


def build_additional_fields(entry, fields, num_3d_structures, fam_ncbi_ids, entry_type):
    """
    This function expands the entry xml field with the additional fields

    entry:  This is the xml.etree.ElementTree at the point of entry
    fields: A list of additional fields to expand the entry with
    """

    add_fields = ET.SubElement(entry, "additional_fields")

    # adding entry type
    ET.SubElement(add_fields, "field", name="entry_type").text = entry_type

    # adding authors
    authors = fields["author"]
    authors = authors.replace(';', ',')
    author_list = get_value_list(authors, rs.AUTH_DEL)

    # to be deleted when author name is corrected on the db
    if author_list.count("Argasinska") > 0:
        author_list.remove("J")
        author_list.remove("Argasinska")
        author_list.append("Argasinska J")

    for author in author_list:
        ET.SubElement(add_fields, "field", name="author").text = author

    if entry_type == "Family":

        # number of species
        ET.SubElement(add_fields, "field", name="num_species").text = str(fields[
                                                                              "num_species"])
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
        rna_types = get_value_list(fields["rna_type"], rs.RNA_TYPE_DEL)

        for rna_type in rna_types:
            ET.SubElement(add_fields, "field", name="rna_type").text = rna_type

        # has 3d structure
        if num_3d_structures > 0:
            ET.SubElement(
                add_fields, "field", name="has_3d_structure").text = "Yes"
        else:
            ET.SubElement(
                add_fields, "field", name="has_3d_structure").text = "No"

        # add popular species if any
        for species in rs.POPULAR_SPECIES:
            if species in fam_ncbi_ids:
                ET.SubElement(
                    add_fields, "field", name="popular_species").text = str(species)

                # build hierarchical_fields tree here...

    # perhaps move this to clan and motif xml builder
    else:
        num_families = None

        if entry_type == "Motif":
            num_families = fetch_value(
                rs.NUM_FAMS_MOTIF, fields["id"])

        elif entry_type == "Clan":
            num_families = fields["num_families"]

        ET.SubElement(add_fields, "field", name="num_families").text = str(
            num_families)

    # returning node
    return add_fields


# ----------------------------------------------------------------------------


def get_value_list(val_str, delimiter=','):
    """
    Splits an input string according to delimiter and returns a list of the
    elements

    val_str:    A string of family specific values. This string is a
                concatenation of multiple values related to a single family
    delimiter:  The delimeter that will be used to split the values' string
    """

    val_str = val_str.strip()
    values = val_str.split(delimiter)

    value_list = map(lambda x: x.strip(), values)
    value_list = filter(lambda x: x != '', value_list)

    return value_list


# ----------------------------------------------------------------------------

def fetch_value_list(rfam_acc, query):
    """
    Retrieves and returns a list of all rfam_acc related values, returned
    by executing the query. Values in list are converted to string format.
    If rfam_acc is None then query is executed without an rfam_acc

    rfam_acc:   A family specific accession
    query:  A string with the MySQL query to be executed
    """

    cnx = RfamDB.connect()

    cursor = cnx.cursor(raw=True)

    if rfam_acc is None:
        cursor.execute(query)

    else:
        cursor.execute(query % rfam_acc)

    values = cursor.fetchall()

    cursor.close()
    cnx.disconnect()

    return map(lambda x: str(x[0]), values)


# ----------------------------------------------------------------------------


def fetch_entry_fields(entry_acc, entry_type):
    """
    Returns a dictionary with the entry's fields.

    entry_acc:  An Rfam associated accession (Motif, Clan, Family)
    entry_type: The type of the entry accession
    """

    # maybe the entry type not required... use rfam_acc[0:2]

    cnx = RfamDB.connect()
    cursor = cnx.cursor(dictionary=True)

    entry_type = entry_type[0].capitalize()

    try:
        if entry_type == rs.FAMILY:
            cursor.execute(rs.FAM_FIELDS % entry_acc)

        elif entry_type == rs.CLAN:
            cursor.execute(rs.CLAN_FIELDS % entry_acc)

        elif entry_type == rs.MOTIF:
            cursor.execute(rs.MOTIF_FIELDS % entry_acc)

        fields = cursor.fetchall()[0]

    except:
        print "Failure retrieving values for entry %s." % entry_acc

    cursor.close()
    cnx.disconnect()

    return fields


# ----------------------------------------------------------------------------


def fetch_value(query, accession):
    """
    Retrieves and returns a value from the database depending to the query
    executed. The query should return a single value

    query:  The query to be executed in the form of string.
    accession:  Rfam specific accession (family, clan, motif)
                to execute the query on
    """

    cnx = RfamDB.connect()

    cursor = cnx.cursor(raw=True)

    cursor.execute(query % accession)

    value = cursor.fetchall()

    cursor.close()
    cnx.disconnect()

    if len(value) > 0:
        return value[0][0]

    return None


# ----------------------------------------------------------------------------


def main(entry_type, rfam_acc, outdir, hfields=True):
    """
    This function puts everything together.

    entry_type: One of the three entry types in Rfam (Motif, Clan, Family)
    rfam_acc: An Rfam associated accession (RF*,CL*,RM*). If rfam_acc is set
              to None, then all data related to the entry type will be
              exported
    hfields: A flag (True/False) indicating whether to add hierarchical
             fields on not. True by default
    outdir: Destination directory
    """

    rfam_accs = None
    entry = ""

    name_object = {}
    name_dict = {}

    try:

        # this will create the output directory if it doesn't exist
        if not os.path.exists(outdir):

            try:
                os.mkdir(outdir)
            except:
                print "Error creating output directory at: ", outdir

        # export all entries
        if rfam_acc is None:

            # Motif accessions
            if entry_type == rs.MOTIF:
                rfam_accs = fetch_value_list(
                    None, rs.MOTIF_ACC)

            # Clan accessions
            elif entry_type == rs.CLAN:
                rfam_accs = fetch_value_list(
                    None, rs.CLAN_ACC)

            # Family accessions
            elif entry_type == rs.FAMILY:

                # load ncbi taxonomy browser here
                name_dict, name_dict_reverse = read_ncbi_names_dmp(
                    rfc.TAX_NAMES_DUMP)
                name_object = read_ncbi_taxonomy_nodes(
                    name_dict, rfc.TAX_NODES_DUMP)

                rfam_accs = fetch_value_list(
                    None, rs.FAM_ACC)

                for entry in rfam_accs:
                    t0 = timeit.default_timer()
                    xml4db_dumper(
                        name_dict, name_object, entry_type, entry, hfields, outdir)
                    print "Execution time for %s: %s" % (entry, str(timeit.default_timer() - t0))

                return

            # Don't build hierarchical references for Clans and Motifs
            for entry in rfam_accs:
                t0 = timeit.default_timer()
                xml4db_dumper(
                    None, None, entry_type, entry, False, outdir)
                print "Execution time for %s: %s" % (entry, str(timeit.default_timer() - t0))

        # export single entry
        else:
            # need to check the validity of an rfam_acc (rfam, motif, clan)
            if entry_type == rs.MOTIF or entry_type == rs.CLAN:
                xml4db_dumper(None, None, entry_type, rfam_acc, False, outdir)

            # export single family entry
            else:

                # load ncbi taxonomy browser here
                name_dict, name_dict_reverse = read_ncbi_names_dmp(
                    rfc.TAX_NAMES_DUMP)
                name_object = read_ncbi_taxonomy_nodes(
                    name_dict, rfc.TAX_NODES_DUMP)

                xml4db_dumper(
                    name_dict, name_object, entry_type, rfam_acc, hfields, outdir)

    except:
        # need to correct this one
        if rfam_acc is None:
            gen_fams = Set([x.partition('.')[0] for x in os.listdir(outdir)])
            loaded_fams = Set(rfam_accs)
            # get remaining families
            rem_fams = loaded_fams - gen_fams

            # open a log file
            logging.basicConfig(
                filename=os.path.join("missing_accs" + ".log"), filemode='w', level=logging.DEBUG)

            # write accessions to log file
            for rfam_acc in rem_fams:
                logging.debug(rfam_acc)
        else:
            print "Error exporting %s." % rfam_acc


# ----------------------------------------------------------------------------


def get_valid_family_tax_ids(name_object, family_tax_ids):
    """
    Returns a list of all family tax ids found in the NCBI dumps

    name_object: NCBI tax browser node dictionary
    family_tax_ids: A list of all family ncbi ids
    """

    valid_family_tax_ids = []

    for taxid in family_tax_ids:
        if (taxid in name_object):
            valid_family_tax_ids.append(taxid)

    return valid_family_tax_ids


# ----------------------------------------------------------------------------


def get_family_tax_tree(name_object, name_dict, family_tax_ids):
    """
    Returns the family genealogy list

    name_object: NCBI tax browser node dictionary
    name_dict: A dictionary with all ncbi names per tax id
    family_tax_ids: A list of family specific ncbi ids
    """

    species_tax_trees = {}

    for taxid in family_tax_ids:

        if (taxid in name_object):
            species_tax_trees[taxid] = name_object[
                taxid].get_lineage(name_object)

    return species_tax_trees


# ----------------------------------------------------------------------------

def usage():
    """
    Parses arguments and displays usage information on screen
    """

    parser = argparse.ArgumentParser(
        description="Rfam Search Xml4db Dumper.", epilog='')

    # group required arguments together
    req_args = parser.add_argument_group("required arguments")

    req_args.add_argument("--type", help="rfam entry type (F: Family, M: Motif, C: Clan)",
                          type=str, choices=['F', 'M', 'C'], required=True)

    parser.add_argument(
        "--acc", help="a valid rfam entry accession (RF*|CL*|RM*)",
        type=str, default=None)

    parser.add_argument(
        "--hfields", help="include hierarchical fields", action="store_true")

    req_args.add_argument(
        "--out", help="path to output directory", type=str, required=True)

    return parser


# ----------------------------------------------------------------------------

if __name__ == '__main__':

    parser = usage()
    args = parser.parse_args()

    # Additional checks
    # Check if export type matches accession
    wrong_input = False
    if (args.acc is not None):

        if (args.type == 'F' and args.acc[0:2] != "RF"):
            wrong_input = True
        elif (args.type == 'M' and args.acc[0:2] != "RM"):
            wrong_input = True
        elif (args.type == 'C' and args.acc[0:2] != "CL"):
            wrong_input = True

    if (wrong_input is True):
        print "\nAccession does not match the export type.\n"
        parser.print_help()
        sys.exit()

    # check output directory
    if (os.path.isdir(args.out) is False):
        print "\nPlease provide a valid output directory.\n"
        parser.print_help()
        sys.exit()

    main(args.type, args.acc, args.out, hfields=args.hfields)
