#!/usr/bin/env python2.7
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
from __future__ import print_function

import argparse
import datetime
import logging
import subprocess
import timeit
import traceback
import xml.etree.ElementTree as ET
from xml.dom import minidom

import django
from config.rfam_config import RFAMLIVE, RFAMREL
from rfam_schemas.RfamLive.models import Genome, Genseq
from sets import Set
from utils import RfamDB
from utils.parse_taxbrowser import *

from config import rfam_config as rfc
from config import rfam_search as rs

"""
Description: This module exports Rfam data for the search engine

TO DO:
       - Optimizations (motif_xml_dumper, family_xml_dumper, clan_xml_dumper)
       - Set release version and date automatically
"""

django.setup()
# settings.configure()

DB_CONFIG = None


def xml4db_dumper(
    name_dict, name_object, entry_type, entry_acc, hfields: bool, outdir: str
):
    """
    Exports query results into EB-eye's XML4dbDUMP format

    name_dict:  A dictionary with all ncbi names per tax id
    name_object: NCBI tax browser node dictionary
    entry_type: Single char signifying the type of the entry accession
                ('M': Motif, 'F': Family, 'C': Clan, 'G': Genome)
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

    entries = ET.SubElement(db_xml, "entries")

    # call family xml builder to add a new family to the xml tree
    if entry_type == rs.FAMILY:
        family_xml_builder(
            name_dict, name_object, entries, rfam_acc=entry_acc, hfields=hfields
        )

    elif entry_type == rs.CLAN:
        clan_xml_builder(entries, clan_acc=entry_acc)

    elif entry_type == rs.MOTIF:
        motif_xml_builder(entries, motif_acc=entry_acc)

    elif entry_type == rs.GENOME:
        genome_xml_builder(entries, gen_acc=entry_acc)

    elif entry_type == rs.MATCH:
        full_region_xml_builder(entries, entry_acc)

    # adding entry_count
    entry_count = len(entries.findall("entry"))

    if entry_count == 0:
        print("No full region entries found for %s" % entry_acc)
        return

    ET.SubElement(db_xml, "entry_count").text = str(entry_count)

    # export xml tree - writes xml tree into a file
    filename = os.path.join(outdir, entry_acc + ".xml")
    if not os.path.exists(filename):
        fp_out = open(filename, "w")

        db_str = ET.tostring(db_xml, "utf-8")
        db_str_reformated = minidom.parseString(db_str)

        fp_out.write(db_str_reformated.toprettyxml(indent="\t"))

        fp_out.close()
    else:
        print("File already exists: ", filename)
    # xmllint(filename)


# ----------------------------------------------------------------------------


def get_taxonomy_info(rfam_acc):
    """
    Get distinct ncbi_ids and tax_strings associated with a family.
    """
    cnx = RfamDB.connect(db_config=DB_CONFIG)
    cursor = cnx.cursor(dictionary=True, buffered=True)
    cursor.execute(rs.NCBI_IDs_QUERY % rfam_acc)
    ncbi_ids = []
    tax_strings = set()  # distinct ncbi_ids can have identical tax_strings
    for row in result_iterator(cursor):
        ncbi_ids.append(row["ncbi_id"])
        tax_strings.add(row["tax_string"])
    cursor.close()
    cnx.disconnect()
    return ncbi_ids, tax_strings


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
    ncbi_ids, tax_strings = get_taxonomy_info(rfam_acc)

    if hfields:
        valid_ncbi_ids = get_valid_family_tax_ids(name_object, ncbi_ids)
    else:
        valid_ncbi_ids = ncbi_ids

    # fetch family specific ncbi_ids
    pdb_ids = fetch_value_list(rfam_acc, rs.PDB_IDs_QUERY)

    # fetch all author orcids associated with a family accession
    orcids = fetch_value_list(rfam_acc, rs.AU_ORCIDS)
    # pass orcid list in fam_fields

    # fetch family upids
    upids = fetch_value_list(rfam_acc, rs.FAMILY_UPIDS)

    # get pubmed ids
    pmids = get_value_list(fam_fields["pmids"], ",")

    # get dbxrefs (GO, SO)
    dbxrefs = get_value_list(fam_fields["dbxrefs"], ",")

    # get associated clan
    clan = fetch_value(rs.FAM_CLAN, rfam_acc)

    # get pseudoknot evidence
    pseudoknots = []

    # check if seed pseudoknot with covariation
    pseudoknot_evidence = int(fetch_value(rs.SEED_PK_WITH_COV, rfam_acc))
    if pseudoknot_evidence > 0:
        pseudoknots.append("seed with covariation support")

    # check if seed pseudoknot with no covariation
    pseudoknot_evidence = int(fetch_value(rs.SEED_PK_NO_COV, rfam_acc))
    if pseudoknot_evidence > 0:
        pseudoknots.append("seed no covariation support")

    # check if rscape pseudoknot with  covariation
    pseudoknot_evidence = int(fetch_value(rs.RSCAPE_PK_WITH_COV, rfam_acc))
    if pseudoknot_evidence > 0:
        pseudoknots.append("predicted with covariation support")

    # check if rscape pseudoknot with no covariation
    pseudoknot_evidence = int(fetch_value(rs.RSCAPE_PK_NO_COV, rfam_acc))
    if pseudoknot_evidence > 0:
        pseudoknots.append("predicted no covariation support")

    fam_fields["pseudoknots"] = pseudoknots

    if len(pseudoknots) > 0:
        fam_fields["has_pseudoknot"] = 1
    else:
        fam_fields["has_pseudoknot"] = 0

    # need a function here to split dbxrefs in a pretty way
    go_ids = [x for x in dbxrefs if x.find("GO") != -1]
    so_ids = [x for x in dbxrefs if x.find("SO") != -1]

    # update cross references dictionary
    cross_refs["ncbi_taxonomy_id"] = valid_ncbi_ids
    cross_refs["PDB"] = pdb_ids
    cross_refs["PUBMED"] = pmids
    cross_refs["GO"] = go_ids
    cross_refs["SO"] = so_ids
    cross_refs["ORCID"] = orcids

    if len(upids) > 0:
        cross_refs["UniProt"] = upids

    if clan is not None:
        cross_refs["RFAM"] = [clan]

    # add a new family entry to the xml tree
    entry = ET.SubElement(entries, "entry", id=rfam_acc)

    # entry name
    ET.SubElement(entry, "name").text = str(fam_fields["name"]).replace("_", " ")

    # entry description
    ET.SubElement(entry, "description").text = str(fam_fields["description"]).replace(
        "_", " "
    )

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
        entry,
        fam_fields,
        len(pdb_ids),
        valid_ncbi_ids,
        entry_type=entry_type,
        tax_strings=tax_strings,
    )

    if hfields is True:
        species_tax_trees = get_family_tax_tree(name_object, name_dict, valid_ncbi_ids)
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
        entry, clan_fields, clan_fields["num_families"], None, entry_type=entry_type
    )


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

    build_additional_fields(entry, motif_fields, 0, None, entry_type=entry_type)


# ----------------------------------------------------------------------------


def genome_xml_builder(entries, gen_acc=None):
    """
    Expands the Xml4dbDump with a Genome entry

    entries:    Entries node on xml tree
    gen_acc:  An Rfam associated motif accession
    """

    entry_type = "Genome"

    cross_ref_dict = {}

    # fetch genome fields
    genome_fields = fetch_entry_fields(gen_acc, rs.GENOME)

    # add a new genome entry to the xml tree
    entry = ET.SubElement(entries, "entry", id=gen_acc)

    if genome_fields["name"] is not None:
        ET.SubElement(entry, "name").text = genome_fields["name"]
    else:
        ET.SubElement(entry, "name").text = ""

    if genome_fields["description"] is not None:
        ET.SubElement(entry, "description").text = genome_fields["description"]
    else:
        ET.SubElement(entry, "description").text = ""

    # entry dates - common to motifs and clans
    dates = ET.SubElement(entry, "dates")

    created = genome_fields["created"].date().strftime("%d %b %Y")
    updated = genome_fields["updated"].date().strftime("%d %b %Y")

    ET.SubElement(dates, "date", value=created, type="created")
    ET.SubElement(dates, "date", value=updated, type="updated")

    # build genome cross references, if there are any
    if genome_fields["num_families"] != 0:
        genome_fams = fetch_value_list(gen_acc, rs.GENOME_FAMS)

        cross_ref_dict["RFAM"] = genome_fams
        cross_ref_dict["ncbi_taxonomy_id"] = [genome_fields["ncbi_id"]]
        build_cross_references(entry, cross_ref_dict)

    # genome additional fields
    build_genome_additional_fields(entry, genome_fields)


# ----------------------------------------------------------------------------


def result_iterator(cursor, arraysize=1000):
    """
    An iterator that uses fetchmany to keep memory usage down
    """
    while True:
        results = cursor.fetchmany(arraysize)
        if not results:
            break
        for result in results:
            yield result


# ----------------------------------------------------------------------------


def format_full_region(
    entries: ET.Element, region, genome, chromosome, rnacentral_ids: ty.Dict[str, str]
):
    """
    Format full regions for a genome. Genome metadata is retrieved only once.
    """
    timestamp = datetime.datetime.now().strftime("%d %b %Y")
    name = "%s/%s:%s" % (region["rfamseq_acc"], region["seq_start"], region["seq_end"])

    scientific_name = None
    if genome is not None:
        scientific_name = genome.scientific_name
    else:
        scientific_name = region["scientific_name"]

    description = "%s %s" % (scientific_name, region["rfam_id"])

    # add a new family entry to the xml tree
    entry_id = "%s_%s_%s" % (
        region["rfamseq_acc"],
        region["seq_start"],
        region["seq_end"],
    )
    entry = ET.SubElement(entries, "entry", id=entry_id)

    ET.SubElement(entry, "name").text = name
    ET.SubElement(entry, "description").text = description
    dates = ET.SubElement(entry, "dates")
    ET.SubElement(dates, "date", value=timestamp, type="created")
    ET.SubElement(dates, "date", value=timestamp, type="updated")

    # additional fields
    build_full_region_additional_fields(entry, region, genome, chromosome)

    # adding cross references
    cross_refs = {}

    ncbi_id = None
    if genome is not None:
        ncbi_id = genome.ncbi_id
    else:
        ncbi_id = region["ncbi_id"]

    # create cross references dictionary
    cross_refs["ncbi_taxonomy_id"] = [ncbi_id]
    cross_refs["RFAM"] = [region["rfam_acc"]]

    ena_accession = ""
    if region["rfamseq_acc"].find(".") != -1:
        ena_accession = region["rfamseq_acc"].partition(".")[0]
        cross_refs["ENA"] = [ena_accession]

    elif region["rfamseq_acc"][0:3] != "URS":
        ena_accession = region["rfamseq_acc"]
        cross_refs["ENA"] = [ena_accession]

    if genome is not None:
        cross_refs["Uniprot"] = [genome.upid]

    if name in rnacentral_ids:
        cross_refs["RNACENTRAL"] = [rnacentral_ids[name] + "_" + str(ncbi_id)]
    else:
        if "/" in name:
            name = name.partition("/")[0]
        cross_refs["RNACENTRAL"] = [name]

    build_cross_references(entry, cross_refs)


# ----------------------------------------------------------------------------


def get_chromosome_metadata():
    """
    Get chromosome metadata as a dictionary.
    """
    chromosomes = {}
    for genseq in Genseq.objects.exclude(chromosome_name__isnull=True).values():
        chromosomes[genseq["rfamseq_acc"]] = genseq
    return chromosomes


# ----------------------------------------------------------------------------


def get_rnacentral_mapping(upid):
    """
    Get RNAcentral mappings for all sequences from a given genome.
    """
    query = """
    SELECT rm.rfamseq_acc, rm.seq_start, rm.seq_end, rm.rnacentral_id
    FROM genseq gs, full_region fr, rnacentral_matches rm
    WHERE
    gs.rfamseq_acc = fr.rfamseq_acc
    AND fr.rfamseq_acc = rm.rfamseq_acc
    AND fr.seq_start = rm.seq_start
    AND fr.seq_end = rm.seq_end
    AND fr.is_significant = 1
    AND rm.rnacentral_id IS NOT NULL
    AND gs.upid = '%s'
    AND gs.version='15.0'
    """
    cnx = RfamDB.connect(db_config=DB_CONFIG)
    cursor = cnx.cursor(dictionary=True, buffered=True)
    cursor.execute(query % upid)
    rnacentral_ids = {}
    for row in result_iterator(cursor):
        rfamseq_acc = row["rfamseq_acc"]
        name = "%s/%s:%s" % (row["rfamseq_acc"], row["seq_start"], row["seq_end"])
        rnacentral_ids[name] = str(row["rnacentral_id"])
    cursor.close()
    cnx.disconnect()
    return rnacentral_ids


# ----------------------------------------------------------------------------


def full_region_xml_builder(entries: ET.Element, upid: str):
    """
    Export full region entries for a genome.

    entries:    Entries node on xml tree
    upid:  Genome identifier.
    """

    tax_id_duplicates = {
        "562": 1,
        "1280": 1,
        "7209": 1,
        "10679": 1,
        "10717": 1,
        "11021": 1,
        "11036": 1,
        "11072": 1,
        "11082": 1,
        "11228": 1,
        "11636": 1,
        "11963": 1,
        "31649": 1,
        "84589": 1,
        "90370": 1,
        "93838": 1,
        "186617": 1,
        "229533": 1,
        "351048": 1,
        "456327": 1,
        "766192": 1,
        "1891747": 1,
    }
    genome = None
    chromosomes = []
    if upid[0:2] == "UP":
        genome = Genome.objects.select_related("ncbi").get(upid=upid)
        chromosomes = get_chromosome_metadata()

    rnacentral_ids = get_rnacentral_mapping(upid=upid)
    cnx = RfamDB.connect(db_config=DB_CONFIG)
    cursor = cnx.cursor(dictionary=True, buffered=True)

    # work on 'full' refions
    cursor.execute(rs.FULL_REGION_FIELDS % upid)
    for row in result_iterator(cursor):
        format_full_region(entries, row, genome, chromosomes, rnacentral_ids)

    # work on 'seed' regions if not already exported
    # cursor.execute(rs.FULL_REGION_SEEDS % upid)

    """
    # if one of the cases of duplicates, work with the flags
    if genome.ncbi_id in tax_id_duplicates:
        if tax_id_duplicates[genome.ncbi_id] == 1:
            for row in result_iterator(cursor):
                format_full_region(entries, row, genome, chromosomes, rnacentral_ids)
            # set flag to 0 to disable export
            tax_id_duplicates[genome.ncbi_id] = 0
    """
    # capture the rest of the cases
    # else:
    cursor.execute(rs.FULL_REGION_SEEDS % upid)
    for row in result_iterator(cursor):
        format_full_region(entries, row, genome, chromosomes, rnacentral_ids)

    cursor.close()
    cnx.disconnect()


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
                ET.SubElement(cross_refs, "ref", dbkey=str(value), dbname=db_name)

    return cross_refs


# ----------------------------------------------------------------------------


def add_hierarchical_fields(xml_tree_node, tax_tree_dict, name_dict):
    """
    Expands the cross references xml tree by adding hierarchical references
    for the ncbi ids in valid_ncbi_ids

    xml_tree_node:  An existing xml tree node to expand with hierarchical
                    fields
    tax_tree_dict:  Species taxonomy tree dictionary as generated by
                    get_family_tax_tree
    name_dict:  NCBI's name dictionary as returned by read_ncbi_names_dmp
    """

    # add a new hierarchical ref for every tax_id in the family
    for tax_id in tax_tree_dict.keys():
        hfields = ET.SubElement(
            xml_tree_node, "hierarchical_field", name="taxonomy_lineage"
        )

        # fetch lineage
        lineage = tax_tree_dict[tax_id]
        tax_tree = lineage[::-1]

        for tax_tree_node in tax_tree:
            # create the root node
            if tax_tree_node == "1":
                # need to create one root node in order for the xml dump to
                # validate
                ET.SubElement(hfields, "root", label="root").text = "1"
            else:
                # skip root while creating child nodes
                if tax_tree_node != "1":
                    ET.SubElement(
                        hfields, "child", label=name_dict[tax_tree_node]
                    ).text = tax_tree_node


# ----------------------------------------------------------------------------


def build_additional_fields(
    entry, fields, num_3d_structures, fam_ncbi_ids, entry_type, tax_strings=None
):
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
    authors = authors.replace(";", ",")
    author_list = get_value_list(authors, rs.AUTH_DEL)

    for author in author_list:
        ET.SubElement(add_fields, "field", name="author").text = author

    if entry_type == "Family":

        # number of species
        ET.SubElement(add_fields, "field", name="num_species").text = str(
            fields["num_species"]
        )
        # number of 3D structures
        ET.SubElement(add_fields, "field", name="num_3d_structures").text = str(
            num_3d_structures
        )
        # num seed
        ET.SubElement(add_fields, "field", name="num_seed").text = str(
            fields["num_seed"]
        )
        # num full
        ET.SubElement(add_fields, "field", name="num_full").text = str(
            fields["num_full"]
        )

        # rna types
        rna_types = get_value_list(fields["rna_type"], rs.RNA_TYPE_DEL)

        for rna_type in rna_types:
            ET.SubElement(add_fields, "field", name="rna_type").text = rna_type

        # has 3d structure
        if num_3d_structures > 0:
            ET.SubElement(add_fields, "field", name="has_3d_structure").text = "Yes"
        else:
            ET.SubElement(add_fields, "field", name="has_3d_structure").text = "No"

        # add popular species if any
        for species in rs.POPULAR_SPECIES:
            if species in fam_ncbi_ids:
                ET.SubElement(add_fields, "field", name="popular_species").text = str(
                    species
                )

        for tax_string in tax_strings:
            ET.SubElement(add_fields, "field", name="tax_string").text = str(tax_string)

        # work on pseudoknots here
        if fields["has_pseudoknot"] == 1:
            pseudoknots = fields["pseudoknots"]

            ET.SubElement(add_fields, "field", name="has_pseudoknot").text = "Yes"

            for pk_evidence in pseudoknots:
                ET.SubElement(
                    add_fields, "field", name="pseudoknot_evidence"
                ).text = pk_evidence
        else:
            ET.SubElement(add_fields, "field", name="has_pseudoknot").text = "No"

            # build hierarchical_fields tree here...

    # perhaps move this to clan and motif xml builder
    else:
        num_families = None

        if entry_type == "Motif":
            num_families = fetch_value(rs.NUM_FAMS_MOTIF, fields["id"])

        elif entry_type == "Clan":
            num_families = fields["num_families"]

        ET.SubElement(add_fields, "field", name="num_families").text = str(num_families)

    # returning node
    return add_fields


# ----------------------------------------------------------------------------


def build_genome_additional_fields(entry, fields):
    """
    Builds additional field nodes for a Genome

    entry:  This is the xml.etree.ElementTree at the point of entry
    fields: A list of additional fields to expand the entry with

    return: void
    """

    # TO DO - Generalize this one by executing a query to fetch additional
    # fields here

    add_fields = ET.SubElement(entry, "additional_fields")

    # adding entry type
    ET.SubElement(add_fields, "field", name="entry_type").text = "Genome"

    if fields["assembly_acc"] is None:
        ET.SubElement(add_fields, "field", name="gca_accession").text = ""
    else:
        ET.SubElement(add_fields, "field", name="gca_accession").text = fields[
            "assembly_acc"
        ]

    ET.SubElement(add_fields, "field", name="length").text = str(fields["total_length"])
    ET.SubElement(add_fields, "field", name="tax_string").text = fields["tax_string"]
    ET.SubElement(add_fields, "field", name="ncbi_taxid").text = str(fields["ncbi_id"])
    ET.SubElement(add_fields, "field", name="num_rfam_hits").text = str(
        fields["num_rfam_regions"]
    )
    ET.SubElement(add_fields, "field", name="num_families").text = str(
        fields["num_families"]
    )
    ET.SubElement(add_fields, "field", name="scientific_name").text = str(
        fields["name"]
    )  # redundant

    # add popular species if any
    species = str(fields["ncbi_id"])
    if species in rs.POPULAR_SPECIES:
        ET.SubElement(add_fields, "field", name="popular_species").text = species

    if fields["common_name"] is not None:
        ET.SubElement(add_fields, "field", name="common_name").text = str(
            fields["common_name"]
        )

    if fields["assembly_level"] is not None:
        ET.SubElement(add_fields, "field", name="assembly_level").text = str(
            fields["assembly_level"]
        )

    if fields["assembly_name"] is not None:
        ET.SubElement(add_fields, "field", name="assembly_name").text = str(
            fields["assembly_name"]
        )

    return add_fields


# ----------------------------------------------------------------------------


def build_full_region_additional_fields(entry, fields, genome, chromosomes):
    """
    Builds additional field nodes for a the full_region xml dump

    entry:  This is the xml.etree.ElementTree at the point of entry
    fields: A list of additional fields to expand the entry with

    return: void
    """

    # TO DO - Generalize this one by executing a query to fetch additional
    # fields here

    add_fields = ET.SubElement(entry, "additional_fields")

    tax_string = ""
    species = ""
    common_name = ""
    scientific_name = ""
    if genome is not None:
        tax_string = genome.ncbi.tax_string
        species = genome.ncbi_id
        common_name = genome.common_name
        scientific_name = genome.scientific_name
    else:
        tax_string = fields["tax_string"]
        species = fields["ncbi_id"]
        scientific_name = fields["scientific_name"]

    # adding entry type
    ET.SubElement(add_fields, "field", name="entry_type").text = "Sequence"
    ET.SubElement(add_fields, "field", name="rfamseq_acc").text = str(
        fields["rfamseq_acc"]
    )
    ET.SubElement(add_fields, "field", name="rfamseq_acc_description").text = str(
        fields["rfamseq_acc_description"]
    )
    ET.SubElement(add_fields, "field", name="seq_start").text = str(fields["seq_start"])
    ET.SubElement(add_fields, "field", name="seq_end").text = str(fields["seq_end"])
    ET.SubElement(add_fields, "field", name="cm_start").text = str(fields["cm_start"])
    ET.SubElement(add_fields, "field", name="cm_end").text = str(fields["cm_end"])

    ET.SubElement(add_fields, "field", name="evalue_score").text = str(
        fields["evalue_score"]
    )
    ET.SubElement(add_fields, "field", name="bit_score").text = str(fields["bit_score"])
    ET.SubElement(add_fields, "field", name="alignment_type").text = str(
        fields["alignment_type"]
    )
    ET.SubElement(add_fields, "field", name="truncated").text = str(fields["truncated"])
    # ET.SubElement(add_fields, "field", name="tax_string").text = genome.ncbi.tax_string
    ET.SubElement(add_fields, "field", name="tax_string").text = tax_string

    if fields["rfamseq_acc"] in chromosomes:
        ET.SubElement(add_fields, "field", name="chromosome_name").text = chromosomes[
            fields["rfamseq_acc"]
        ]["chromosome_name"]
        ET.SubElement(add_fields, "field", name="chromosome_type").text = chromosomes[
            fields["rfamseq_acc"]
        ]["chromosome_type"]

    # add popular species if any
    # species = genome.ncbi_id
    if species in rs.POPULAR_SPECIES:
        ET.SubElement(add_fields, "field", name="popular_species").text = species

    # rna types
    rna_types = get_value_list(fields["rna_type"], rs.RNA_TYPE_DEL)

    for rna_type in rna_types:
        ET.SubElement(add_fields, "field", name="rna_type").text = rna_type

    if common_name:
        # ET.SubElement(add_fields, "field", name="common_name").text = genome.common_name
        ET.SubElement(add_fields, "field", name="common_name").text = common_name
    # ET.SubElement(add_fields, "field", name="scientific_name").text = genome.scientific_name
    ET.SubElement(add_fields, "field", name="scientific_name").text = scientific_name

    return add_fields


# ----------------------------------------------------------------------------


def get_value_list(val_str, delimiter=","):
    """
    Splits an input string using delimiter and returns a list of the
    elements

    val_str:    A string of family specific values. This string is a
                concatenation of multiple values related to a single family
    delimiter:  The delimeter that will be used to split the values' string
    """

    val_str = val_str.strip()
    # split string
    values = val_str.split(delimiter)
    # filter values
    value_list = [x.strip() for x in values if x != ""]

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

    cnx = RfamDB.connect(db_config=DB_CONFIG)

    cursor = cnx.cursor(raw=True)

    if rfam_acc is None:
        cursor.execute(query)
    else:
        cursor.execute(query % rfam_acc)

    values = cursor.fetchall()

    cursor.close()
    cnx.disconnect()

    if len(values) > 0:
        if isinstance(values[0], tuple):
            return [str(x[0]) for x in values]
        else:
            return [str(x) for x in values]

    return []


# ----------------------------------------------------------------------------


def fetch_entry_fields(entry_acc, entry_type):
    """
    Returns a dictionary with the entry's fields

    entry_acc:  An Rfam associated accession (Motif, Clan, Family)
    entry_type: The type of the entry accession
    """

    # maybe the entry type not required... use rfam_acc[0:2]

    cnx = RfamDB.connect(db_config=DB_CONFIG)
    cursor = cnx.cursor(dictionary=True)

    entry_type = entry_type[0].capitalize()

    fields = None

    try:
        if entry_type == rs.FAMILY:
            cursor.execute(rs.FAM_FIELDS % entry_acc)

        elif entry_type == rs.CLAN:
            cursor.execute(rs.CLAN_FIELDS % entry_acc)

        elif entry_type == rs.MOTIF:
            cursor.execute(rs.MOTIF_FIELDS % entry_acc)

        elif entry_type == rs.GENOME:
            cursor.execute(rs.GENOME_FIELDS % entry_acc)

        fields = cursor.fetchall()[0]

    except:
        print("Failure retrieving values for entry %s." % entry_acc)

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

    cnx = RfamDB.connect(db_config=DB_CONFIG)

    cursor = cnx.cursor(raw=True)

    if accession is not None:
        cursor.execute(query % accession)

    else:
        cursor.execute(query)

    value = cursor.fetchall()

    cursor.close()
    cnx.disconnect()

    if len(value) > 0:
        return value[0][0]

    return None


# ----------------------------------------------------------------------------


def main(entry_type, rfam_acc, outdir, hfields=False):
    """
    This function puts everything together

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
                print("Error creating output directory at: ", outdir)

        # export all entries
        if rfam_acc is None:
            # Motif accessions
            if entry_type == rs.MOTIF:
                rfam_accs = fetch_value_list(None, rs.MOTIF_ACC)

            # Clan accessions
            elif entry_type == rs.CLAN:
                rfam_accs = fetch_value_list(None, rs.CLAN_ACC)

            # Genome accessions
            elif entry_type == rs.GENOME:
                rfam_accs = fetch_value_list(None, rs.GENOME_ACC)

            # Genome accessions required for exporting full region
            elif entry_type == rs.MATCH:
                rfam_accs = fetch_value_list(None, rs.GENOME_ACC)

            # Family accessions
            elif entry_type == rs.FAMILY:
                if hfields:
                    # load ncbi taxonomy browser here
                    name_dict, name_dict_reverse = read_ncbi_names_dmp(
                        rfc.TAX_NAMES_DUMP
                    )
                    name_object = read_ncbi_taxonomy_nodes(
                        name_dict, rfc.TAX_NODES_DUMP
                    )

                rfam_accs = fetch_value_list(None, rs.FAM_ACC)

                for entry in rfam_accs:
                    print(entry)
                    t0 = timeit.default_timer()
                    xml4db_dumper(
                        name_dict, name_object, entry_type, entry, hfields, outdir
                    )
                    print(
                        "%s execution time: %.1fs"
                        % (entry, timeit.default_timer() - t0)
                    )

                return

            # Don't build hierarchical references for Clans and Motifs
            for entry in rfam_accs:
                print(entry)
                entry_file = os.path.join(outdir, entry + ".xml")
                if os.path.exists(entry_file):
                    print(
                        "Skipping {acc}. File already exists {f}".format(
                            acc=entry, f=entry_file
                        )
                    )
                else:
                    t0 = timeit.default_timer()
                    xml4db_dumper(None, None, entry_type, entry, False, outdir)
                    print(
                        "%s execution time: %.1fs"
                        % (entry, timeit.default_timer() - t0)
                    )

        # export single entry
        else:
            # need to check the validity of an rfam_acc (rfam, motif, clan)
            if (
                entry_type == rs.MOTIF
                or entry_type == rs.CLAN
                or entry_type == rs.GENOME
            ):
                xml4db_dumper(None, None, entry_type, rfam_acc, False, outdir)

            # export single family entry
            else:
                if hfields:
                    # load ncbi taxonomy browser here
                    name_dict, name_dict_reverse = read_ncbi_names_dmp(
                        rfc.TAX_NAMES_DUMP
                    )
                    name_object = read_ncbi_taxonomy_nodes(
                        name_dict, rfc.TAX_NODES_DUMP
                    )

                xml4db_dumper(
                    name_dict, name_object, entry_type, rfam_acc, hfields, outdir
                )

    except:
        traceback.print_exc()
        # need to correct this one
        if rfam_acc is None:
            gen_fams = Set([x.partition(".")[0] for x in os.listdir(outdir)])
            loaded_fams = Set(rfam_accs)
            # get remaining families
            rem_fams = loaded_fams - gen_fams

            # open a log file
            logging.basicConfig(
                filename=os.path.join("missing_accs" + ".log"),
                filemode="w",
                level=logging.DEBUG,
            )

            # write accessions to log file
            for rfam_acc in rem_fams:
                logging.debug(rfam_acc)
        else:
            print("Error exporting %s." % rfam_acc)


# ----------------------------------------------------------------------------


def get_valid_family_tax_ids(name_object, family_tax_ids):
    """
    Returns a list of all family tax ids found in the NCBI dumps

    name_object: NCBI tax browser node dictionary
    family_tax_ids: A list of all family ncbi ids
    """

    valid_family_tax_ids = []

    for taxid in family_tax_ids:
        if taxid in name_object:
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

        if taxid in name_object:
            species_tax_trees[taxid] = name_object[taxid].get_lineage(name_object)

    return species_tax_trees


# ----------------------------------------------------------------------------


def xmllint(filepath):
    """
    Validate xml files against EBI Search schema.
    Run xmllint on the output file and print the resulting report.
    """
    schema_url = "http://www.ebi.ac.uk/ebisearch/XML4dbDumps.xsd"
    cmd = ("xmllint {filepath} --schema {schema_url} --noout --stream").format(
        filepath=filepath, schema_url=schema_url
    )
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print("ERROR: xmllint validation failed")
        print(e.output)
        print(e.cmd)
        print("Return code {0}".format(e.returncode))
        sys.exit(1)


# ----------------------------------------------------------------------------


def usage():
    """
    Parses arguments and displays usage information on screen
    """

    parser = argparse.ArgumentParser(
        description="Rfam Search Xml4db Dumper.", epilog=""
    )

    # group required arguments together
    req_args = parser.add_argument_group("required arguments")

    req_args.add_argument(
        "--type",
        help="rfam entry type (F: Family, M: Motif, C: Clan, G: Genome, R: Regions)",
        type=str,
        choices=["F", "M", "C", "G", "R"],
        required=True,
    )

    parser.add_argument(
        "--acc",
        help="a valid rfam entry accession (RF*|CL*|RM*)",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--hfields", help="include hierarchical fields", action="store_true"
    )

    req_args.add_argument(
        "--out", help="path to output directory", type=str, required=True
    )

    parser.add_argument(
        "--db", help="database to use - rel or live", type=str, default=None
    )

    return parser


# ----------------------------------------------------------------------------

if __name__ == "__main__":

    parser = usage()
    args = parser.parse_args()

    # Additional checks
    # Check if export type matches accession
    wrong_input = False
    if args.acc is not None:
        if args.type == "F" and args.acc[0:2] != "RF":
            wrong_input = True
        elif args.type == "M" and args.acc[0:2] != "RM":
            wrong_input = True
        elif args.type == "C" and args.acc[0:2] != "CL":
            wrong_input = True
        elif args.type == "G" and (args.acc[0:2] != "UP" and args.acc[0:2] != "RG"):
            wrong_input = True

    if wrong_input is True:
        print("\nAccession does not match the export type.\n")
        parser.print_help()
        sys.exit()

    # check output directory
    if os.path.isdir(args.out) is False:
        print("\nPlease provide a valid output directory.\n")
        parser.print_help()
        sys.exit()

    if not args.db:
        DB_CONFIG = RFAMLIVE
    else:
        if args.db == "rfamrel":
            DB_CONFIG = RFAMREL
        elif args.db == "rfamlive":
            DB_CONFIG = RFAMLIVE
        else:
            print("\nPlease provide a valid database option: rfamrel or rfamlive\n")

    main(args.type, args.acc, args.out, hfields=args.hfields)
