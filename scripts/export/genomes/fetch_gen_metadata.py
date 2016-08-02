"""
Created on 12 Jul 2016

@author: ikalvari

TO DO:   - Rename functions
         - Split fetch_gca_data
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import datetime
import xml.etree.ElementTree as ET
import requests

# Config Files
import genome_fetch as gf
from config import gen_config as gc

# -----------------------------------------------------------------------------


def fetch_gca_data(upid, assembly_acc, kingdom):
    """
    Parses ENA GCA accession xml, and returns the accession's data in the
    form of a dictionary

    upid: A valid Uniprot proteome accession (e.g. UP000005640 - Homo Sapiens)
    assembly_acc: A valid ENA GCA accession
    kingdom: The corresponding species kingdom
    """

    genome_entry = {}
    fields = {}
    tmp_acc = assembly_acc

    '''
    if tmp_acc.find('.') != -1:
        tmp_acc = tmp_acc.partition('.')[0]
    '''

    assembly_xml = requests.get(gc.ENA_XML_URL % tmp_acc).content

    root = ET.fromstring(assembly_xml)
    assembly = root.find("ASSEMBLY")

    primary_id = assembly.find("IDENTIFIERS").find("PRIMARY_ID")

    if primary_id is not None:
        fields["gca_acc"] = primary_id.text
    else:
        fields["gca_acc"] = primary_id

    version = fields["gca_acc"].partition('.')[2]
    fields["gca_version"] = int(version)

    # add ensembl fields as a post-processing step
    fields["ensembl_id"] = None
    fields["ensembl_source"] = None
    fields["assembly_name"] = assembly.find("NAME").text

    ass_links = None
    ass_level = assembly.find("ASSEMBLY_LEVEL").text
    ass_links = assembly.find("ASSEMBLY_LINKS")

    fields["assembly_level"] = ass_level

    if ass_level == "contig" and ass_links is None:
        wgs_fields = assembly.find("WGS_SET")
        wgs_acc = gf.get_wgs_set_accession(
            wgs_fields.find("PREFIX").text, wgs_fields.find("VERSION").text)

        fields["wgs_acc"] = wgs_acc
        fields["wgs_version"] = int(wgs_fields.find("VERSION").text)
    else:
        fields["wgs_acc"] = None
        fields["wgs_version"] = None

    fields["study_ref"] = assembly.find("STUDY_REF").find(
        "IDENTIFIERS").find("PRIMARY_ID").text

    # description can be very long, using title instead as short description
    fields["description"] = assembly.find("TITLE").text

    attributes = fetch_assembly_attributes(
        assembly.find("ASSEMBLY_ATTRIBUTES"))

    fields["total_length"] = int(attributes["total-length"])
    fields["ungapped_length"] = int(attributes["ungapped-length"])

    genome_desc = assembly.find("DESCRIPTION").text

    if genome_desc is not None:
        if genome_desc.find("circular") != -1:
            fields["circular"] = 1
        else:
            fields["circular"] = 0
    else:
        fields["circular"] = None

    taxid = assembly.find("TAXON")
    fields["ncbi_id"] = int(taxid.find("TAXON_ID").text)

    fields["scientific_name"] = taxid.find("SCIENTIFIC_NAME").text

    common_name = None
    common_name = taxid.find("COMMON_NAME")
    if common_name is not None:
        fields["common_name"] = taxid.find("COMMON_NAME").text
    fields["kingdom"] = kingdom

    fields["num_gen_regions"] = int(attributes["count-regions"])
    fields["num_rfam_regions"] = None
    fields["num_families"] = None

    # this takes the date of the entry is created
    entry_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fields["created"] = entry_date
    fields["updated"] = entry_date

    genome_entry["model"] = gc.GENOME_MODEL
    genome_entry["pk"] = upid
    genome_entry["fields"] = fields

    return genome_entry

# -----------------------------------------------------------------------------


def fetch_assembly_accessions(upid, gca_acc, acc_ftp_link, reg_ftp_link=None):
    """
    Parses assembly report file and exports all assembly accessions in a
    django dict format..

    upid: A valid Uniprot proteome accession (e.g. UP000005640 - Homo Sapiens)
    gca_acc: A valid ENA GCA accession
    acc_ftp_link: Assembly report file ftp url (as retrieved from GCA xml file)
    reg_ftp_link: Assembly region file ftp url (as retrieved from GCA xml file)
    """

    assembly_accs = []
    fields = {}
    entry = {}
    regions = None

    # load accession regions
    if reg_ftp_link is not None:
        regions = region_loader(reg_ftp_link)

    http_link = acc_ftp_link.replace("ftp://", "http://")
    response = requests.get(http_link).content

    accessions = response.strip().split('\n')
    # remove header
    accessions.pop(0)

    accession = ''
    for acc_line in accessions:
        accession = acc_line.strip().split('\t')

        if (accession[0].find('.') != -1):
            entry["model"] = gc.GENSEQ_MODEL
            entry["pk"] = str(accession[0])  # seq_acc

            if len(accession) == 7:
                acc_meta = fetch_gca_acc_metadata(accession[0])

                fields["ncbi_id"] = acc_meta["ncbi_id"]
                fields["description"] = acc_meta["description"]
                fields["upid"] = upid
                fields["gen_acc"] = gca_acc
                fields["seq_version"] = int(accession[0].partition('.')[2])
                fields["seq_length"] = int(accession[2])
                fields["seq_role"] = accession[3]  # patch, loci etc

                # need to parse the regions file for these fields
                if reg_ftp_link is not None and (accession[0] in regions.keys()):
                    fields["seq_start"] = int(regions[accession[0]][0])
                    fields["seq_end"] = int(regions[accession[0]][1])
                else:
                    fields["seq_start"] = 0
                    fields["seq_end"] = 0

                fields["mol_type"] = acc_meta["mol_type"]
                # primary, PATCHES, ALT_REF_LOCI_1
                fields["assembly_unit"] = accession[6]
                # this takes the date of the entry is created
                entry_date = datetime.datetime.now().strftime(
                    "%Y-%m-%d %H:%M:%S")
                fields["created"] = entry_date
                fields["updated"] = entry_date

            elif(len(accession) < 7 and accession[0].find('.') != -1):
                fields = fetch_wgs_acc_metadata(accession[0])

        else:
            # skip cases no accession is provided
            pass

        entry["fields"] = fields
        assembly_accs.append(entry)
        fields = {}
        entry = {}

    return assembly_accs

# -----------------------------------------------------------------------------


def region_loader(reg_ftp_link):
    """
    Parses an assembly's region file and builds a dictionary of accessions and
    start-end coordinates which are stored in a tuple format

    reg_ftp_link: The ftp url of the corresponding region file
    """

    region_dict = {}

    http_link = reg_ftp_link.replace("ftp://", "http://")
    response = requests.get(http_link).content

    regions = response.strip().split('\n')

    # remove header
    regions.pop(0)

    for region in regions:
        region = region.strip().split('\t')
        coords = region[5].split('-')
        region_dict[region[1]] = (int(coords[0]), int(coords[1]))

    return region_dict

# -----------------------------------------------------------------------------


def fetch_wgs_metadata(upid, assembly_acc, kingdom):
    """
    Parses ENA WGS accession xml, and returns the accession's data in the
    form of a dictionary

    upid: A valid Uniprot's proteome id
    assembly_acc: A valid ENA's WGS accession
    kingdom: A string representing the kingdom a species belongs to (e.g.viruses)
    """

    wgs_entry = {}
    fields = {}

    response = requests.get(gc.ENA_XML_URL % assembly_acc).content
    assembly_xml = ET.fromstring(response)  # root
    entry = assembly_xml.find("entry")

    fields["gca_acc"] = None
    fields["gca_version"] = None

    fields["wgs_acc"] = entry.get("accession")
    fields["wgs_version"] = int(entry.get("version"))
    fields["ensembl_id"] = None  # post processing

    # fetching assembly name out of entry's comment lines

    entry_cmnt = None
    entry_cmnt = entry.find("comment")

    if entry_cmnt is not None:
        assembly_name = get_wgs_assembly_name(
            entry.find("comment").text)
        fields["assembly_name"] = assembly_name

    else:
        fields["assembly_name"] = None

    if entry.get("topology") == "linear":
        fields["circular"] = 0
    else:
        fields["circular"] = 1

    fields["assembly_level"] = "contig"

    fields["study_ref"] = entry.find("projectAccession").text
    fields["description"] = entry.find("reference").find("title").text
    fields["kingdom"] = kingdom

    taxon = entry.find("feature").find("taxon")
    fields["scientific_name"] = taxon.get("scientificName")
    fields["common_name"] = None
    fields["ncbi_id"] = int(taxon.get("taxId"))

    fields["total_length"] = None
    fields["ungapped_length"] = None

    fields["num_gen_regions"] = None

    fields["num_rfam_regions"] = None  # post process
    fields["num_families"] = None  # post process

    # this takes the date of the entry is created
    entry_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fields["created"] = entry_date
    fields["updated"] = entry_date

    wgs_entry["model"] = gc.GENOME_MODEL
    wgs_entry["pk"] = upid
    wgs_entry["fields"] = fields

    return wgs_entry

# -----------------------------------------------------------------------------


def fetch_wgs_accs_metadata(upid, assembly_acc, wgs_range):
    """
    This function runs over the wgs range and exports all metadata for every
    accession in the provided wgs range

    wgs_range: A valid ENA-WGS set range
    """

    wgs_entries = []
    entry = {}
    fields = {}
    wgs_accs = gf.fetch_wgs_range_accs(wgs_range)

    for acc in wgs_accs:
        entry["model"] = gc.GENSEQ_MODEL
        entry["pk"] = acc
        fields = fetch_wgs_acc_metadata(acc)
        entry["fields"] = fields
        # adding upid and wgs_acc in entry fields
        entry["fields"]["upid"] = upid
        entry["fields"]["gen_acc"] = assembly_acc
        wgs_entries.append(entry)
        entry = {}

    return wgs_entries

# -----------------------------------------------------------------------------


def fetch_wgs_acc_metadata(wgs_acc):
    """
    Return a fields dictionary

    wgs_acc: A valid ENA wgs accession
    """

    fields = {}
    response = requests.get(gf.ENA_XML_URL % wgs_acc).content
    acc_xml = ET.fromstring(response)
    entry = acc_xml.find("entry")

    fields["seq_version"] = int(entry.get("entryVersion"))
    fields["seq_length"] = int(entry.get("sequenceLength"))
    fields["seq_role"] = "contig"
    fields["mol_type"] = entry.get("moleculeType")

    taxon_node = None
    taxon_node = entry.find("feature").find("taxon")
    fields["ncbi_id"] = int(taxon_node.get("taxId"))

    location = None
    refs = None
    refs = entry.findall("reference")

    for ref in refs:
        location = ref.get("location")
        if location is not None:
            break

    if location is not None:
        location = location.strip().split('-')
        fields["seq_start"] = int(location[0])
        fields["seq_end"] = int(location[1])
    else:
        fields["seq_start"] = 0
        fields["seq_end"] = 0

    fields["assembly_unit"] = None  # need to check this one
    fields["description"] = entry.find("description").text

    # this takes the date of the entry is created
    entry_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fields["created"] = entry_date
    fields["updated"] = entry_date  # when adding a new genome

    return fields
# -----------------------------------------------------------------------------


def get_wgs_assembly_name(wgs_comment):
    """
    Search the wgs comment text and return the Assembly name

    wgs_comment: WGS comment text
    """

    assembly_name = ''
    cmnt_lines = wgs_comment.strip().split('\n')

    for line in cmnt_lines:
        if line.find("Assembly Name") != -1:
            assembly_name = line
            break

    assembly_name = assembly_name.strip().partition("::")[2].strip()

    return assembly_name
# -----------------------------------------------------------------------------


def fetch_assembly_attributes(attrs_node):
    """
    Runs over the attributes node in the xml file and returns all attribute
    tag-value pairs in a dictionary format

    attrs_node: GCA xml's attributes node
    """

    attribute_values = {}
    attributes = attrs_node.findall("ASSEMBLY_ATTRIBUTE")

    for attr in attributes:
        tag = attr.find("TAG").text
        value = attr.find("VALUE").text
        attribute_values[tag] = value

    return attribute_values

# -----------------------------------------------------------------------------


def fetch_gca_acc_metadata(accession):
    """
    Fetch accession metadata and return a dictionary with ncbi_id, molecule's
    type, description

    accession: A valid GCA related ENA accession
    """

    metadata = {}
    xml_str = requests.get(gc.ENA_XML_URL % accession).content

    xml_root = ET.fromstring(xml_str)

    entry = xml_root.find("entry")

    '''
    mol_type = None
    mol_type = entry.get("moleculeType")
    '''

    # None if no moleculeType found
    metadata["mol_type"] = entry.get("moleculeType")
    metadata["description"] = entry.find("description").text
    metadata["ncbi_id"] = int(entry.find("feature").find("taxon").get("taxId"))

    return metadata

# -----------------------------------------------------------------------------


def fetch_assembly_links(gca_acc):
    """
    Retrieves and returns a dictionary with all ftp links found in the GCA xml
    file.

    gca_acc: A valid ENA GCA accession
    """

    gca_ftp_links = {}

    response = requests.get(gc.ENA_XML_URL % gca_acc).content

    xml = ET.fromstring(response)

    # fetch assembly links node
    assembly_links = xml.find("ASSEMBLY").find(
        "ASSEMBLY_LINKS").findall("ASSEMBLY_LINK")

    # loop over all available links
    for link in assembly_links:
        label = link.find("URL_LINK").find("LABEL").text.replace(' ', '_')
        url = link.find("URL_LINK").find("URL").text

        gca_ftp_links[label] = url

        label = None
        url = None

    return gca_ftp_links

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    pass
