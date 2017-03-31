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

#TO DO:   - Rename functions
#         - Split fetch_gca_data

# ---------------------------------IMPORTS-------------------------------------

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

    assembly_links = None
    assembly_level = assembly.find("ASSEMBLY_LEVEL").text
    assembly_links = assembly.find("ASSEMBLY_LINKS")

    if assembly_level == "complete genome":
        assembly_level = assembly_level.replace(' ', '-')

    fields["assembly_level"] = assembly_level

    if assembly_level == "contig" and assembly_links is None:
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
    dict format to be easily loaded via Django ORM

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

        if accession[0].find('.') != -1:
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

                fields["mol_type"] = None
                if acc_meta.has_key("mol_type"):
                    fields["mol_type"] = acc_meta["mol_type"]

                # primary, PATCHES, ALT_REF_LOCI_1
                fields["assembly_unit"] = accession[6]
                # this takes the date of the entry is created
                entry_date = datetime.datetime.now().strftime(
                    "%Y-%m-%d %H:%M:%S")
                fields["created"] = entry_date
                fields["updated"] = entry_date

            elif len(accession) < 7 and accession[0].find('.') != -1:
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

    if entry is None:
        return fields

    fields["seq_version"] = 0
    if entry.get("entryVersion") is not None:
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
    type, description and tax id

    accession: A valid GCA accession from ENA
    """

    metadata = {}
    xml_str = requests.get(gc.ENA_XML_URL % accession).content

    xml_root = ET.fromstring(xml_str)

    entry = xml_root.find("entry")

    # None if no moleculeType found
    mol_type = entry.find("moleculeType")

    if mol_type is not None:
        metadata["mol_type"] = entry.get("moleculeType")

    # get molecule description and tax id
    metadata["description"] = entry.find("description").text
    metadata["ncbi_id"] = int(entry.find("feature").find("taxon").get("taxId"))

    return metadata


# -----------------------------------------------------------------------------


def fetch_assembly_links(gca_acc):
    """
    Retrieves and returns a dictionary with all ftp links found in the GCA xml
    file

    gca_acc: A valid GCA accession from ENA
    """

    gca_ftp_links = {}

    response = requests.get(gc.ENA_XML_URL % gca_acc).content

    xml = ET.fromstring(response)

    # fetch assembly links nodet
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
def get_general_accession_metadata(upid, accession_list):
    """
    Extracts the metadata for a list of accessions which correspond to a
    specific genome. This function is for the case in which accessions are
    extracted from Uniprot's proteome rdf files

    upid: Uniprot's proteome id
    accession_list: A list of valid sequence accessions to be used for
    extracting the metadata

    return: A json like structure to be used for generating the json dumps
    """

    genome_entries = []
    entry = {}
    fields = {}

    for acc in accession_list:
        entry["model"] = gc.GENSEQ_MODEL
        entry["pk"] = acc
        fields = fetch_wgs_acc_metadata(acc)
        entry["fields"] = fields
        # adding upid and wgs_acc in entry fields
        entry["fields"]["upid"] = upid
        entry["fields"]["gen_acc"] = None
        genome_entries.append(entry)
        entry = {}

    return genome_entries

# -----------------------------------------------------------------------------

def extract_uniprot_genome_metadata(upid):
    """
    Parses a proteome's xml file from Uniprot and converts it to a json
    like object which is returned

    upid: A valid Uniprot Proteome identifier

    returns: A dictionary
    """

    proteome_dict = {}

    # namespace prefix # or register a namespace in the ET
    prefix = "{http://uniprot.org/uniprot}%s"

    response = requests.get(gc.PROTEOME_XML_URL % upid).content
    prot_tree_root = ET.fromstring(response)

    proteome = prot_tree_root.find(prefix % "proteome")

    upid = proteome.find(prefix % "upid")
    proteome_dict["upid"] = upid.text
    proteome_dict["ncbi_id"] = int(proteome.find(prefix % "taxonomy").text)
    proteome_dict["description"] = proteome.find(prefix % "description").text
    proteome_dict["scientific_name"] = proteome.find(prefix % "name").text

    # get sequence accessions
    acc_dict = {}
    other_accs = []
    accession_nodes = proteome.findall(prefix % "component")

    for node in accession_nodes:

        if node.get("name").find("WGS") != -1:
            acc_dict["WGS"] = node.find(prefix % "genome_accession").text

        elif node.get("name").find("Chloroplast") != -1:
            acc_dict["Chloroplast"] = node.find(prefix % "genome_accession").text

        elif node.get("name").find("Mitochondrion") != -1:
            acc_dict["Mitochondrion"] = node.find(prefix % "genome_accession").text

        else:
            other_accs.append(node.find(prefix % "genome_accession").text)

    if len(other_accs) > 0:
        acc_dict["other"] = other_accs

    proteome_dict["accessions"] = acc_dict

    return proteome_dict

# -----------------------------------------------------------------------------

def dump_uniprot_genome_metadata(upid, kingdom):
    """
    Parses ENA GCA accession xml, and returns the accession's data in the
    form of a dictionary

    proteome_dict: A proteome dict built from Uniprot's proteome xml files
    upid: A valid Uniprot proteome accession (e.g. UP000005640 - Homo Sapiens)
    assembly_acc: A valid ENA GCA accession
    kingdom: The corresponding species kingdom
    """

    genome_entry = {}
    fields = {}

    # extract the values from the corresponding xml
    proteome_dict = extract_uniprot_genome_metadata(upid)

    fields["upid"] = upid
    fields["gca_acc"] = None
    fields["gca_version"] = None

    # add ensembl fields as a post-processing step
    fields["ensembl_id"] = None
    fields["ensembl_source"] = None

    fields["assembly_name"] = None
    fields["assembly_level"] = None

    # check if there is a WGS accession
    if "WGS" in proteome_dict["accessions"].keys():
        fields["wgs_acc"] = proteome_dict["accessions"]["WGS"]
        fields["wgs_version"] = int(proteome_dict["accessions"]["WGS"][4:6])
    else:
        fields["wgs_acc"] = None
        fields["wgs_version"] = None

    fields["study_ref"] = None

    # description can be very long, using title instead as short description
    fields["description"] = proteome_dict["description"]

    fields["total_length"] = None
    fields["ungapped_length"] = None
    fields["circular"] = None

    fields["ncbi_id"] = proteome_dict["ncbi_id"]

    fields["scientific_name"] = proteome_dict["scientific_name"]
    fields["common_name"] = None
    fields["kingdom"] = kingdom

    fields["num_gen_regions"] = None
    fields["num_rfam_regions"] = None
    fields["num_families"] = None

    # this takes the date of the entry is created
    entry_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fields["created"] = entry_date
    fields["updated"] = entry_date

    genome_entry["model"] = gc.GENOME_MODEL
    genome_entry["pk"] = upid
    genome_entry["fields"] = fields

    # return genome dictionary (in Django format)
    return genome_entry

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    pass
