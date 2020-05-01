import xml.etree.ElementTree as ET
import requests
import datetime

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
    created_date =

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

if __name__ == "__main__":

    fp = open("/Users/ikalvari/Desktop/SARS-CoV-2/release_data/corona_genomes_from_kevin/sars_accessions.txt", 'r')

    for line in fp:
        line = line.strip()
        tax_id = fetch_taxid_from_ncbi(line)

        print "%s\t%s" % (line, tax_id)

"""

	fp = open(sys.argv[1], 'r')
	rfam_taxids = [x.strip() for x in fp]
	fp.close()

	fp = open(sys.argv[2], 'r')
	accessions = [x.strip() for x in fp]
	fp.close()

	new_taxids = []
	taxid_seqacc_mappings = {}
	
	# fetch all available taxids from ncbi 
	for accession in accessions:

		taxid = fetch_taxid_from_ncbi(accession)	
		
		if taxid is not None:
			new_taxids.append(taxid)
			print taxid
			taxid_seqacc_mappings[taxid] = accession
	
	# fetch all novel taxids
	novel_taxids = list(set(new_taxids).difference(set(rfam_taxids)))

	fp = open("ncbi_novel_corona_taxids.txt", 'w')
	fp1 = open("ncbi_novel_corona_seqaccs.txt", 'w')

	for taxid in novel_tax_ids:
		fp.write(taxid+'\t')
		fp1.write(taxid_seqacc_mappings[taxid]+'\n')

	fp.close()
	fp1.close()
"""
