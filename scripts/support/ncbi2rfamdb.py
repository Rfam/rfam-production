import xml.etree.ElementTree as ET

import requests


def fetch_taxid_from_ncbi(accession):
    """

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
