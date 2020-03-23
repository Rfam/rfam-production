import os
import sys
import requests
import xml.etree.ElementTree as ET


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

# test

# ------------------------------------------------------------------------------

if __name__=="__main__":

	fp = open(sys.argv[1], 'r')
	rfam_taxids = [x.strip() for x in fp]
	fp.close()

	fp = open(sys.argv[2], 'r')
	accessions = [x.strip() for x in fp]
	fp.close()

	new_taxids = []

	# fetch all available taxids from ncbi 
	for accession in accessions:

		taxid = fetch_taxid_from_ncbi(accession)	
		
		if taxid is not None:
			new_taxids.append(taxid)

	# fetch all novel taxids
	novel_taxids = list(set(new_taxids).difference(set(rfam_taxids)))

	print len(novel_taxids)
