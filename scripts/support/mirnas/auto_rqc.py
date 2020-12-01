import os
import sys
import argparse
import subprocess
import json

from datetime import date
from subprocess import Popen, PIPE

search_dirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]

# ---------------------------------------------------------------------------------------------


def check_desc_ga(DESC, cut_ga):
	"""
	"""

	process = Popen(['grep', "GA", DESC], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = process.communicate()
	
        if output.find("%.2f"%float(cut_ga)) == -1:
                return False

        return True 
	

# ---------------------------------------------------------------------------------------------


def check_family_passes_qc(family_dir):
	
	dir_elements = os.path.split(family_dir)
	search_dir = dir_elements[0]

	os.chdir(search_dir)

	process = Popen(["rqc-all.pl", dir_elements[1]], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = process.communicate()[1]

    	if output.find("Family passed with no serious errors") == -1:
        	return False

    	return True 

# ---------------------------------------------------------------------------------------------


def find_rqc_error(rqc_output):

	error_types = {"CM": 0, "FORMAT": 0, "OVERLAP": 0, "STRUCTURE": 0,
			"MISSING": 0, "SEQUENCE": 0, "NON-CODING": 0}
	
	rqc_lines = [x.strip() for x in rqc_output.split("\n") if x!='']
	
	error_type_count = len(error_types.keys())

	success_string = "--%s check completed with no major errors"

	for error_type in error_types.keys():
		for rqc_line in rqc_lines:
			if rqc_line == success_string % error_type:
				error_types[error_type] = 1
				break

	return error_types

# ---------------------------------------------------------------------------------------------


def fetch_rqc_output(family_dir):

	dir_elements = os.path.split(family_dir)
        parent_dir = dir_elements[0]

        os.chdir(parent_dir)

	process = Popen(["rqc-all.pl", dir_elements[1]], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = process.communicate()[1]

	return output

# ---------------------------------------------------------------------------------------------


def generate_rqc_report(rqc_error_types):

	# develop functionality here 
	pass

# ---------------------------------------------------------------------------------------------


def parse_arguments():
	
	parser = argparse.ArgumentParser()

	parser.add_argument("--mirna-ids", help="A .json file with miRNAs to commit", action="store")
	parser.add_argument("--log-dir", help="Log file destination", action="store", default=os.getcwd())
	
	return parser

# ---------------------------------------------------------------------------------------------


if __name__=='__main__':


	parser = parse_arguments()
	args = parser.parse_args()

	fp = open(args.mirna_ids, 'r')	
	miRNA_accessions = json.load(fp)
	fp.close()


	#fp = open(os.path.join(args.log_dir, 'successful_mirna_commits.log'), 'w')

	for accession in miRNA_accessions.keys():
		dir_label = ''
		if accession.find("_relabelled")==-1:
			dir_label = accession+"_relabelled"
			
		for search_dir in search_dirs:
			family_dir_loc = os.path.join(search_dir, dir_label)
			if os.path.exists(family_dir_loc):
				#print family_dir_loc
				rqc_output = fetch_rqc_output(family_dir_loc)
				#print find_rqc_error(rqc_output)
