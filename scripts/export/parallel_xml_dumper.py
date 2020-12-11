"""
The purpose of this script is to parallelize the launching process of xml_dumps 
in order to speed up Rfam exports. It uses rfam_xml_dumper and submits an individual 
lsf job for each entry provided in the input file  
"""

# ----------------------------------------------------------------------------------------

import os
import subprocess
import argparse

# ----------------------------------------------------------------------------------------

def parse_arguments():

	parser = argparse.ArgumentParser()
	
	parser.add_argument("--accessions", help="Rfam list of accessions to index", 
			action="store", required=True)
	parser.add_argument("--acc-type", help="rfam entry type (F: Family, M: Motif, C: Clan, G: Genome, R: Regions)",
                          type=str, choices=['F', 'M', 'C', 'G', 'R'], required=True)

	parser.add_argument("--dest-dir", help="Destination directory for the output", 
			action="store", required=True)

	return parser

# ----------------------------------------------------------------------------------------


if __name__ == '__main__':

	parser = parse_arguments()
	args = parser.parse_args()

	path_to_xml_dump = os.path.join(os.getcwd(), "rfam_xml_dumper.py")

	rfam_accession_file = args.accessions
	acc_type = args.acc_type
	dest_dir = args.dest_dir

	fp = open(rfam_accession_file, 'r')
	
	accession_list = [x.strip() for x in fp]

	for accession in accession_list:
		cmd = "bsub -M 16384 -q production-rh74 -g /rfam_xml_dumps -R \"rusage[mem=16384]\" -F 1000000 source /nfs/production/xfam/users/rfamprod/code/env2/bin/activate && export DJANGO_SETTINGS_MODULE=\'rfam_schemas.rfam_schemas.settings\' && python %s --type %s --acc %s --out %s" % (path_to_xml_dump, acc_type, accession, dest_dir)
		subprocess.call(cmd, shell=True)
		
	

	
