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
	parser.add_argument("--acc-type", help="Rfam entry type (F: Family, M: Motif, C: Clan, G: Genome, R: Regions)",
                          type=str, choices=['F', 'M', 'C', 'G', 'R'], required=True)
	parser.add_argument("--memory", help="Memory to reserve for running the job", type=str, default='16384', required=False)

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
	mem = args.memory

	fp = open(rfam_accession_file, 'r')
	
	accession_list = [x.strip() for x in fp]

	for accession in accession_list:
		cmd = "bsub -M %s -q production-rh74 -g /rfam_xml_dumps -R \"rusage[mem=%s]\" -F 1000000 source /nfs/production/xfam/users/rfamprod/code/env2/bin/activate && export DJANGO_SETTINGS_MODULE=\'rfam_schemas.rfam_schemas.settings\' && python %s --type %s --acc %s --out %s" % (mem, mem, path_to_xml_dump, acc_type, accession, dest_dir)
		subprocess.call(cmd, shell=True)
		
	

	
