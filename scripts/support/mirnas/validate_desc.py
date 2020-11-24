import os
import json
import argparse
import subprocess

from subprocess import Popen, PIPE

# ------------------------------------------------------------------------------------------


searchdirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]

MEMORY = 8000
CPU = 4
LSF_GROUP = "/family_srch"


REF_STRING = """RN   [1]
RM   30423142
RT   miRBase: from microRNA sequences to function.
RA   Kozomara A, Birgaoanu M, Griffiths-Jones S;
RL   Nucleic Acids Res. 2019;47:D155."""

# ------------------------------------------------------------------------------------------


def check_desc_reference_is_valid(desc_loc, ref_string):

	process = Popen(["grep", ref_string, desc_loc], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = process.communicate()[1]

        if output.find(ref_string) == -1:
                return False

	return True

# ------------------------------------------------------------------------------------------

def parse_arguments():

	parser = argparse.ArgumentParser()
	
	parser.add_argument("--mirna-list", help="A .json file containing all miRNAs to validate", action="store")
	parser.add_argument("--log", help="Creates a log file with all validated DESC files", action="store_true", default=False)

	return parser	

# ------------------------------------------------------------------------------------------

def get_mirna_directory_location(mirna_id):
	
	
        if accession.find("_relabelled")==-1:
		dir_label = accession+"_relabelled"

        for search_dir in search_dirs:
        	family_dir_loc = os.path.join(search_dir, dir_label)
                if os.path.exists(family_dir_loc):
			return family_dir_loc

	return None

# ------------------------------------------------------------------------------------------

if __name__=='__main__':

	parser = parse_arguments()
	args = parser.parse_args()

	fp = open(args.mirna_list, 'r')
	mirnas = json.load(fp)
	fp.close()

	for mirna in mirnas:
		mirna_dir_loc = get_mirna_directory_location(mirna_id)
		if mirna_dir_loc is not None:
			desc_loc = os.path.join(mirna_dir_loc, "DESC")
			if os.path.exists(desc_loc):
				check = check_desc_reference_is_valid(desc_loc, REF_STRING)
				if check is False:
					print (mirna_dir_loc)
					
		
	


