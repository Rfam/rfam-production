import os
import json
import argparse
import urllib
import logging

# ------------------------------------------------------------------------------------------


search_dirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
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

	fp = open(desc_loc, 'r')
	desc_lines = fp.read()
	fp.close()
	
	# check if we can find the reference lines in DESC
        if desc_lines.find(REF_STRING) != -1:
		return True	

	return False
	
# ------------------------------------------------------------------------------------------

def parse_arguments():

	parser = argparse.ArgumentParser()
	
	parser.add_argument("--mirna-list", 
		help="A .json file containing all miRNAs to validate", action="store")
	parser.add_argument("--desc", help="Only perform DESC validation", 
		action="store_true", default=False)
	parser.add_argument("--svn", help="Check family exists in the SVN repository", 
		action="store_true", default=False)
	parser.add_argument("--log", help="Creates a log file with all validated DESC files", 
		action="store_true", default=False)

	return parser	

# ------------------------------------------------------------------------------------------

def get_mirna_directory_location(mirna_id):
	
	
        if mirna_id.find("_relabelled")==-1:
		dir_label = mirna_id+"_relabelled"

        for search_dir in search_dirs:
        	family_dir_loc = os.path.join(search_dir, dir_label)
                if os.path.exists(family_dir_loc):
			return family_dir_loc

	return None

# ------------------------------------------------------------------------------------------


def check_family_exists_in_svn(rfam_acc):
	
	svn_url = "https://xfamsvn.ebi.ac.uk/svn/data_repos/trunk/Families/%s"
	
	status = False
	# Check if entry existis on SVN repo; status=True

	return status

# ------------------------------------------------------------------------------------------


if __name__=='__main__':

		
	parser = parse_arguments()
	args = parser.parse_args()

	fp = open(args.mirna_list, 'r')
	mirnas = json.load(fp)
	fp.close()

	
#	if args.log is True:

	for mirna in mirnas:
		mirna_dir_loc = get_mirna_directory_location(mirna)
		if mirna_dir_loc is not None:
			if args.desc is True:
				desc_loc = os.path.join(mirna_dir_loc, "DESC")
				if os.path.exists(desc_loc):
					check = check_desc_reference_is_valid(desc_loc, REF_STRING)
					if check is False:
						print (mirna_dir_loc)
				
		
	
