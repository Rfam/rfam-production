import os
import sys
import json
import argparse
import subprocess


# ------------------------------------------------------------------------------------


searchdirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches", 
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]

MEMORY = 8000
CPU = 4
LSF_GROUP = "/family_srch"


# ------------------------------------------------------------------------------------


def auto_addref(thresholds_file, reference):
    
    cmd = ("bsub -M %s -o %s -e %s -n %s -g %s -q production-rh74 "
               "-J %s \"cd %s && add_ref.pl %s\"")

    fp = open(thresholds_file, 'r')
    thresholds = json.load(fp)

    for family in thresholds.keys():
        for searchdir in searchdirs:
            family_dir = ""

            if family.find("relabelled") == -1:
                family_dir = os.path.join(searchdir, family + "_relabelled")
            else:
                family_dir = os.path.join(searchdir, family)

            if os.path.exists(family_dir):
                lsf_err_file = os.path.join(family_dir, "auto_add_ref.err")
                lsf_out_file = os.path.join(family_dir, "auto_add_ref.out")
		job_name = family

  		subprocess.call(cmd % (MEMORY, lsf_out_file, lsf_err_file, CPU, LSF_GROUP, job_name, family_dir, reference), shell=True)
            else:
                continue


# ------------------------------------------------------------------------------------


def add_ref_sequentially(thresholds_file, reference):

	cmd = ("add_ref.pl %s" % reference)

	fp = open(thresholds_file, 'r')
	thresholds = json.load(fp)
	fp.close()

	for family in thresholds.keys():
        	for searchdir in searchdirs:
            	family_dir = ""

            	if family.find("relabelled") == -1:
                	family_dir = os.path.join(searchdir, family + "_relabelled")
            	else:
                	family_dir = os.path.join(searchdir, family)

            	if os.path.exists(family_dir):
			# move to family directory
			os.chdir(family_dir)
			subprocess.call(cmd, shell=True)
		
# ------------------------------------------------------------------------------------


def parse_arguments():

	parser = argparse.ArgumentParser()
	
	parser.add_argument("--mirna-list", help="A json file with all miRNA candidates", action="store")
	parser.add_argument("--ref", help="A string indicating the PubMed id to use for reference", action="store", default="30423142")

	return parser

# ------------------------------------------------------------------------------------


if __name__=='__main__':


	parser = parse_arguments()
	args = parser.parse_args()
	
	auto_addref(args.mirna_list, args.ref)
