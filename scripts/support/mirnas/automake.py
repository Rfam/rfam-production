import os
import json
import subprocess
import argparse


# -------------------------------------------------------------------------------------------------------


searchdirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches", 
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]

MEMORY = 15000
CPU = 4
LSF_GROUP = "/family_srch"

# -------------------------------------------------------------------------------------------------------


def autorfmake(thresholds_file, serial=False):
    
    cmd = ("bsub -M %s -o %s -e %s -n %s -g %s -q production-rh74 "
               "-J %s \"cd %s && rfmake.pl -t %s -forcethr -a\"")

    serial_cmd = "rfmake.pl -t %s -forcethr"

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

		if serial is True:
			os.chdir(family_dir)
			subprocess.call(serial_cmd % str(thresholds[family]), shell=True)

		else:
			lsf_err_file = os.path.join(family_dir, "auto_rfmake.err")
                	lsf_out_file = os.path.join(family_dir, "auto_rfmake.out")
                	job_name = family
  		
			subprocess.call(cmd % (MEMORY, lsf_out_file, lsf_err_file, CPU, LSF_GROUP, job_name, family_dir, str(thresholds[family])), shell=True)
            else:
                continue
    fp.close()

# -------------------------------------------------------------------------------------------------------


def parse_arguments():

	parser = argparse.ArgumentParser()
	
	parser.add_argument("--thresholds", help="A json file with miRNA : threshold pairs", action="store")
	parser.add_argument("--serial", help="Serial execution of rfmake", action="store_true", default=False)
	
	return parser

# -------------------------------------------------------------------------------------------------------

if __name__=='__main__':


	parser = parse_arguments()
	args = parser.parse_args()

	json_file = args.thresholds

	autorfmake(json_file, args.serial)
