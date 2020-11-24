import os
import subprocess
import argparse

MEMORY = 8000
CPU = 4
LSF_GROUP = "/rfam_srch"
REQUIRED_FILES = ["SEED", "DESC", "species", "outlist", "seedoutlist"]

# ---------------------------------------------------------------------------------------------

def launch_new_rfsearch(family_dir, cpu=4):
    """
    Launches a new LSF job

    family_dir: The location of a valid family directory
    cpus: Number of CPUs to use per thread

    return: void
    """

    lsf_err_file = os.path.join(family_dir, "auto_rfsearch.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfsearch.out")

    job_name = os.path.basename(family_dir)

    cmd = ''
    if os.path.exists(os.path.join(family_dir, "DESC")) is False:
        # LSF command to be executed
        cmd = ("bsub -M %s -o %s -e %s -n %s -g %s -q production-rh74 "
               "-J %s \"cd %s && rfsearch.pl -t 25 -cnompi -q production-rh74 -relax -nodesc\"")
    else:
        cmd = ("bsub -M %s -o %s -e %s -n %s -g %s -q production-rh74 "
               "-J %s \"cd %s && rfsearch.pl -t 25 -cnompi -q production-rh74 -relax -ignoresm\"")
    
    file_extensions = ["cmsearch", "tbl", "err"]
    
    for file_type in file_extensions:
	
	subprocess.call("rm %s/*.%s" % (family_dir, file_type), shell=True)
	
    # call command
    #print (cmd % (MEMORY, lsf_out_file, lsf_err_file, cpu, LSF_GROUP, job_name, family_dir))
    subprocess.call(cmd % (MEMORY, lsf_out_file, lsf_err_file,
                          cpu, LSF_GROUP, job_name, family_dir), shell=True)

# ---------------------------------------------------------------------------------------------

def remove_all_gaps(family_dir):
	"""
	"""

	seed = os.path.join(family_dir, "SEED")
	old_seed = os.path.join(family_dir, "OLD_SEED") 
	
	# rename SEED to OLD_SEED
	os.rename(seed, old_seed)
	# regenerate SEED file removing all gaps
	cmd = "esl-reformat --mingap stockholm %s > %s" % (old_seed, seed)

	subprocess.call(cmd, shell=True)

	if os.path.exists(seed):
		return True

	return False

# ---------------------------------------------------------------------------------------------

def parse_arguments():

	"""
	"""
	parser = argparse.ArgumentParser(description='Script to precompute families')

    	required_arguments = parser.add_argument_group("required arguments")
    	required_arguments.add_argument("--input", help="A file listing all family directories", type=str)
	
	parser.add_argument('--no-gap', help='Remove all gap columns from the alignment', action="store_true", default=False)
	
	return parser

# ---------------------------------------------------------------------------------------------

if __name__ == '__main__':

	
	"""
	dirs = [x for x in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(),x))]
	err = {"MIPF0001046__mir-1012_relabelled": "", "MIPF0000411__mir-343_relabelled": "", "MIPF00000243__mir-240_relabelled": ""}
	for fam_dir in dirs:
		if fam_dir not in err:
			fam_dir_loc = os.path.join(os.getcwd(), fam_dir)
			launch_new_rfsearch(fam_dir_loc, cpu=4)
	"""

	parser = parse_arguments()
    	args = parser.parse_args()

	fp = open(args.input, 'r')
	family_dirs = [x.strip() for x in fp]
	fp.close()

	for family_dir in family_dirs:
		if args.no_gap:
			check = remove_all_gaps(family_dir)
			if check is True:
				launch_new_rfsearch(family_dir, cpu=4)

