import os
import sys
import argparse
import shutil
import subprocess

# ---------------------------- Variable initialization ---------------------------

MEMORY = 2000
CPU = 8
LSF_GROUP = "/family_srch"

# --------------------------------------------------------------------------------

def check_seed_format(seed_alignment, format="stockholm"):
	"""
	Checks if a SEED alignment is in a specific format

	return: void
	"""
	
	pass

# --------------------------------------------------------------------------------

def create_family_directory(dirname, seed_alignment, dest_dir):
	"""
	Creates a new Rfam family directory and sets all the files required
	to launch rfsearch.

	dirname: A name for the new Rfam family in the format str_str
	seed_alignment: A valid seed alignment in stockholm format
	dest_dir: The base directory where the new family directory will be 
	created

	return: void
	"""
	
	# check if destination directory existis
	if not os.path.exists(dest_dir):
		sys.exit("\nDestination directory does not exist!\n")

	# go on with creating a new family directory
	family_dir = os.path.join(dest_dir, dirname)
	os.mkdir(family_dir)
	
	shutil.copy(seed_alignment, os.path.join(family_dir, "SEED"))
	
# --------------------------------------------------------------------------------

def launch_new_rfsearch(family_dir):
	"""
	family_dir: The location of a valid family directory
	
	return: void 
	"""

	lsf_err_file = os.path.join(family_dir, "auto_rfsearch.err")
	lsf_out_file = os.path.join(family_dir, "auto_rfsearch.out")

	# LSF command to be executed
	cmd = ("bsub -M %s -R \"rusage[mem=%s]\" -o %s -e %s -n %s -g %s -q production-rh74 "
          "-J %s \"cd %s && rfsearch.pl -t 30 -cnompi -q production-rh74 -relax\"")

	# call command
	subprocess.call(cmd % (MEMORY, MEMORY, lsf_out_file, lsf_err_file,
                         CPU, LSF_GROUP, rfam_acc, family_dir), shell=True)
	

# --------------------------------------------------------------------------------

def parse_arguments():
	"""
	Does basic command line argument parsing 
	
	return: Argparse object
	"""
	
	parser = argparse.ArgumentParser(description='Rfam family Auro-Builder')
	# group required arguments together

	req_args = parser.add_argument_group("required arguments")
	req_args.add_argument('--input', help='a directory with multiple SEEDs or a single SEED',
				type=str, required=True)
	req_args.add_argument('--dest-dir', help='destination directory where create new \
				family directories', type=str, required=True)

	return parser

# --------------------------------------------------------------------------------	


if __name__ == '__main__':

	parser = parse_arguments()
	args = parser.parse_args()
	
	if (os.path.isdir(args.input)):

		seeds = os.listdir(args.input)
	else:
		seeds = [args.input]
	
	for seed in seeds:
		dirname = seed.partition('.')[0]
		create_family_directory(dirname, seed, args.dest_dir) 

