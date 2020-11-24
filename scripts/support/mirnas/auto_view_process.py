import os
import sys
import subprocess

MEMORY = 20000
CPU = 4
LSF_GROUP = "/family_srch"

# ---------------------------------------------------------------------------


def run_view_processes(rfam_acc, uuid, log_dir):

	cmd = ("bsub -M %s -o %s -e %s -n %s -g %s -q production-rh74 "
		"-J %s \"rfam_family_view.pl -id %s -f %s family\"")	

	err_file = os.path.join(log_dir, rfam_acc+'.err')
	out_file = os.path.join(log_dir, rfam_acc+'.out')

	subprocess.call(cmd %(MEMORY, out_file, err_file, CPU, 
			LSF_GROUP,rfam_acc, uuid, rfam_acc), shell=True)


# ---------------------------------------------------------------------------


if __name__=='__main__':

	log_dir = sys.argv[2]
	id_file = sys.argv[1]

	fp = open(id_file, 'r')
	
	for line in fp:
		line = line.strip().split("\t")
		run_view_processes(line[0], line[1], log_dir)

	fp.close()
