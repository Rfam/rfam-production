import os
import sys
from subprocess import Popen, PIPE

# --------------------------------------------------------------------------------------------------

def get_job_run_time(lsf_output_file, time='s'):
	"""
	"""

	fp = open(lsf_output_file, 'r')
	process = Popen(['grep', 'Run time', lsf_output_file], stdin=PIPE, stdout=PIPE, stderr=PIPE)
	output, err = process.communicate()
	
	run_time = int(output.split(" ")[-2])	

	if time == 'm':
		run_time = run_time/60	

	return run_time

# --------------------------------------------------------------------------------------------------
	
def get_job_max_used_memory(lsf_output_file):
	"""
	"""
	fp = open(lsf_output_file, 'r')
        process = Popen(['grep', 'Max Memory', lsf_output_file], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = process.communicate()

        max_memory = int(output.split(" ")[-2])

	return max_memory

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

	source_dir = sys.argv[1]
	families = [x for x in os.listdir(source_dir) if os.path.isdir(os.path.join(source_dir, x))]

	for family in families:
		family_dir = (os.path.join(source_dir, family))
		lsf_output_file = os.path.join(family_dir, "auto_rfsearch.out")
		run_time = get_job_run_time(lsf_output_file, time='m')
		memory = get_job_max_used_memory(lsf_output_file)
 
		print "%s\t%s\t%s" % (family, run_time, memory)

