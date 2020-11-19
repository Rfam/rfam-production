import os
import sys
import json
import subprocess

if __name__ == '__main__':

    cmd = "bsub -M 10000 -R \"rusage[mem=10000]\" -o %s -e %s python /nfs/production/xfam/users/ikalvari/rfamp/code/production/rfam-production/support/restore_genome.py %s %s --uniprot"

    project_dir = sys.argv[1]
    input = sys.argv[2]

    if os.path.isfile(input):
        fp = open(input, 'r')
        accs = json.load(fp)
	fp.close()
        for upid in accs.keys():
            subdir = os.path.join(project_dir, upid[-3:])
            updir = os.path.join(subdir, upid)

            out_file = os.path.join(updir, "download.out")
            err_file = os.path.join(updir, "download.err")

            upaccs_file = os.path.join(updir, upid+"_accessions.json")
    
            job_cmd = cmd % (out_file, err_file,
                             updir, upaccs_file)

            subprocess.call(job_cmd, shell=True)
    else:

        upid = input
        subdir = os.path.join(project_dir, upid[-3:])
        updir = os.path.join(subdir, upid)

        out_file = os.path.join(updir, "download.out")
        err_file = os.path.join(updir, "download.err")

        upaccs_file = os.path.join(updir, upid + "_accessions.json")

        job_cmd = cmd % (out_file, err_file,
                         updir, upaccs_file)

        subprocess.call(job_cmd, shell=True)
