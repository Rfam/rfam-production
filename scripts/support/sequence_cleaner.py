import os
import sys
import json
import subprocess

if __name__ == '__main__':

    cmd = "bsub -M 6000 -R \"rusage[mem=6000]\" -o %s -e %s python /nfs/production/xfam/users/ikalvari/rfamp/code/production/rfam-production/support/cleanup_sequences.py %s"

    project_dir = sys.argv[1]
    upid_input = sys.argv[2]

    if os.path.isfile(upid_input):
        fp = open(upid_input, 'r')
        accs = json.load(fp)
        fp.close()

        for upid in accs.keys():
            subdir = os.path.join(project_dir, upid[-3:])
            updir = os.path.join(subdir, upid)

            out_file = os.path.join(updir, "clean.out")
            err_file = os.path.join(updir, "clean.err")

            job_cmd = cmd % (out_file, err_file, updir)
            print job_cmd
            # subprocess.call(job_cmd, shell=True)
    else:
        upid = upid_input
        subdir = os.path.join(project_dir, upid[-3:])
        updir = os.path.join(subdir, upid)

        out_file = os.path.join(updir, "clean.out")
        err_file = os.path.join(updir, "clean.err")

        job_cmd = cmd % (out_file, err_file, updir)
        print job_cmd
        #subprocess.call(job_cmd, shell=True)
