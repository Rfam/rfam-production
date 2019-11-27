import os
import sys
import subprocess
from config import rfam_local as rl

if __name__=='__main__':

    project_dir = sys.argv[1]
    upid_list = sys.argv[2]

    fp = open(upid_list, 'r')
    upids = [x.strip() for x in fp]
    fp.close()

    cmd = "bsub -M 6000 -R \"rusage[mem=6000]\" -g /rfam_search python %s %s %s"

    for upid in upids:
        sub_cmd = cmd % (rl.SCANNER, project_dir, upid)
        subprocess.call(sub_cmd, shell=True)
