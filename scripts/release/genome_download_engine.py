import os
import sys
import subprocess

from config import gen_config as gc
from config import rfam_local as rl

# ----------------------------------------------------------------------


def launch_genome_download(project_dir, upid_list):

    fp = open(upid_list, 'r')
    upids = [x.strip() for x in fp]
    fp.close()

    if not os.path.exists(project_dir):
        os.mkdir(project_dir)
        os.chmod(project_dir, 0777)

    for upid in upids:
        # get subdir index and create if it does not exist
        subdir_idx = upid[-3:]
        subdir = os.path.join(project_dir, subdir_idx)

        if not os.path.exists(subdir):
            os.mkdir(subdir)
            os.chmod(subdir, 0777)

        prot_dir = os.path.join(subdir, upid)

        if not os.path.exists(prot_dir):
            os.mkdir(prot_dir)
            os.chmod(prot_dir, 0777)

        cmd = ("bsub -M %s "
               "-R \"rusage[mem=%s,tmp=%s]\" "
               "-o \"%s\" "
               "-e \"%s\" "
               "-u \"%s\" "
               "-n 4 "
               "-Ep \"rm -rf luigi\" "
               "-g %s "
               "python %s %s %s") % (
                  gc.MEM, gc.MEM, gc.TMP_MEM,
                  os.path.join(prot_dir, "download.out"),
                  os.path.join(prot_dir, "download.err"),
                  gc.USER_EMAIL, gc.LSF_GEN_GROUP,
                  rl.DWL_SCRIPT, prot_dir, upid)

        subprocess.call(cmd, shell=True)

# ----------------------------------------------------------------------

if __name__ == '__main__':

    project_dir = sys.argv[1]
    upid_list = sys.argv[2]

    launch_genome_download(project_dir, upid_list)