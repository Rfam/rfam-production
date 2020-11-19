import os
import sys
import json
import subprocess

# ------------------------------------------------------------------------


def bsub_redownload_genome_from_gca_report(updir, gca_report_file):
    """
    Launch the genome downloads on LSF

    :return:
    """

    cmd = "bsub -M 6000 -R \"rusage[mem=6000]\" -o %s -e %s python %s %s %s"

    bsub_cmd = cmd % (os.path.join(updir, 'restore.out'),
                      os.path.join(updir, 'restore.err'),
                      os.path.join(os.getcwd(), "restore_genome.py"),
                      updir, gca_report_file)

    subprocess.call(bsub_cmd, shell=True)
    #print bsub_cmd
# ------------------------------------------------------------------------

if __name__ == '__main__':

    upid_list = sys.argv[1]
    project_dir = sys.argv[2]

    # load upids
    fp = open(upid_list, 'r')
    upids = [x.strip() for x in fp]
    fp.close()

    gca_acc_fp = open(os.path.join(project_dir, "upid_gca_dict.json"), 'r')
    gca_dict = json.load(gca_acc_fp)
    gca_acc_fp.close()

    for upid in upids:
        suffix = upid[-3:]
        subdir = os.path.join(project_dir, suffix)
        updir = os.path.join(subdir, upid)
        gca_acc = str(gca_dict[upid]["GCA"])
        gca_report_file = os.path.join(updir, gca_acc+"_sequence_report.txt")

        if os.path.exists(gca_report_file):
            bsub_redownload_genome_from_gca_report(updir, gca_report_file)
        else:
            print upid




