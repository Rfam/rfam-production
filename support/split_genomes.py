import os
import sys
import shutil
import subprocess

from config import rfam_local as conf
from config import gen_config as gc
from utils import genome_search_utils as gsu

# ------------------------------------------------------------------------


def split_genome_to_chunks(updir, upid):
    """

    updir:
    upid:

    return:
    """
    # get updir location
    upid_fasta = os.path.join(updir, upid + '.fa')
    seq_chunks_dir = os.path.join(updir, "search_chunks")

    if not os.path.exists(seq_chunks_dir):
        os.mkdir(seq_chunks_dir)
        os.chmod(seq_chunks_dir, 0777)

        # check if we need to split the seq_file
        if gsu.count_nucleotides_in_fasta(upid_fasta) >= gc.SPLIT_SIZE:
            # split sequence file into smalled chunks
            gsu.split_seq_file(upid_fasta, gc.SPLIT_SIZE, dest_dir=seq_chunks_dir)

            # now index the fasta files
            seq_files = os.listdir(seq_chunks_dir)
            for seq_file in seq_files:
                seq_file_loc = os.path.join(seq_chunks_dir, seq_file)
                cmd = "%s --index %s" % (conf.ESL_SFETCH, seq_file_loc)
                subprocess.call(cmd, shell=True)

        # for input consistency if the sequence file is small, copy it in the
        # search_chunks directory
        else:
            # copy file
            shutil.copyfile(upid_fasta, os.path.join(seq_chunks_dir,
                                                     upid + '.fa'))
            # index file
            cmd = "%s --index %s" % (conf.ESL_SFETCH, os.path.join(seq_chunks_dir,
                                                                   upid + '.fa'))
            subprocess.call(cmd, shell=True)

# ------------------------------------------------------------------------

if __name__ == '__main__':

    project_dir = sys.argv[1]
    # this can be a file of upids or a upid string UPXXXXXXXX
    upid_input = sys.argv[2]

    if os.path.isfile(upid_input):
        fp = open(upid_input, 'r')
        upids = [x.strip() for x in fp]
        fp.close()

        for upid in upids:
            suffix = upid[-3:]
            subdir_loc = os.path.join(project_dir, suffix)
            updir_loc = os.path.join(subdir_loc, upid)

            split_genome_to_chunks(updir_loc, upid)
    else:
        # get updir location and subdir
        suffix = upid_input[-3:]
        subdir_loc = os.path.join(project_dir, suffix)
        updir_loc = os.path.join(subdir_loc, upid_input)

        split_genome_to_chunks(updir_loc, upid_input)



