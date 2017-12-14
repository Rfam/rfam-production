import os
import sys
import gzip
import subprocess

# -----------------------------------------------------------------------------------------------------------
# TODO expand this to work with project directory and possibly wrap it up in a luigi pipeline for the merge to be executed in parallel

GUNZIP_CMD = "gunzip %s"

# -----------------------------------------------------------------------------------------------------------


def merge_genome_files(upid_dir):
    """
    Merge all sequence files of a genome in a single file

    upid_dir: The path to a genome directory

    :return:
    """

    sequence_dir_loc = os.path.join(upid_dir, "sequences")

    seq_dir_contents = os.listdir(sequence_dir_loc)

    upid = os.path.split(upid_dir)[1]

    # open a new sequence file for the genome
    genome_fasta = open(os.path.join(upid_dir, upid + '.fa'), 'w')

    # check if it's divided in subdirs
    if os.path.isfile(os.path.join(sequence_dir_loc, seq_dir_contents[0])):
        for seq_file in seq_dir_contents:
            seq_file_loc = os.path.join(sequence_dir_loc, seq_file)

            seq_file_fp = None
            # decompress file
            if seq_file.endswith(".gz"):
                seq_file_fp = gzip.open(seq_file_loc, 'rb')
            else:
                seq_file_fp = open(seq_file_loc, 'r')

            for line in seq_file_fp:
                genome_fasta.write(line)

            seq_file_fp.close()

    # work with multiple subdirectories
    else:
        for subdir in seq_dir_contents:
            subdir_loc = os.path.join(sequence_dir_loc, subdir)
            seq_files = os.listdir(subdir_loc)

            for seq_file in seq_files:
                seq_file_loc = os.path.join(subdir_loc, seq_file)

                # decompress file
                seq_file_fp = None
                if seq_file.endswith(".gz"):
                    seq_file_fp = gzip.open(seq_file_loc, 'rb')
                else:
                    seq_file_fp = open(seq_file_loc, 'r')

                seq_file_fp = open(seq_file_loc, 'r')

                for line in seq_file_fp:
                    genome_fasta.write(line)

                seq_file_fp.close()

    genome_fasta.close()

# -----------------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    updir = sys.argv[1]

    merge_genome_files(updir)

