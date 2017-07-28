#!/usr/bin/python
"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import shutil
import subprocess
import sys

from config import gen_config as gc
from scripts.validation import genome_search_validator as gsv
from utils import genome_search_utils as gsu

# -----------------------------------------------------------------------------

SPLIT_SIZE = 5427083
SRCH_MEM = 36000
SCAN_MEM = 36000
RFAM13_SIZE = 451031.997884  # size of rfamseq13 in Mb
CM_NO = 2588  # number of cms in Rfam.cm file
CPU_NO = 5
SGROUP_IDX = 65  # subgroup index - Ascii for A
LSF_GROUP = "/rfam_srch_mpi/%s"
SUB_DIRS = 26  # max 26 subdirs as the number of the alphabet

# filtered - add this as a command line option

CMD_TEMPLATE_MPI = ("bsub -q mpi-rh7 -M %s -R \"rusage[mem=%s]\" -o %s -e %s -n %s -g %s -R \"span[hosts=1]\" "
                    "-a openmpi mpiexec -mca btl ^openib -np %s "
                    "%s -o %s --tblout %s --acc --cut_ga --rfam --notextw --nohmmonly -Z %s --mpi %s %s")

"""
# unfiltered
CMD_TEMPLATE_MPI = ("bsub -q mpi-rh7 -M %s -R \"rusage[mem=%s]\" -o %s -e %s -n %s -g %s -R \"span[hosts=1]\" "
                    "-a openmpi mpiexec -mca btl ^openib -np %s "
                    "%s -o %s --tblout %s --acc --notextw --nohmmonly -Z %s --mpi %s %s")
"""


# -------------------------------------------------------------------------

def calculate_required_memory():
    """

    :return:
    """

    memory = 0

    return memory


# -------------------------------------------------------------------------

def gather_genome_search_results(input_dir, dest_dir=None):
    """
    Gather and merge results from input_dir perhaps a recursive function to
    loop over all subdirs eg copytree etc look into available python
    libraries
    input_dir:
    dest_dir:
    return:
    """

    memory = 0

    # to be implemented

    return memory


# -------------------------------------------------------------------------

def genome_scan_advanced(input_dir, dest_dir, tool="cmsearch"):
    """
    using cmsearch by default unless defined otherwise

    :param input_dir:
    :param dest_dir:
    :param tool:
    :return:
    """

    # initialize output space
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir, 0775)

    # create subdirs
    i = 0
    while i < SUB_DIRS:
        subdir = os.path.join(dest_dir, chr(SGROUP_IDX + i))
        if not os.path.exists(subdir):
            os.mkdir(subdir)
            os.chmod(subdir, 0775)
        i += 1

    # get candidate sequence files
    seq_files = [x for x in os.listdir(input_dir)
                 if x.endswith('.fa') or x.endswith('.fasta')]

    out_idx = 0

    for seq_file in seq_files:
        # 1. get size
        seq_file_loc = os.path.join(input_dir, seq_file)
        filename = seq_file.partition('.')[0]

        gen_output_dir = os.path.join(os.path.join(dest_dir, chr(SGROUP_IDX + out_idx)),
                                      filename)

    if not os.path.exists(gen_output_dir):
        os.mkdir(gen_output_dir, 0775)

        gen_input_dir = ''
        if gsu.count_nucleotides_in_fasta(seq_file_loc) >= SPLIT_SIZE:
            gen_input_dir = os.path.join(input_dir, filename)

            # create a distinct directory for the genome
            if not os.path.exists(gen_input_dir):
                os.mkdir(gen_input_dir, 0775)

            shutil.move(seq_file_loc, os.path.join(gen_input_dir, os.path.basename(seq_file_loc)))

            # new input file location
            seq_file_loc = os.path.join(gen_input_dir, os.path.basename(seq_file_loc))

            # split sequence file into smalled chunks and store under destination directory
            gsu.split_seq_file(seq_file_loc, SPLIT_SIZE, dest_dir=gen_input_dir)

            # list all smaller files
            genome_chunks = [x for x in os.listdir(gen_input_dir) if x.find(".fa.") != -1]

            group_idx = out_idx

            # group_idx = 0
            for genome_chunk in genome_chunks:
                # index all sequence files
                chunk_loc = os.path.join(gen_input_dir, genome_chunk)

                gsu.index_sequence_file(chunk_loc)
                chunk_name = genome_chunk
                lsf_out_file = os.path.join(gen_output_dir, chunk_name + ".out")
                lsf_err_file = os.path.join(gen_output_dir, chunk_name + ".err")
                inf_tbl_file = os.path.join(gen_output_dir, chunk_name + ".tbl")
                inf_out_file = os.path.join(gen_output_dir, chunk_name + ".inf")

                if tool == 'cmsearch':
                    cmd = CMD_TEMPLATE_MPI % (SRCH_MEM, SRCH_MEM, lsf_out_file, lsf_err_file,
                                              CPU_NO, LSF_GROUP % (chr(SGROUP_IDX + group_idx)),
                                              CPU_NO, gc.CMSEARCH, inf_out_file,
                                              inf_tbl_file, RFAM13_SIZE,
                                              gc.CMFILE, chunk_loc)
                else:
                    cmd = CMD_TEMPLATE_MPI % (SRCH_MEM, SRCH_MEM, lsf_out_file, lsf_err_file,
                                              CPU_NO, LSF_GROUP % (chr(SGROUP_IDX + group_idx)),
                                              CPU_NO, gc.CMSCAN, inf_out_file,
                                              inf_tbl_file, RFAM13_SIZE,
                                              gc.CMFILE, chunk_loc)

                subprocess.call(cmd, shell=True)

        # small enough genomes that don't need splitting
        else:
            gsu.index_sequence_file(seq_file_loc)

            lsf_out_file = os.path.join(gen_output_dir, filename + ".out")
            lsf_err_file = os.path.join(gen_output_dir, filename + ".err")
            inf_tbl_file = os.path.join(gen_output_dir, filename + ".tbl")
            inf_out_file = os.path.join(gen_output_dir, filename + ".inf")

            cmd = CMD_TEMPLATE_MPI % (SRCH_MEM, SRCH_MEM, lsf_out_file, lsf_err_file,
                                      CPU_NO, LSF_GROUP % (chr(SGROUP_IDX + out_idx)),
                                      CPU_NO, gc.CMSEARCH, inf_out_file,
                                      inf_tbl_file, RFAM13_SIZE,
                                      gc.CMFILE, seq_file_loc)

            subprocess.call(cmd, shell=True)

        out_idx += 1
        if out_idx == SUB_DIRS:
            out_idx = 0


# -------------------------------------------------------------------------


def restore_io_paths(input_dir, output_dir):
    """
    Restore searches from and to the same directories.
    ** Will need to simplify this

    input_dir: The path to the input directory as organized by genome_scanner
    output_dir: The path to the output directory of the run we want to restore

    return: A list of tuples with input and output paths
    """

    io_path_list = []
    err_cases = gsv.check_search_err_files(output_dir)

    # loop over all output subdirs
    for subdir in err_cases.keys():
        dest_dir = os.path.join(output_dir, subdir)

        # get genome i/o directories
        for genome in err_cases[subdir].keys():
            out_gen_dir = os.path.join(dest_dir, genome)
            in_gen_dir = os.path.join(input_dir, genome)

            # get sequence file i/o paths -
            # need to handle the case in which we don't have a dir,
            # or just create one for all cases in input dir
            for err_file in err_cases[subdir][genome]:
                # add a check here to look for an index number
                out_file_loc = os.path.join(out_gen_dir, err_file)
                in_file_loc = os.path.join(in_gen_dir, err_file)
                io_path_list.append((in_file_loc, out_file_loc))

                # remove files from previous runs

    return io_path_list


# -------------------------------------------------------------------------

def restore_jobs_with_multi_cms(cm_dir, input_dir, output_dir):
    """
    Restore search jobs by scanning using smaller cm files

    cm_dir: The path to the CM directory. If None, use default
    input_dir: The path to the input directory as organized by genome_scanner
    output_dir: The path to the output directory of the run we want to restore
    returns: void
    """

    io_path_list = restore_io_paths(input_dir, output_dir)

    cms = [x for x in os.listdir(cm_dir) if x.find('.') == -1]

    for cm_file in cms:
        cm_file_loc = os.path.join(cm_dir, cm_file)

        for err_case in io_path_list:
            seq_file = err_case[0]
            out_file = err_case[1]

            # cleanup old files
        if os.path.exists(out_file + ".out"):
            os.remove(out_file + ".out")
        if os.path.exists(out_file + ".err"):
            os.remove(out_file + ".err")
        if os.path.exists(out_file + ".tbl"):
            os.remove(out_file + ".tbl")
        if os.path.exists(out_file + ".inf"):
            os.remove(out_file + ".inf")

        lsf_out_file = out_file + '_' + cm_file + ".out"
        lsf_err_file = out_file + '_' + cm_file + ".err"
        inf_tbl_file = out_file + '_' + cm_file + ".tbl"
        inf_out_file = out_file + '_' + cm_file + ".inf"

        lsf_subgroup = out_file.split('/')[-3]

        cmd = CMD_TEMPLATE_MPI % (SRCH_MEM, SRCH_MEM, lsf_out_file,
                                  lsf_err_file, CPU_NO,
                                  LSF_GROUP % lsf_subgroup,
                                  CPU_NO, gc.CMSEARCH,
                                  inf_out_file, inf_tbl_file,
                                  RFAM13_SIZE, cm_file_loc,
                                  seq_file)

        subprocess.call(cmd, shell=True)  # re-initialization

        lsf_out_file = ''
        lsf_err_file = ''
        inf_tbl_file = ''
        inf_out_file = ''


# -------------------------------------------------------------------------


if __name__ == '__main__':

    input_dir = sys.argv[1]
    dest_dir = sys.argv[2]

    # restore searches if --restore option is provided
    if '--restore' in sys.argv:
        cm_dir = ''  # provide this as an option or add in a config file
        # cm_dir = "/nfs/production/xfam/rfam/RELEASE_FILES/RFAM_12_2/cm_groups"
        restore_jobs_with_multi_cms(cm_dir, input_dir, dest_dir)

    else:
        genome_scan_advanced(input_dir, dest_dir, tool="cmsearch")
