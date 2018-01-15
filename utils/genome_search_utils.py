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

"""
Genome search support functions
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import copy
import math
import subprocess
import logging
import gzip
import shutil

from config import rfam_local as rfl

# ---------------------------------GLOBALS-------------------------------------

ESL_SEQSTAT = rfl.ESL_SEQSTAT
# could import this from a configuration file
ESL_SSPLIT = rfl.ESL_SSPLIT

MB = 1000000

# -----------------------------------------------------------------------------


def split_seq_file(seq_file, size, dest_dir=None):
    """
    Splits a fasta sequence file of size X into chunks of specified size using
    Bio-Easel's esl-ssplit.pl

    seq_file (string): A string representing the path to the sequence file
    size (int): An integer specifying the size of the file chunks
    dest_dir (string): A string representing the path to the output directory
    """

    # get seq_file's size
    seq_file_size = os.path.getsize(seq_file)

    # calculate number of files to split seq_file to
    chunks_no = math.ceil(seq_file_size / size)

    # execute command
    try:
        cmd = ''
        filename = os.path.basename(seq_file).partition('.')[0]

        if dest_dir is None:
            cmd = "%s -n -r %s %s" % (ESL_SSPLIT,
                                      seq_file,
                                      str(chunks_no))
        else:
            cmd = "%s -n -r -oroot %s -odir %s %s %s" % (ESL_SSPLIT,
                                                          filename,
                                                          dest_dir,
                                                          seq_file,
                                                          str(chunks_no))
        subprocess.call(cmd, shell=True)

    except:
        raise IOError


# -----------------------------------------------------------------------------


def extract_job_stats(lsf_output_file):
    """
    Loops over the out_dir which contains all .out LSF job files, parses the files and returns job
    details such as, start and end dates, max required memory etc.

    out_dir: A directory where job .out files have been stored
    """

    gen_exec_stats = {}

    # open lsf output file and read contents
    fp = open(os.path.join(input, file), 'r')
    content = fp.readlines()
    fp.close()

    # get reference proteome id
    upid = file.partition('.')[0]
    stats = {}

    for line in content:
        if line.find("Started") != -1:
            line = line.strip().split(' ')
            stats["start_date"] = ' '.join(line[2:])

        elif line.find("Results") != -1:
            line = line.strip().split(' ')
            stats["end_date"] = ' '.join(line[3:])

        elif line.find("CPU time") != -1:
            line = line.strip().split(' ')
            stats["cpu_time"] = float(line[17])  # /60.0

        elif line.find("Max Memory") != -1:
            line = line.strip().split(' ')
            stats["mem"] = line[15]

        elif line.find("Max Processes") != -1:
            line = line.strip().split(' ')
            stats["max_proc"] = line[len(line) - 1]

        elif line.find("Max Threads") != -1:
            line = line.strip().split(' ')
            stats["max_threads"] = line[len(line) - 1]

    gen_exec_stats[upid] = stats

    return gen_exec_stats

# -----------------------------------------------------------------------------


def get_max_required_memory(lsf_output_dir):
    """

    lsf_output_dir:
    """

    max_mem = 0

    lsf_out_files = [x for x in os.listdir(lsf_output_dir) if x.endswith(".out")]

    for lsf_out_file in lsf_out_files:
        output_file_loc = os.path.join(lsf_output_dir, lsf_out_file)
        all_stats = extract_job_stats(output_file_loc)


# -----------------------------------------------------------------------------


def extract_project_stats(lsf_output_dir):
    """
    Loops over the out_dir which contains all .out LSF job files, parses the files and returns job
    details such as, start and end dates, max required memory etc.

    lsf_output_dir: A directory where job .out files have been stored
    """

    all_stats = {}

    # move all files in a single dir and get from there
    output_files = os.listdir(lsf_output_dir)

    total_exec_time = 0.0

    for output_file in output_files:
        job_stats = extract_job_stats(output_file)
        upid = job_stats.keys()
        total_exec_time = total_exec_time + float(job_stats[upid]["cpu_time"])
        all_stats[upid] = copy.deepcopy(job_stats[upid])

    avg_exec_time = total_exec_time / len(all_stats.keys())

    return all_stats, total_exec_time, avg_exec_time

# -----------------------------------------------------------------------------


def index_sequence_file(seq_file):
    """
    Uses esl-sfetch to index a sequence file. The sequence file must be in
    fasta format

    seq_file (string): A string representing the path to the sequence file

    output: An indexed file X.fa.ssi
    returns: void
    """
    esl_sfetch = ""

    # call command to index sequence file
    cmd = esl_sfetch + " --index %s" % seq_file
    subprocess.call(cmd, shell=True)

# -----------------------------------------------------------------------------


def calculate_genome_size(genome):
    """
    Uses Infernal's esl-seqstat to calculate the size of a genome

    genome: This can be either a directory containing multiple fasta files or
    a single fasta file
    returns: The size of the genome as a number of nt
    """
    # call count_nucleotides
    pass

# -----------------------------------------------------------------------------


def count_nucleotides_in_fasta(fasta_file):
    """
    Uses Infernal's esl-seqstat to get the number of nucleotides in a fasta a
    given fasta file

    param fasta_file (string): A string representing the path to a valid fasta
    file
    returns (int): The number of nucleotides in the given fasta file
    """

    # some sanity checks
    if os.path.exists(fasta_file):
        fasta_file_dir = os.path.split(fasta_file)[0]
        # decompress fasta file first
        if fasta_file.endswith(".gz"):
            filename = fasta_file.partition('.')[0]
            with gzip.open(fasta_file, 'r') as fasta_in, open(os.path.join(fasta_file_dir,
                                                                           filename+'.fa'), 'w') as fasta_out:
                shutil.copyfileobj(fasta_in, fasta_out)
            fasta_in.close()
            fasta_out.close()
            # delete compressed version
            os.remove(fasta_file)
            fasta_file = os.path.join(fasta_file_dir, filename+'.fa')

        # create a process channel and execute the command in cmd_args
        cmd_args = [ESL_SEQSTAT, '--dna', "-c", fasta_file]
        channel = subprocess.Popen(cmd_args, stdout=subprocess.PIPE)

        # fetch esl-seqstat result
        proc_output = channel.communicate()

        # process the response
        esl_seqstat_out = proc_output[0].split('\n')

        # extract the total number of residues from the response
        total_residues = 0
        for item in esl_seqstat_out:
            if item.find("Total") != -1:
                item = item.split(' ')
                total_residues = item[len(item)-1]

    else:
        sys.exit("\nInvalid input. The file provided does not exist!!\n")

    # cast to int and return
    return int(total_residues)

# -----------------------------------------------------------------------------


def get_genome_file_sizes(genome_dir, output_file=True):
    """
    Generates a list of accessions each corresponding to a sequence file in the
    genome directory the number of nucleotides in that particular file. Can be
    used during metadata extraction as well

    genome_dir (path): The path to a genome directory as organized by
    genome_downloader.py

    returns: A dictionary with the number of nucleotides per genome file
    {acc : nucleotides}
    """

    # list all sequence files - we need a sanity check here...
    genome_files = os.listdir(genome_dir)

    genome_sizes = {}

    for gen_file in genome_files:
        gen_file_loc = os.path.join(genome_dir, gen_file)

        accession = gen_file.partition('.')[0]
        file_size = count_nucleotides_in_fasta(gen_file_loc)
        genome_sizes[accession] = file_size

    if output_file is True:
        genome_sizes_fp = open(os.path.join(genome_dir, "genome_sizes.tsv"),'w')
        erroneous_files = open(os.path.join(genome_dir, "erroneous_files.tsv"),'w')

        for acc in genome_sizes.keys():
                genome_sizes_fp.write(acc + '\t' + str(genome_sizes[acc]) + '\n')
                if genome_sizes[acc] == 0:
                    erroneous_files.write(acc + '\n')

        genome_sizes_fp.close()
        erroneous_files.close()

    return genome_sizes

# -----------------------------------------------------------------------------


def calculate_seqdb_size(project_dir, mb=True):
    """
    Loops over all genome directories in the project dir, as organized by
    genome_downloader.py and calculates the size of the new seqdb

    project_dir (path): The path to the project_dir (result of genome_downloader.py)
    mb (boolean): If True convert nucleotides to megabases. Default True

    return: The size of the seqdb (nt)
    """

    seqdb_size = 0

    # list domain directories
    domain_dirs = [x for x in os.listdir(project_dir)
                   if os.path.isdir(os.path.join(project_dir, x))]

    # loop over domain directory
    for domain_dir in domain_dirs:
        domain_dir_loc = os.path.join(project_dir, domain_dir)
        genome_dirs = os.listdir(domain_dir_loc)

        # loop over genome directory
        for genome in genome_dirs:
            genome_dir_loc = os.path.join(domain_dir_loc, genome)
            # get genome seq sizes
            gen_size_dict = get_genome_file_sizes(genome_dir_loc, output_file=False)

            # sum genome seq sizes and seqdb_size
            for acc in gen_size_dict.keys():
                seqdb_size += gen_size_dict[acc]

    # convert seqdb_size to Megabases and return
    if mb is True:
        return float(seqdb_size)/float(MB)

    return seqdb_size

# -----------------------------------------------------------------------------


if __name__ == '__main__':

   pass