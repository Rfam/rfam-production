import os
import sys
import subprocess

from config import gen_config as gc
from scripts.export.genomes import genome_fetch as gf

"""
A collection of functions to validate files with a new genome download
"""


# -----------------------------------------------------------------------------


def validate_sequence_file(seq_file, seq_type='dna'):
    """
    Validating sequence file file using esl-seqstat

    seq_file: The sequence file to validate
    seq_type: The type of the sequences in the file (e.g. dna, rna, amino)

    return: True if valid, False if invalid
    """

    # command string should look like esl-seqstat --seq_type seq_file
    cmd_args = []
    esl_tool = os.path.join(gc.LSF_RFAM_BIN, 'esl-seqstat')
    seq_type_arg = "--%s" % seq_type

    cmd_args = [esl_tool, seq_type_arg, seq_file]
    popen_obj = subprocess.Popen(cmd_args, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)

    message = popen_obj.communicate()

    if message[0].find("failed") == -1:
        return True

    return False

# -----------------------------------------------------------------------------


def check_compressed_file(filename):
    """
    Checks if the provided file is in one of the compressed formats

    filename: The path to input file
    returns: Boolean - True if the file is compressed, False otherwise
    """

    magic_dict = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        "\x50\x4b\x03\x04": "zip"
    }

    max_len = max(len(x) for x in magic_dict)

    with open(filename) as fp_in:
        file_start = fp_in.read(max_len)

    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return True  # can also return filetype

    fp_in.close()

    return False


# -----------------------------------------------------------------------------


def check_genome_download_status(lsf_out_file):
    """
    Opens LSF output file and checks whether the job's status is success

    lsf_out_file: LSF platform's output file generated by -o option
    returns: status 1 if the download was successful, otherwise 0
    """

    infile_fp = open(lsf_out_file, 'r')

    status = 0

    for line in infile_fp:
        if line.find("Success") != -1:
            status = 1

    infile_fp.close()

    return status

# -----------------------------------------------------------------------------


def check_all_files_downloaded(gca_report_file, genome_dir):
    """
    Parse genome GCA file and check that all files have been downloaded and
    report any missing files

    gca_report_file: A genome assembly report file provided by ENA
    genome_dir: The path to a specific genome directory

    return: True if all files were downloaded, alternatively a list of
    the accessions of the missing files
    """

    missing_files = []
    gca_accessions = gf.assembly_report_parser(gca_report_file)

    downloaded_files = [x for x in os.listdir(genome_dir) if x.endswith(".fa")]

    for accession in gca_accessions:
        file_path = os.path.join(genome_dir, accession + ".fa")

        if not os.path.exists(file_path):
            missing_files.append(accession)

    if len(missing_files) > 0:
        return missing_files

    return True
        
# -----------------------------------------------------------------------------

if __name__ == '__main__':

    pass