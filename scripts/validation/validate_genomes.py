import os
import sys
import json
import copy
import subprocess

from config import rfam_config as rc


# -----------------------------------------------------------------------------------------------------------


def check_genome_download_status(lsf_out_file):
    """
    Opens LSF output file and checks whether the job's status is success

    lsf_out_file: LSF platform's output file generated by -o option
    returns: status 1 if the download was successful, otherwise 0
    """

    infile_fp = open(lsf_out_file, 'r')

    status = False

    for line in infile_fp:
        if line.find("Success") != -1:
            status = True

    infile_fp.close()

    return status

# -----------------------------------------------------------------------------------------------------------


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


# -----------------------------------------------------------------------------------------------------------


def check_file_format(seq_file):
    """
    Performs some sanity checks on the sequence file. Checks if file is
    compressed and if not validates the format using esl-seqstat. It will also
    check if the sequence file provided is empty or not

    seq_file: The path to a valid sequence file
    returns: True if file passed validation checks, False otherwise
    """

    status = True

    # compressed file
    if seq_file.endswith(".gz"):
        if not os.path.exists(seq_file):
            return False
        else:
            return check_compressed_file(seq_file)

    # uncompressed fasta format
    elif seq_file.endswith(".fa"):

        # check that file exists
        if not os.path.exists(seq_file):
            return False

        # check content
        else:
            cmd_args = [rc.ESL_SEQSTAT, '--informat', "fasta", "--dna", seq_file]
            channel = subprocess.Popen(cmd_args, stdout=subprocess.PIPE)

            # fetch esl-seqstat result
            proc_output = channel.communicate()

            # process the response
            esl_seqstat_out = proc_output[0].split('\n')
            # check only first line of response
            if esl_seqstat_out[0].find("Parse failed") != -1 \
                or esl_seqstat_out[0].find("Format") == -1:
                return False

    # check the size of the file
    if os.path.getsize(seq_file) == 0:
        return False

    return status

# -----------------------------------------------------------------------------------------------------------


def check_seq_file(upid, seq_file):

    err_messages = []
    # check if fasta file exists
    accession = os.path.basename(seq_file).partition('.')[0]

    if os.path.exists(seq_file):
        # we may want to delete the incorrect files to clean up the sequence directory
        if not check_file_format(seq_file):
            err_messages.append("%s\t%s\tIncorrect file format" % (upid, accession))

        if os.path.getsize(seq_file) == 0:
            err_messages.append("%s\t%s\tEmpty file" % (upid, accession))
    else:
        err_messages.append("%s\t%s\tMissing file" % (upid, accession))

    return err_messages

# ---------------------------------------------------------------------------------------------------------


def validate_genome_download_project(project_dir):
    """
    This function will loop over all directories in a genome download project and run a couple of
    validation steps to evaluate the state of the download and report any issues identified

    project_dir: A valid path to a genome download project directory as generated by the
    genome_downloader pipeline

    return:
    """

    # TODO Split this function into smaller fragments

    upid_gca_file = os.path.join(project_dir, "upid_gca_dict.json")

    # load genome accessions
    fp = open(upid_gca_file)
    upid_gca_dict = json.load(fp)
    fp.close()

    upids = upid_gca_dict.keys()

    # validate all genomes in project directory
    for upid in upids:
        # 1. get updir loc
        upid_idx = upid[-3:]
        subdir_loc = os.path.join(project_dir, upid_idx)
        updir_loc = os.path.join(subdir_loc, upid)

        # check if genome directory does not exist, print message and move to the next one
        if not os.path.exists(updir_loc):
            print "%s\tOutput directory is missing" % upid
            continue

        else:
            # create a genome err file to report any issues
            upid_err_fp = open(os.path.join(updir_loc, upid + '.err'), 'w')
            # check download status was successful
            if check_genome_download_status(os.path.join(updir_loc, "download.out")):
                # check the sequence dir
                seq_dir_loc = os.path.join(updir_loc, "sequences")

                # check if the sequence directory has been created
                if not os.path.exists(seq_dir_loc):
                    print "%s\tSequences directory is missing" % upid
                    upid_err_fp.write("%s\tSequences directory is missing\n" % upid)
                    # close err file and move to next genome
                    upid_err_fp.close()
                    continue

                # if it passes the previous step, check if sequence directory has contents
                seq_items = os.listdir(seq_dir_loc)  # dir or file
                # check if dir is not empty
                if len(seq_items) == 0:
                    print "%s\tEmpty sequence directory" % upid
                    upid_err_fp.write("%s\tEmpty sequence directory\n" % upid)
                    upid_err_fp.close()
                    continue

                # went over all structure issues, work with accessions
                # check if there's a json file with all proteome accessions - generated for all
                uniprot_acc_file = os.path.join(updir_loc, upid + '_accessions.json')

                if not os.path.exists(uniprot_acc_file):
                    print "%s\tMissing proteome_accessions file" % upid
                    upid_err_fp.write("%s\tMissing proteome_accessions file\n" % upid)
                    upid_err_fp.close()
                    continue

                # load genome accessions retrieved from Uniprot to check for any missing sequence files
                proteome_accs = {}
                fp = open(uniprot_acc_file, 'r')
                proteome_accs = json.load(fp)
                fp.close()

                # now look if there is a gca accession
                if upid_gca_dict[upid]["GCA"] != -1:
                    # get GCA report file path
                    gca_report_file = os.path.join(updir_loc,
                                                   upid_gca_dict[upid]["GCA"] + "_sequence_report.txt")

                    # check if gca file exists and load accessions
                    if os.path.exists(gca_report_file):
                        # load accessions from gca report file - without versions
                        fp = open(gca_report_file, 'r')
                        ena_accs = [x.strip().split('\t')[0].partition('.')[0] for x in fp]
                        fp.close()

                        # generate a list of unique genome accessions from Uniprot and ENA
                        genome_unique_accs = []
                        ena_accs.pop(0)  # remove header
                        cleaned_gca_accs = set(ena_accs)
                        upid_other_accs = set(proteome_accs["OTHER"].values())
                        genome_unique_accs = upid_other_accs.union(cleaned_gca_accs)

                        upid_restore_fp = open(os.path.join(updir_loc, "restore.json"), 'w')

                        # default values for restoration procedure
                        accessions_to_restore = {"WGS": proteome_accs["WGS"], "OTHER": [], "MULTI": 0}

                        # check fasta files
                        # case A - sequence directory only contains fasta files
                        if os.path.isfile(os.path.join(seq_dir_loc, seq_items[0])):

                            for accession in genome_unique_accs:
                                seq_file = os.path.join(seq_dir_loc, accession + ".fa")
                                # check if fasta file exists
                                err_messages = check_seq_file(upid, seq_file)

                                if len(err_messages) > 0:
                                    restore = False
                                    for err_message in err_messages:
                                        if err_message.find("format") != -1:
                                            restore = True
                                        upid_err_fp.write(err_message+'\n')

                                    if restore is True:
                                        proteome_accs["OTHER"].append(accession)

                            upid_err_fp.close()
                            json.dump(proteome_accs, upid_restore_fp)

                        # case B - multiple subdirectories
                        else:
                            # to speed up the process we will work with subdirs instead
                            # and remove the accession from the list if everything is ok
                            # with the file. Any missing sequence files will be reported at the end
                            accessions_to_restore["MULTI"] = 1
                            temp_accession_list = list(copy.deepcopy(genome_unique_accs))

                            for subdir in seq_items:
                                subdir_loc = os.path.join(seq_dir_loc, subdir)
                                seq_files = os.listdir(subdir_loc)

                                for seq_file in seq_files:
                                    print "%s\t%s\t%s" % (upid, subdir_loc, seq_file)
                                    sec_file_loc = os.path.join(subdir_loc, seq_file)
                                    err_messages = check_seq_file(upid, seq_file)

                                    if len(err_messages) != 0:
                                        restore = False
                                        for err_message in err_messages:
                                            if err_message.find("format") != -1:
                                                restore = True
                                            upid_err_fp.write(err_message + '\n')

                                        if restore is True:
                                            accessions_to_restore["OTHER"].append(accession)

                                    # delete accession from list
                                    acc = seq_file.partition('.')[0]
                                    temp_accession_list.remove(acc)

                                # Done with all file checks, report any remaining files not found
                                # in the sequence sub directories
                                if len(temp_accession_list) > 0:
                                    for accession in temp_accession_list:
                                        # print "%s\t%s\tMissing file" % (upid, acc)
                                        upid_err_fp.write("%s\t%s\tMissing file\n" % (upid, acc))
                                        accessions_to_restore["OTHER"].append(accession)

                            upid_err_fp.close()
                            json.dump(proteome_accs, upid_restore_fp)

                            if restore is True:
                                print "%s\tRestore download"

                    # no gca report file found, but check if there is a WGS set available
                    # and if the corresponding fasta file has been copied
                    else:
                        if proteome_accs["WGS"] != -1:
                            wgs_prefix = proteome_accs["WGS"][0:6]
                            wgs_file_loc = os.path.join(seq_dir_loc,
                                                        wgs_prefix + ".fasta.gz")

                            if not os.path.exists(wgs_file_loc):
                                print "%s\t%s\tWGS file has not been copied" % (upid,
                                                                                wgs_prefix + ".fasta.gz")
                                upid_err_fp.write("%s\t%s\tWGS file has not been copied\n" % (upid,
                                                                                              wgs_prefix + ".fasta.gz"))
                                upid_err_fp.close()
                                continue

                # check if WGS set has been copied
                else:
                    if proteome_accs["WGS"] != -1:
                        # get wgs prefix and the file location
                        wgs_prefix = proteome_accs["WGS"][0:6]
                        wgs_file_loc = os.path.join(seq_dir_loc,
                                                    wgs_prefix + ".fasta.gz")

                        if not os.path.exists(wgs_file_loc):
                            print "%s\t%s\tWGS file has not been copied" % (upid,
                                                                            wgs_prefix + ".fasta.gz")
                            upid_err_fp.write("%s\t%s\tWGS file has not been copied\n" % (upid,
                                                                                          wgs_prefix + ".fasta.gz"))

                    # check if uniprot accessions match all the accessions in sequence dir
                    else:
                        # check that
                        seq_files = [x.partition('.')[0] for x in os.listdir(seq_dir_loc) if x.endswith('.fa')]
                        uniprot_accs = proteome_accs["OTHER"].values()

                        if len(seq_files) != len(uniprot_accs):
                            for accession in uniprot_accs:
                                if accession not in seq_files:
                                    upid_err_fp.write("%s\t%s\tMissing file\n" % (upid, accession))

                    upid_err_fp.close()

            else:
                # everything ok with lsf, and files but run some luigi checks here??
                # run some luigi checks here
                print "%s\tUnsuccessful download" % upid
                # check if WGS available and use it in restore function or launch download
                # pipeline again


# -----------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    """
    Usage python validate_genomes.py /path/to/project/dir > genome_validation_report.txt
    """
    project_dir = sys.argv[1]
    validate_genome_download_project(project_dir)

