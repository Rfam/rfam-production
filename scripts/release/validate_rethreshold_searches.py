import argparse
import os


def load_rfam_accessions_from_file(file_of_accessions):
    """
    Parses a .txt file containing Rfam accessions

    :param file_of_accessions: .txt file containing a list of Rfam accessions
    :return: list of Rfam family accessions
    """
    with open(file_of_accessions, 'r') as fp:
        accs = [x.strip() for x in fp]

    return accs


def validate(file_name, all_families_dir):
    """
    Check the lsf output file and log files to validate if the process completed successfully

    :param file_name: .txt with list of Rfam accessions
    :param all_families_dir: destination directory to checkout the families
    """
    validation_file = os.path.join(all_families_dir, "validation.log")
    accessions = load_rfam_accessions_from_file(file_name)
    with open(validation_file, 'w') as vf:
        for acc in accessions:
            if not is_valid(all_families_dir, acc):
                vf.write(acc + '\n')

    if os.path.getsize(validation_file) == 0:
        print ("Validation process completed. All searches completed successfully.")

    else:
        print ("Validation process found errors! Check validation.log for families with errors")


def is_valid(all_families_dir, rfam_acc):
    """
    Checks if the job ran successfully by checking if .err file is empty and
    that 'Success' keyword exists in .out file. As an additional sanity check, we
    look for the rfsearch.log file as an indication that rfsearch actually ran.

    :return: True if the family passes validation, False otherwise
    """

    family_dir = os.path.join(all_families_dir, rfam_acc)

    if not os.path.exists(os.path.join(family_dir, "rfsearch.log")):
        return False

    if not os.path.getsize(os.path.join(family_dir, "auto_rfsearch.err")) == 0:
        return check_rfsearch_log_success(family_dir)

    lsf_out_file = os.path.join(family_dir, "auto_rfsearch.out")
    with open(lsf_out_file, 'r') as outfile:
        contents = outfile.read()
        if "Successfully completed." not in contents:
            return False
    return True


def check_rfsearch_log_success(family_dir):
    """
    Checks if the rfsearch.log file contains the success string # [ok] in
    order to mark the family as successfully completed.

    :param family_dir: Rfam family directory
    :return: True if the family passes validation, False otherwise
    """

    rfsearch_log_file = os.path.join(family_dir, "rfsearch.log")
    with open(rfsearch_log_file, 'r') as log:
        contents = log.read()
        if "# [ok]" not in contents:
            return False

    return True


def parse_arguments():
    """
    Parse the command line arguments

    :return: Argparse parser object
    """
    parser = argparse.ArgumentParser(description='Re-threshold the families for new rfamseq')

    req_args = parser.add_argument_group("required arguments")
    req_args.add_argument('--dir', help='directory where families have been checked out',
                          type=str, required=True)
    req_args.add_argument('-f', help='a file containing a list of Rfam family accessions', type=str)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    validate(args.f, args.dir)
