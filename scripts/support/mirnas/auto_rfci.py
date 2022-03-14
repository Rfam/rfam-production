import os
import argparse
import subprocess
import time

from scripts.support.mirnas.mirna_config import UPDATE_DIR
from scripts.support.mirnas.update_mirnas_helpers import MEMORY, CPU, LSF_GROUP

checked_in = []
not_checked_in = []


def check_svn_error(family):
    """
    Check the output to decide if `svn add` needs to be applied
    :param family: accession number of family
    :return: True if `svn add` applied, else False
    """

    family_dir = os.path.join(UPDATE_DIR, family)
    out_file = os.path.join(family_dir, "auto_rfci.err")

    with open(out_file) as rqc_output:
        read_output = rqc_output.read()
        if "Filesystem has no item: '/trunk/Families/{0}/SEEDTBLOUT'".format(family) in read_output:
            print('svn add for this family')
            os.chdir(family_dir)
            subprocess.call('svn add SEEDSCORES SEEDTBLOUT', shell=True)
            return True

    return False


def check_successful(rfam_acc):
    """
    Check if the output contains a success message
    :param rfam_acc: accession number of family
    :return: True if family has been checked in, else False
    """

    family_dir = os.path.join(UPDATE_DIR, rfam_acc)
    out_file = os.path.join(family_dir, "auto_rfci.err")

    with open(out_file) as rqc_output:
        read_output = rqc_output.read()
        if "Successfully checked family in" in read_output:
            return True

    return False


def check_in(acc, pre_seed=False, ignore_seed=False):
    """
    Call rfci.pl to check in the family
    :param acc: accession number of family
    :param pre_seed: True if using -preseed option, else False
    :param ignore_seed: True if using -i seed option, else False
    """
    family_dir = os.path.join(UPDATE_DIR, acc)
    lsf_err_file = os.path.join(family_dir, "auto_rfci.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfci.out")
    cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} -J {job_name} "
           "\"cd {update_dir} && yes | rfci.pl -m {msg} {option} {rfam_acc}\"")
    message = "\'Update using miRBase seed\'"
    option = ''
    if pre_seed:
        option += "-preseed "
    if ignore_seed:
        option += "-i seed "
    cmd = cmd.format(mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP,
                     job_name=acc,
                     update_dir=UPDATE_DIR, msg=message, option=option, rfam_acc=acc)
    subprocess.call(cmd, shell=True)


def call_check_in(families, ignore_seed):
    """
    Call check in and determine if family has been successfully checked in
    :param families: accession numbers of the families to check in
    :param ignore_seed: True if using -i seed option, else False
    :return:
    """
    for family in families:
        check_in(family, ignore_seed=ignore_seed)
        time.sleep(20)
        if check_successful(family):
            print('{0} has been checked in'.format(family))
            checked_in.append(family)
        else:
            svn_add = check_svn_error(family)
            if svn_add:
                print('svn add has been applied for {0}. Trying to check in again using -preseed'.format(family))
                check_in(family, pre_seed=True, ignore_seed=ignore_seed)
                if check_successful(family):
                    print('{0} has been checked in'.format(family))
                    checked_in.append(family)
            else:
                not_checked_in.append(family)
                print('{0} HAS NOT been checked in. Please check manually {1}'.format(
                    family, os.path.join(UPDATE_DIR, family, "auto_rfci.err")))

    return checked_in, not_checked_in


def write_result_to_files(checked_in, not_checked_in):
    """
    Write to files, the list of families checked in and not checked in
    """

    with open(os.path.join(UPDATE_DIR, 'checked_in.txt'), 'a') as outfile:
        for family in checked_in:
            outfile.write(family)
    with open(os.path.join(UPDATE_DIR, 'not_checked_in.txt'), 'a') as outfile:
        for family in not_checked_in:
            outfile.write(family)
    if len(not_checked_in) > 0:
        print('The following families could not be checked in and require manual intervention: {0}'.format(
            not_checked_in))


def check_in_all_families(families, ignore_seed):
    """
    Start the process of checking in all given families
    :param families: list of accession numbers for families to check in
    :param ignore_seed: True if using -i seed option, else False
    """
    checked, not_checked = call_check_in(families, ignore_seed)
    write_result_to_files(checked, not_checked)


def get_families_from_files(qc_passed_file):
    """
    Get the list of families to check in from the given text file
    :param qc_passed_file: text file
    :return: list of families to check in
    """
    with open(qc_passed_file) as infile:
        file_data = infile.read()
        families = [line for line in file_data.split('\n') if line.strip() != '']

    return families


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Text file of the families that have passed QC checks, and can be checked in",
                        action="store")
    parser.add_argument("--ignore-seed", help="True if families should be checked in with -i seed",
                        action="store_true", default=False)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    families = get_families_from_files(args.input)
    check_in_all_families(families, ignore_seed=args.ignore_seed)
