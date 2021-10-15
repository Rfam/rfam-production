import os
import argparse
import subprocess
import time

from scripts.support.mirnas.config import UPDATE_DIR


def check_svn_error(family):
    family_dir = os.path.join(UPDATE_DIR, family)
    err_file = os.path.join(family_dir, "auto_rfci.err")

    with open(err_file) as rqc_output:
        read_output = rqc_output.read()
        if "Filesystem has no item: '/trunk/Families/{0}/SEEDTBLOUT'".format(family) in read_output:
            print('svn add for this family')
            os.chdir(family_dir)
            subprocess.call('svn add SEEDSCORES SEEDTBLOUT', shell=True)
            return True

    return False


def check_successful(rfam_acc):
    family_dir = os.path.join(UPDATE_DIR, rfam_acc)
    err_file = os.path.join(family_dir, "auto_rfci.err")

    with open(err_file) as rqc_output:
        read_output = rqc_output.read()
        if "Successfully checked family in" in read_output:
            return True

    return False


def check_in(rfam_acc, preseed=False):
    family_dir = os.path.join(UPDATE_DIR, rfam_acc)
    out_file = os.path.join(family_dir, "auto_rfci.out")
    message = "\'Update using miRBase seed\'"
    option = ''
    if preseed:
        option = "-preseed"
    if rfam_acc.ignore_seed:
        option = "-i seed"
    os.chdir(UPDATE_DIR)
    subprocess.call("rfci.pl -m {msg} {option} {rfam_acc} > {out_file}".format(
        msg=message, option=option, rfam_acc=rfam_acc, out_file=out_file), shell=True)


def get_families_from_files(qc_passed_file, ignore_seed_file):
    with open(qc_passed_file) as infile:
        file_data = infile.read()
        families = [line for line in file_data.split('\n') if line.strip() != '']

    with open(ignore_seed_file) as infile:
        file_data = infile.read()
        families_ignore_seed = [line for line in file_data.split('\n') if line.strip() != '']

    return families, families_ignore_seed


def check_in_all_families(families_to_ci, families_ignore_seed):
    checked_in = []
    not_checked_in = []
    for family in families_ignore_seed:
        family.ignore_seed = True
    for family in families_to_ci+families_ignore_seed:
        check_in(family)
        time.sleep(60)
        if check_successful(family):
            print('{0} has been checked in'.format(family))
            checked_in.append(family)
        else:
            svn_add = check_svn_error(family)
            if svn_add:
                print('svn add has been applied for {0}. Trying to check in again using -preseed'.format(family))
                check_in(family, preseed=True)
                if check_successful(family):
                    print('{0} has been checked in'.format(family))
                    checked_in.append(family)
            else:
                not_checked_in.append(family)
                print('{0} HAS NOT been checked in. Please check manually {1}'.format(
                    family, os.path.join(UPDATE_DIR, family, "auto_rfci.err")))

    with open(os.path.join(UPDATE_DIR, 'checked_in.txt'), 'w') as outfile:
        for family in checked_in:
            outfile.write(family)
    with open(os.path.join(UPDATE_DIR, 'not_checked_in.txt'), 'w') as outfile:
        for family in not_checked_in:
            outfile.write(family)
    if len(not_checked_in) > 0:
        print('The following families could not be checked in and require manual intervention: {0}'.format(
            not_checked_in))


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc_passed", help="list of families to check in", action="store")
    parser.add_argument("--ignore-seed", help="list of families to check in with -i seed", action="store")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    families, families_iseed = get_families_from_files(args.qc_passed, args.ignore_seed)
    check_in_all_families(families, families_iseed)