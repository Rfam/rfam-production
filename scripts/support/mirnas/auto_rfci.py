import os
import argparse
import subprocess
import time

from scripts.support.mirnas.update_mirnas_helpers import UPDATE_DIR, MEMORY, CPU, LSF_GROUP


def check_error_file(family):
    family_dir = os.path.join(UPDATE_DIR, family)
    err_file = os.path.join(family_dir, "auto_rfci.err")
    with open(err_file) as f:
        f_read = f.read()
        if "Filesystem has no item: '/trunk/Families/{0}/SEEDTBLOUT'".format(family) in f_read:
            print('svn add for this family')
            os.chdir(family_dir)
            subprocess.call('svn add SEEDSCORES SEEDTBLOUT', shell=True)
            return True

    return False


def check_successful(rfam_acc):
    family_dir = os.path.join(UPDATE_DIR, rfam_acc)
    lsf_err_file = os.path.join(family_dir, "auto_rfci.err")

    with open(lsf_err_file) as rqc_output:
        read_output = rqc_output.read()
        if "Successfully checked family in" in read_output:
            return True

    return False


def check_in(rfam_acc, preseed=False):
    family_dir = os.path.join(UPDATE_DIR, rfam_acc)
    lsf_err_file = os.path.join(family_dir, "auto_rfci.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfci.out")
    job_name = rfam_acc
    if preseed:
        option = 'Update using miRBase seed -preseed'
    else:
        option = 'Update using miRBase seed'
    cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} -J {job_name} "
           "\"cd {update_dir} && rfci.pl -m {option} {rfam_acc}\"")
    subprocess.call(
        cmd.format(
            mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP,
            job_name=job_name, update_dir=UPDATE_DIR, option=option, rfam_acc=rfam_acc), shell=True)


def check_in_all_families(families):
    checked_in = []
    not_checked_in = []
    with open(families) as infile:
        file_data = infile.read()
        all_families = [line for line in file_data.split('\n') if line.strip() != '']
        for family in all_families:
            check_in(family)
            time.sleep(120)
            if check_successful(family):
                print('{0} has been checked in'.format(family))
                checked_in.append(family)
            else:
                svn_add = check_error_file(family)
                if svn_add:
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
    parser.add_argument("--input", help="list of families to check in/ passed QC", action="store")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    check_in_all_families(args.input)
