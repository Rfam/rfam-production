import argparse
import os
import shutil
import subprocess
import sys

# ---------------------------- Variable initialization ---------------------------

MEMORY = 2000
CPU = 4
LSF_GROUP = "/family_srch"
REQUIRED_FILES = ["SEED", "DESC", "species", "outlist", "seedoutlist"]


# --------------------------------------------------------------------------------


def check_seed_format(seed_alignment, format="stockholm"):
    """
    Checks if a SEED alignment is in a specific format
    return: void
    """
    pass


# --------------------------------------------------------------------------------


def create_family_directory(dirname, seed_alignment, dest_dir):
    """
    Creates a new Rfam family directory and sets all the files required
    to launch rfsearch.

    dirname: A name for the new Rfam family in the format str_str
    seed_alignment: A valid seed alignment in stockholm format
    dest_dir: The base directory where the new family directory will be

    return: The path to the family directory, None otherwise
    """

    # check if destination directory existis
    if not os.path.exists(dest_dir):
        sys.exit(
            "\nDestination directory does not exist!\n"
        )  # go on with creating a new family directory

    family_dir = os.path.join(dest_dir, dirname)
    if not os.path.exists(family_dir):
        os.mkdir(family_dir)

    # check if the family directory was created successfully
    if os.path.exists(family_dir):

        # copy SEED to family dir
        shutil.copyfile(seed_alignment, os.path.join(family_dir, "SEED"))
        return family_dir

    return None


# --------------------------------------------------------------------------------


def launch_new_rfsearch(family_dir, cpu=4):
    """
    Launches a new LSF job

    family_dir: The location of a valid family directory
    cpus: Number of CPUs to use per thread

    return: void
    """

    lsf_err_file = os.path.join(family_dir, "auto_rfsearch.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfsearch.out")

    job_name = os.path.basename(family_dir)

    # LSF command to be executed
    cmd = (
        "bsub -M %s -o %s -e %s -n %s -g %s -q production-rh74 "
        '-J %s "cd %s && rfsearch.pl -t 30 -q production-rh74 -relax -nodesc"'
    )

    # call command
    subprocess.call(
        cmd
        % (MEMORY, lsf_out_file, lsf_err_file, cpu, LSF_GROUP, job_name, family_dir),
        shell=True,
    )


# --------------------------------------------------------------------------------


def parse_arguments():
    """
    Does basic command line argument parsing

    return: Argparse object
    """

    parser = argparse.ArgumentParser(description="Rfam family Auro-Builder")
    # group required arguments together

    req_args = parser.add_argument_group("required arguments")
    req_args.add_argument(
        "--input",
        help="a directory with multiple SEEDs or a single SEED",
        type=str,
        required=True,
    )
    req_args.add_argument(
        "--dest-dir",
        help="destination directory where create new \
                    family directories",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--cpu", help="number of CPUs to use per thread", type=int, action="store"
    )

    return parser


# --------------------------------------------------------------------------------


def check_required_files_exist(family_dir):
    """
    Checks if all required family files were produced.

    return: True if all files exist, False otherwise
    """

    rfam_files = dict.fromkeys(os.listdir(family_dir))

    for file_type in REQUIRED_FILES:
        if file_type not in rfam_files:
            return False

    # if it passes the check then all required files exists
    return True


# --------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = parse_arguments()
    args = parser.parse_args()

    source_path = ""
    if os.path.isdir(args.input):
        seeds = [x for x in os.listdir(args.input) if x[0] != "."]
        source_path = args.input
    else:
        seeds = [args.input]
        source_path = os.path.split(args.input)[0]

    for seed in seeds:
        dirname = seed.partition(".")[0]
        seed_path = os.path.join(source_path, seed)

        family_dir = os.path.join(args.dest_dir, dirname)
        if not os.path.exists(family_dir):
            family_dir = create_family_directory(dirname, seed_path, args.dest_dir)

        if family_dir is not None:
            launch_new_rfsearch(family_dir, args.cpu)
