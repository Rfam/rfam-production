import argparse
import os
import shutil
import subprocess
import sys

MEMORY = 2000
COPIES_DIR = "/nfs/production/agb/rfam/RELEASES/15.0/test_rfamseq/copies/"


def run_searches(file_name, dest_dir, rfmake):
    """
    Call processes to checkout the family and start the rfsearch (and rfmake, if required) jobs

    :param file_name: .txt file with list of Rfam accessions
    :param dest_dir: destination directory to checkout the families
    :param rfmake: default False, if True run rmake on completion of rfsearch
    """
    if not os.path.isfile(file_name):
        print ("The file location you provided does not exist!\n")
        sys.exit()

    accessions = load_rfam_accessions_from_file(file_name)
    for rfam_acc in accessions:
        checkout_and_search_family(rfam_acc, dest_dir, rfmake)


def load_rfam_accessions_from_file(file_of_accessions):
    """
    Parses a .txt file containing Rfam accessions

    :param file_of_accessions: .txt file containing a list of Rfam accessions
    :return: list of Rfam family accessions
    """
    with open(file_of_accessions, 'r') as fp:
        accs = [x.strip() for x in fp]

    return accs


def checkout_and_search_family(rfam_acc, dest_dir, rfmake=False):
    """
    This function combines family checkout (rfco.pl) and re-scoring of hits
    using rfsearch.pl. If the family directory already exists, then the
    checkout step will be ignored

    :param rfam_acc: A valid Rfam family accession (RFXXXXX)
    :param dest_dir: A valid destination directory, where to checkout the family
    :param rfmake: If True, run rfmake after rfsearch completes. Default False
    """

    family_dir = os.path.join(dest_dir, rfam_acc)

    if not os.path.exists(family_dir):
        os.chdir(dest_dir)
        checkout_family(rfam_acc)
        # copy checked out family dir to copies dir for comparison
        if not os.path.exists(COPIES_DIR + "/" + rfam_acc):
            shutil.copytree(family_dir, COPIES_DIR)

    submit_new_rfsearch_job(family_dir, rfmake)


def checkout_family(rfam_acc):
    """
    Checks out a family from Rfam based on a valid Rfam accession.

    :param rfam_acc: A valid Rfam accession
    """
    cmd = "rfco.pl %s" % rfam_acc
    subprocess.call(cmd, shell=True)


def submit_new_rfsearch_job(family_dir, rfmake=False):
    """
    Submits a new lsf job that runs rfsearch to update SCORES.
    If no threshold is set with rfsearch.pl, it uses existing thresholds by default.
    Submit rfmake job if rfmake is True.

    :param family_dir: The physical location of the family directory
    :param rfmake: If True, run rfmake after rfsearch completes. Default False
    """

    lsf_err_file = os.path.join(family_dir, "auto_rfsearch.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfsearch.out")
    cmd = ("bsub -M %s -R \"rusage[mem=%s]\" -o %s -e %s -q bigmem "
           "\"cd %s && rfsearch.pl --scpu 0 -cnompi -q short -relax\"")

    if rfmake is True:
        cmd = ("bsub -M %s -R \"rusage[mem=%s]\" -o %s -e %s -q bigmem "
               "\"cd %s && rfsearch.pl --scpu 0 -cnompi -q short -relax && rfmake.pl -local\"")

    subprocess.call(cmd % (MEMORY, MEMORY, lsf_out_file, lsf_err_file, family_dir), shell=True)


def parse_arguments():
    """
    Parse the command line arguments

    :return: Argparse parser object
    """
    parser = argparse.ArgumentParser(description='Re-threshold the families for new rfamseq')

    req_args = parser.add_argument_group("required arguments")
    req_args.add_argument('--dest-dir', help='destination directory where to checkout families',
                          type=str, required=True)
    req_args.add_argument('-f', help='a file containing a list of Rfam family accessions', type=str)
    parser.add_argument('--rfmake', help='run rfmake after rfsearch completion', action="store_true")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    run_searches(args.f, args.dest_dir, args.rfmake)
