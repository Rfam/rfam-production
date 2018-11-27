import os
import sys
import subprocess

# -----------------------------------------------------------------

MEMORY = 6000

CMD = ("bsub -M %s cd %s && "
       "rfsearch.pl -t 30 -cnompi && rfmake.pl -t %s "
        "-a -forcethr && mkdir rscape-seed && R-scape "
        "--outdir rscape-seed --cyk align && mkdir rscape-align && "
        "R-scape --outdir rscape-align --cyk align && "
        "cd .. && rqc-all.pl %s")


# -----------------------------------------------------------------

def get_threshold_from_DESC(desc_file):

    """
    Parses a family DESC file and returns the gathering
    threshold

    desc_file: The path to a valid Rfam family DESC file

    return: A string which is the corresponding gathering threshold
    """

    threshold = 0
    desc_fp = open(desc_file, 'r')

    for line in desc_fp:
        # look for the GA line
        if line[0:2] == 'GA':
            # extract threshold which after splitting
            # is the last element in the list
            threshold = line.strip().split(' ')[-1]

    return threshold

# -----------------------------------------------------------------

def main(family_dir, multi=False):
    """
    Launches LSF jobs to re-threshold Rfam families

    family_dir: A directory to a single Rfam family or multiple
    Rfam family subdirectories

    multi: A boolean variable specifying a single on multi-launch

    return: None
    """

    if multi is False:
        desc_loc = os.path.join(family_dir, "DESC")
        threshold = get_threshold_from_DESC(desc_loc)

        cmd = CMD % (MEMORY, family_dir, threshold, os.path.basename(family_dir))

        subprocess.call(cmd, shell=True)

    else:
        family_sub_dirs = os.listdir(family_dir)

        for subdir in family_sub_dirs:
            desc_loc = os.path.join(family_dir, os.path.join(subdir, "DESC"))
            threshold = get_threshold_from_DESC(desc_loc)

            cmd = CMD % (MEMORY, subdir, threshold, subdir)

            # this needs to be replaced with a command to submit a job to the
            # cluster
            subprocess.call(cmd, shell=True)


# -----------------------------------------------------------------

if __name__ == '__main__':

    family_dir = sys.argv[1]

    if len(sys.argv) == 3:
        if '-m' in sys.argv:
            main(family_dir, multi=True)
        else:
            print "\nSecond parameter is wrong!\n"
    else:
        main(family_dir, multi=False)
