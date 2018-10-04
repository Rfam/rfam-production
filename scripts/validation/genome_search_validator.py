"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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
Genome search validation modules
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import sys

# -----------------------------------------------------------------------------

def check_genome_search_success(error_file):
    """
    Checks whether genome search was successful by checking lsf error file
    size. If the file is empty return success, otherwise return 0

    error_file (string): A string representing the path to LSF's job error
    file (-e)
    """
    success = 1

    if os.path.getsize(error_file) == 0:
        return success

    return 0

# -----------------------------------------------------------------------------

def get_search_recovery_list(input_dir, recovery_file=True):
    """
    Loops over all genome search lsf error files and outputs and list all
    genomes that crashed, so that we can restore those searches only.
    If destination directory is None then the function will return a list of
    Uniprot's UPID accessions

    input_dir (string): A string representing the path to a domain directory where
    all genome search output files are located (e.g. eukaryota/UPXXXXXXXXX).
    recovery_file (boolean): If True the function will generate a recovery file in
    input directory, otherwise it will return a list of accessions
    """

    recovery_list = []

    # list directory contents
    if os.path.exists(input_dir):
        genome_dirs = [x for x in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, x))]

        # loop over all genome/proteome directories
        for genome in genome_dirs:
            genome_dir = os.path.join(input_dir, genome)

            # get lsf job status
            exec_status = check_genome_search_success(os.path.join(genome_dir, genome+".err"))

            # if the job was unsuccessful list the upid for recovery
            if exec_status == 0:
                recovery_list.append(genome)

        if recovery_file is True:
            fp_out = open(os.path.join(input_dir, "recovery_list.txt"), 'w')
            for upid in recovery_list:
                fp_out.write(upid+'\n')

            fp_out.close()

    else:
        sys.exit("Error message here.")

    # a list of upids to recover
    return recovery_list

# -----------------------------------------------------------------------------

def check_search_err_files(search_output_dir):
    """
    Lookup all output subdirectories and check for cases that .err files are
    not empty

    search_output_dir: search output directory as organised by genome_search

    returns: A dictionary with all erroneous cases
    """

    search_err_cases = {}

    output_subdirs = os.listdir(search_output_dir)

    for subdir in output_subdirs:
        subdir_loc = os.path.join(search_output_dir, subdir)
        updirs = os.listdir(subdir_loc)
        for updir in updirs:
            updir_loc = os.path.join(subdir_loc, updir)
            err_files = [x for x in os.listdir(updir_loc) if x.endswith(".err")]
            gen_err_cases = []
            for err_file in err_files:
                if os.path.getsize(os.path.join(updir_loc,err_file))>0:
                    gen_err_cases.append(err_file[0:-4])

            if len(gen_err_cases) > 0:
                if subdir not in search_err_cases.keys():
                    search_err_cases[subdir] = {updir: gen_err_cases}
                else:
                    search_err_cases[subdir][updir] = gen_err_cases

    return search_err_cases

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    pass
