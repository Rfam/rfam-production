import os
import sys

from config import gen_config as gc
from scripts.export.genomes import genome_fetch as gflib

# ------------------------------------------------------------------------

def setup_genome_directory(upid_dir):
    """
    Setup proteome directory.
    """

    # finally, create the directory if it does not exist
    if not os.path.exists(upid_dir):
        os.mkdir(upid_dir)
        os.chmod(upid_dir, 0777)

    # create a sequence directory where all fasta files will be downloaded
    sequence_dir = os.path.join(upid_dir, "sequences")

    if not os.path.exists(sequence_dir):
        os.mkdir(sequence_dir)
        os.chmod(sequence_dir, 0777)

# ------------------------------------------------------------------------


def download_genome(project_dir, upid):
    """

    project_dir:
    upid:

    return:
    """
    missing_accessions = []

    sub_dir_index = upid[8:]

    # Construct the path for the new proteome directory
    upid_dir = os.path.join(
        os.path.join(project_dir, sub_dir_index), upid)

    # Step 1 - Setup the genome directory structure
    setup_genome_directory(upid_dir)
    sequence_dir = os.path.join(upid_dir, 'sequences')

    max_combinations = 999

    # Fetch proteome accessions, this will also copy GCA file if available
    genome_accessions = gflib.get_genome_unique_accessions(upid, to_file=True,
                                                           output_dir=upid_dir)
    # fetch all other accessions other than WGS and GCA
    wgs_set = None
    other_accessions = genome_accessions["OTHER"]

    # 1. check for assembly report file
    if genome_accessions["GCA"] != -1:
        # fetch wgs set from ENA
        if len(other_accessions) == 0 and genome_accessions["WGS"] == -1:
            wgs_set = gflib.extract_wgs_acc_from_gca_xml(genome_accessions["GCA"])

        if wgs_set is not None or genome_accessions["GCA_NA"] == 1:
            if genome_accessions["GCA_NA"] == 1:
                wgs_set = genome_accessions["WGS"]

            gflib.copy_wgs_set_from_ftp(wgs_set, sequence_dir)

    elif genome_accessions["WGS"] != -1 and genome_accessions["GCA"] == -1:
        # First copy WGS set in upid dir
        gflib.copy_wgs_set_from_ftp(genome_accessions["WGS"], sequence_dir)

    # this should be done in all cases
    # download genome accessions in proteome directory
    if len(other_accessions) > 0:
        if len(other_accessions) < gc.MAX_ALLOWED_FILES:
            for acc in other_accessions:
                status = gflib.fetch_ena_file(acc, "fasta", sequence_dir, compressed=False)
                # add accession to list if file not downloaded
                if status is False:
                    missing_accessions.append(acc)

        # split fasta files in multiple directories
        else:
            # get ranges for
            subdir_ranges = gflib.get_genome_subdirectory_ranges(other_accessions)

            # generate subdirs
            for subdir_index in subdir_ranges:
                if not os.path.exists(os.path.join(sequence_dir, str(subdir_index))):
                    os.mkdir(os.path.join(sequence_dir, str(subdir_index)))

            for accession in other_accessions:
                idx = 0
                acc_index = accession[-3:]
                i = 0

                # find directory index
                while i < len(subdir_ranges) and subdir_ranges[i] < acc_index:
                    i += 1

                if i < len(subdir_ranges):
                    idx = subdir_ranges[i]
                else:
                    idx = max_combinations

                subdir = os.path.join(sequence_dir, str(idx))

                status = gflib.fetch_ena_file(accession, "fasta", subdir, compressed=False)

                if status is False:
                    missing_accessions.append(accession)

    # report any missing accessions
    if len(missing_accessions) > 0:
        fp = open(os.path.join(upid_dir, "unavailable_accessions.txt"), 'w')

        for accession in missing_accessions:
            fp.write(accession+'\n')

        fp.close()

# ------------------------------------------------------------------------

if __name__ == "__main__":

    project_dir = sys.argv[1]
    upid = sys.argv[2]

    download_genome(project_dir, upid)