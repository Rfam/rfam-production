import os
import sys
import subprocess

import config.gen_config as gc

# -----------------------------------------------------------------------------


def fetch_genome_from_ENA(genome_accession, dest_dir):
    """
    Uses ENAs enaBrowserTools to download a specific directory

    genome_accession: A valid ENA accession to Download
    dest_dir: Destination directory where genome will be downloaded

    return: void
    """

    exec_path = os.path.join(gc.ENA_TOOLKIT, 'enaDataGet')
    cmd = "%s -f fasta -m -d %s %s" % (exec_path, dest_dir, genome_accession)

    subprocess.call(cmd, shell=True)

# -----------------------------------------------------------------------------


def main(genome_accession_file, project_dir):

    """
    Parses a file of genome accessions and downloads genomes from ENA. Genomes
    are downloaded in fasta format. It is a requirement that genome_accession
    file containes a GCA or WGS accession per genome

    genome_accession_file: A file with a list of upid\tGCA\tdomain or
    upid\tWGS\tdomain pairs
    project_dir: The path to the directory where all genomes will be downloaded
    It will be created if it does not exist

    return: void
    """

    # generating project_directory
    if not os.path.exists(project_dir):
        os.mkdir(project_dir)

    # create a file handle and read genome id pairs
    input_fp = open(genome_accession_file, 'r')

    for genome in input_fp:
        genome_data = genome.strip().split('\t')

        # create domain directory
        domain_dir = os.path.join(project_dir, genome_data[2])
        if not os.path.exists(domain_dir):
            os.mkdir(domain_dir)

        # create genome directory e.g. project_dir/domain/updir
        updir = os.path.join(project_dir, os.path.join(domain_dir, genome_data[0]))
        if not os.path.exists(updir):
            os.mkdir(updir)

        fetch_genome_from_ENA(genome_data[1], updir)


# -----------------------------------------------------------------------------

if __name__ == "__main__":

    genome_accession_file = sys.argv[1]
    project_dir = sys.argv[2]

    main(genome_accession_file, project_dir)