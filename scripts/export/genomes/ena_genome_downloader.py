"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import os
import sys
import subprocess

import config.gen_config as gc

# -----------------------------------------------------------------------------

LSF_GROUP = "/rfam_gen/%s"
ENA_TOOL = gc.ENA_TOOLKIT

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


def main(genome_accession_file, project_dir, lsf=True):

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
        upid = genome_data[0]
        domain = genome_data[2]
        gen_acc = genome_data[1]

        # create domain directory
        domain_dir = os.path.join(project_dir, domain)
        if not os.path.exists(domain_dir):
            os.mkdir(domain_dir)

        # create genome directory e.g. project_dir/domain/updir
        updir = os.path.join(project_dir, os.path.join(domain_dir, upid))
        if not os.path.exists(updir):
            os.mkdir(updir)

        if lsf is not True:
            fetch_genome_from_ENA(gen_acc, updir)

        # submit a job
        else:
            err_file = os.path.join(updir, upid+'.err')
            out_file = os.path.join(updir, upid+'.err')
            bsub_cmd = "bsub -o %s -e %s -g %s %s -f fasta -m -d %s %s" % (out_file,
                                                                      err_file,
                                                                      LSF_GROUP % domain,
                                                                      os.path.join(ENA_TOOL, 'enaDataGet'),
                                                                      updir,
                                                                      gen_acc)
            subprocess.call(bsub_cmd, shell=True)


# -----------------------------------------------------------------------------

if __name__ == "__main__":

    genome_accession_file = sys.argv[1]
    project_dir = sys.argv[2]

    if 'lsf' in sys.argv:
        main(genome_accession_file, project_dir, lsf=True)

    else:
        main(genome_accession_file, project_dir, lsf=False)
