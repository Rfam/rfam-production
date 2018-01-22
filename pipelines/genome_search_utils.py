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
import json
import shutil
import subprocess
import luigi

from config import rfam_local as conf
from config import gen_config as gc
from support import merge_fasta as mf

from utils import genome_search_utils as gsu

# add parent directory to path
if __name__ == '__main__' and __package__ is None:
    os.sys.path.append(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ----------------------------------TASKS--------------------------------------


class SplitGenomeFasta(luigi.Task):
    """
    Luigi task to split genome files in smaller chunks for efficient searching
    """
    updir = luigi.Parameter()
    upid = luigi.Parameter()

    def run(self):
        """
        Main function that organises genome search directories based on
        genome size
        """
        # get updir location
        upid_fasta = os.path.join(self.updir, self.upid + '.fa')
        seq_chunks_dir = os.path.join(self.updir, "search_chunks")

        if not os.path.exists(seq_chunks_dir):
            os.mkdir(seq_chunks_dir)
            os.chmod(seq_chunks_dir, 0777)

            # check if we need to split the seq_file
            if gsu.count_nucleotides_in_fasta(upid_fasta) >= gs.SPLIT_SIZE:
                # split sequence file into smalled chunks
                gsu.split_seq_file(upid_fasta, gc.SPLIT_SIZE, dest_dir=seq_chunks_dir)

            # for input consistency if the sequence file is small, copy it in the
            # search_chunks directory
            else:
                # copy file
                shutil.copyfile(upid_fasta, os.path.join(seq_chunks_dir,
                                                         self.upid + '.fa'))
                # index file
                cmd = "%s --index %s" % (conf.ESL_SFETCH, os.path.join(seq_chunks_dir,
                                                                       self.upid + '.fa'))
                subprocess.call(cmd, shell=True)

# -----------------------------------------------------------------------------


class MergeGenomeFasta(luigi.Task):
    """
    Task to merge genome fasta files
    """
    updir = luigi.Parameter()

    def run(self):
        """
        Merge all fasta files in updir
        """
        mf.merge_genome_files(self.updir)

    def output(self):
        """
        Check genome fasta file has been generated
        """
        upid = os.path.split(self.updir)[1]
        genome_fasta = os.path.join(self.updir, upid + '.fa')

        return luigi.LocalTarget(genome_fasta)


# -----------------------------------------------------------------------------


class GenomeSearchUtilsEngine(luigi.Task):
    """
    A pipeline providing various utilities for manipulating Genome search data
    """
    project_dir = luigi.Parameter(description="Genome download project directory")

    genome_list = luigi.Parameter(default=None,
                                  description="A list of upids to process")

    utility = luigi.Parameter(default=None,
                             description="Utility to use to a genome (e.g famerge)")

    lsf = luigi.BoolParameter(default=True,
                              description="If used then run on lsf, otherwise run locally")

    def run(self):
        """
        Call load_upid_gca_file to export all upid_gca accession pairs in
        json format.
        """
        id_pairs = None

        upids = []
        if self.genome_list is None:
            upid_gca_file_loc = os.path.join(self.project_dir, "upid_gca_dict.json")
            upid_gca_fp = open(upid_gca_file_loc, 'r')
            accessions = json.load(upid_gca_fp)
            upid_gca_fp.close()

            upids = accessions.keys()

        else:
            upid_list_fp = open(self.genome_list, 'r')
            upids = [x.strip() for x in upid_list_fp]
            upid_list_fp.close()

        cmd = ''
        for upid in upids:
            subdir = os.path.join(self.project_dir, upid[-3:])
            updir = os.path.join(subdir, upid)

            if os.path.exists(updir):
                # Merge Genomes
                if self.method.lower() == 'famerge':
                    if self.lsf is True:
                        cmd = "bsub -M %s -R \"rusage[mem=%s,tmp=%s]\" -o %s -e %s -u %s -Ep \"rm -rf luigi\" " \
                              "-g %s python %s MergeGenomeFasta --updir %s" % (gc.MEM, gc.MEM, gc.TMP_MEM,
                                                                          os.path.join(updir, "merge.out"),
                                                                          os.path.join(updir, "merge.err"),
                                                                          gc.USER_EMAIL, gc.LSF_GEN_GROUP,
                                                                          os.path.realpath(__file__), updir)
                    else:
                        cmd = "python \"{this_file}\" MergeGenomeFasta --updir {upid}".format(
                            this_file=os.path.realpath(__file__),
                            updir=updir)

                elif self.method.lower() == 'fasplit':
                    if self.lsf is True:
                        cmd = "bsub -M %s -R \"rusage[mem=%s,tmp=%s]\" -o %s -e %s -u %s -Ep \"rm -rf luigi\" " \
                              "-g %s python %s SplitGenomeFasta --updir %s --upid %s" % (gc.MEM, gc.MEM, gc.TMP_MEM,
                                                                      os.path.join(updir, "merge.out"),
                                                                      os.path.join(updir, "merge.err"),
                                                                      gc.USER_EMAIL, gc.LSF_GEN_GROUP,
                                                                      os.path.realpath(__file__), updir, upid)

                    else:
                        cmd = "python \"{this_file}\" SplitGenomeFasta --updir {upid} --upid {upid}".format(
                            this_file=os.path.realpath(__file__),
                            updir=updir,
                            upid=upid)

            subprocess.call(cmd, shell=True)

            cmd = ''

# -----------------------------------------------------------------------------


if __name__ == '__main__':
    luigi.run()
