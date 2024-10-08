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
from scripts.support import merge_fasta as mf

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
            os.chmod(seq_chunks_dir, 0o777)

            # check if we need to split the seq_file
            if gsu.count_nucleotides_in_fasta(upid_fasta) >= gc.SPLIT_SIZE:
                # split sequence file into smalled chunks
                gsu.split_seq_file(upid_fasta, gc.SPLIT_SIZE, dest_dir=seq_chunks_dir)

                # now index the fasta files
                seq_files = os.listdir(seq_chunks_dir)
                for seq_file in seq_files:
                    seq_file_loc = os.path.join(seq_chunks_dir, seq_file)
                    cmd = "%s --index %s" % (conf.ESL_SFETCH, seq_file_loc)
                    subprocess.call(cmd, shell=True)

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


class MergeGenomeTBLOUT(luigi.Task):
    """
    Task to merge genome fasta files
    """
    updir = luigi.Parameter()
    upid = luigi.Parameter()

    def run(self):
        """
        Merge all tbl files in updir
        """

        upid = os.path.basename(self.updir)
        results_dir = os.path.join(self.updir, "search_output")
        res_files = [x for x in os.listdir(results_dir) if x.endswith('.tbl')]

        up_tbl = open(os.path.join(self.updir, self.upid +'.tbl'), 'w')

        for res_file in res_files:
            fp_in = open(os.path.join(results_dir, res_file), 'r')
            for line in fp_in:
                up_tbl.write(line)

            fp_in.close()

        up_tbl.close()

    def output(self):
        """
        Check genome tbl file has been generated
        """
        genome_tbl = os.path.join(self.updir, self.upid + '.tbl')

        return luigi.LocalTarget(genome_tbl)

# -----------------------------------------------------------------------------


class RewriteCleanFasta(luigi.Task):
    """
    Task to merge genome fasta files
    """
    updir = luigi.Parameter()
    upid = luigi.Parameter()

    def run(self):
        """
        Merge all fasta files in updir
        """
        up_fasta = os.path.join(self.updir, self.upid)
        gsu.cleanup_illegal_lines_from_fasta(up_fasta, dest_dir=self.updir)

    def output(self):
        """
        Check that a clean fasta file has been generated
        """
        upid = os.path.split(self.updir)[1]
        clean_fasta = os.path.join(self.updir, upid + '_cleaned.fa')

        return luigi.LocalTarget(clean_fasta)

# -----------------------------------------------------------------------------


class GenomeSearchUtilsEngine(luigi.Task):
    """
    A pipeline providing various utilities for manipulating Genome search data
    """
    project_dir = luigi.Parameter(description="Genome download project directory")

    genome_list = luigi.Parameter(default=None,
                                  description="A list of upids to process")

    tool = luigi.Parameter(default=None,
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
                if self.tool.lower() == 'famerge':
                    if self.lsf is True:
                        cmd = "bsub -M %s -R \"rusage[mem=%s,tmp=%s]\" -o %s -e %s -u %s -Ep \"rm -rf luigi\" " \
                              "-g %s python %s MergeGenomeFasta --updir %s" % (gc.MEM, gc.MEM, gc.TMP_MEM,
                                                                          os.path.join(updir, "merge.out"),
                                                                          os.path.join(updir, "merge.err"),
                                                                          gc.USER_EMAIL, gc.SRCH_GROUP,
                                                                          os.path.realpath(__file__), updir)
                    else:
                        cmd = "python \"{this_file}\" MergeGenomeFasta --updir {upid}".format(
                            this_file=os.path.realpath(__file__),
                            updir=updir)

                elif self.tool.lower() == 'fasplit':
                    if self.lsf is True:
                        cmd = "bsub -M %s -R \"rusage[mem=%s,tmp=%s]\" -o %s -e %s -u %s -Ep \"rm -rf luigi\" " \
                              "-g %s python %s SplitGenomeFasta --updir %s --upid %s" % (gc.MEM, gc.MEM, gc.TMP_MEM,
                                                                      os.path.join(updir, "split.out"),
                                                                      os.path.join(updir, "split.err"),
                                                                      gc.USER_EMAIL, gc.SRCH_GROUP,
                                                                      os.path.realpath(__file__), updir, upid)

                    else:
                        cmd = "python \"{this_file}\" SplitGenomeFasta --updir {upid} --upid {upid}".format(
                            this_file=os.path.realpath(__file__),
                            updir=updir,
                            upid=upid)

                elif self.tool.lower() == 'faclean':
                    if self.lsf is True:
                        cmd = "bsub -M %s -R \"rusage[mem=%s,tmp=%s]\" -o %s -e %s -u %s -Ep \"rm -rf luigi\" " \
                              "-g %s python %s RewriteCleanFasta --updir %s --upid %s" % (gc.MEM, gc.MEM, gc.TMP_MEM,
                                                                                         os.path.join(updir,
                                                                                                      "clean.out"),
                                                                                         os.path.join(updir,
                                                                                                      "clean.err"),
                                                                                         gc.USER_EMAIL, gc.SRCH_GROUP,
                                                                                         os.path.realpath(__file__),
                                                                                         updir, upid)

                    else:
                        cmd = "python \"{this_file}\" RewriteCleanFasta --updir {upid} --upid {upid}".format(
                            this_file=os.path.realpath(__file__),
                            updir=updir,
                            upid=upid)

                elif self.tool.lower() == 'tblmerge':
                    if self.lsf is True:
                        cmd = "bsub -M %s -R \"rusage[mem=%s,tmp=%s]\" -o %s -e %s -u %s -Ep \"rm -rf luigi\" " \
                              "-g %s python %s MergeGenomeTBLOUT --updir %s --upid %s" % (gc.MEM, gc.MEM, gc.TMP_MEM,
                                                                               os.path.join(updir, "tbl_merge.out"),
                                                                               os.path.join(updir, "tbl_merge.err"),
                                                                               gc.USER_EMAIL, gc.SRCH_GROUP,
                                                                               os.path.realpath(__file__),
                                                                                          updir, upid)
                    else:
                        cmd = "python \"{this_file}\" MergeGenomeTBLOUT --updir {upid} --upid {upid}".format(
                            this_file=os.path.realpath(__file__),
                            updir=updir,
                            upid=upid)

            subprocess.call(cmd, shell=True)

            cmd = ''

# -----------------------------------------------------------------------------


if __name__ == '__main__':
    luigi.run()
