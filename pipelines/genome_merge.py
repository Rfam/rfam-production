import os
import luigi
import json
import subprocess

from config import gen_config as gc
from support import merge_fasta as mf

# add parent directory to path
if __name__ == '__main__' and __package__ is None:
    os.sys.path.append(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ----------------------------------TASKS--------------------------------------


class MergeGenome(luigi.Task):
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
        genome_fasta = os.path.join(self.updir, upid+'.fa')

        return luigi.LocalTarget(genome_fasta)

# -----------------------------------------------------------------------------


class GenomeMergeEngine(luigi.Task):
    """
    Task to load accessions from file
    """
    project_dir = luigi.Parameter(description="Genome download project directory")
    genome_list = luigi.Parameter(default=None,
                                  description="A list of upids for which to merge the sequence files")
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
                if self.lsf is True:

                    cmd = """
                          bsub -M %s -R "rusage[mem=%s,tmp=%s]" -o %s -e %s -u %s
                          -Ep "rm -rf luigi"
                          -g %s
                          python %s MergeGenome --updir %s
                          """ % (gc.MEM, gc.MEM, gc.TMP_MEM,
                              os.path.join(updir, "merge.out"),
                              os.path.join(updir, "merge.err"),
                              gc.USER_EMAIL, gc.LSF_GEN_GROUP,
                              os.path.realpath(__file__), updir)

                else:
                    cmd = "python \"{this_file}\" MergeGenome --updir {upid}".format(
                        this_file=os.path.realpath(__file__),
                        updir=updir)

                subprocess.call(cmd, shell=True)

            cmd = ''

# -----------------------------------------------------------------------------


if __name__ == '__main__':

    luigi.run()




