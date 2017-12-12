import os
import luigi
import json
import subprocess

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
    updir = luigi.parameter()

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
    project_dir = luigi.Parameter()
    lsf = luigi.BoolParameter(default=True,
                              description="If used then run on lsf, otherwise run locally")

    def run(self):
        """
        Call load_upid_gca_file to export all upid_gca accession pairs in
        json format.
        """
        id_pairs = None
        upid_gca_file_loc = os.path.join(self.project_dir, "upid_gca_dict.json")
        upid_gca_fp = open(upid_gca_file_loc, 'r')
        accessions = json.load(upid_gca_fp)
        upid_gca_fp.close()

        upids = accessions.keys()

        cmd = ''
        for upid in upids:

            subdir = os.path.join(self.project_dir, upid[-3:])
            updir = os.path.join(subdir, updir)

            if os.path.exists(updir):

                if self.lsf is True:
                    # TODO launch an lsf job to merge the genome fasta files
                    pass

                else:
                    cmd = "python \"{this_file}\" MergeGenome --updir {upid}".format(
                        this_file=os.path.realpath(__file__),
                        updir=updir)

                    subprocess.call(cmd, shell=True)

            cmd = ''

# -----------------------------------------------------------------------------


if __name__ == '__main__':

    luigi.run()




