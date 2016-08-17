"""
Copyright [2009-2016] EMBL-European Bioinformatics Institute
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
Download reference genome sequences.

Usage:

python genome_downloader.py GenomesDownloadEngine --project-name <project name>
--upid-file <upid_gca file path> --lsf <Bool> [--local-scheduler]

--project-name: A string indicating the project's name
--upid-file: The path to Uniprot's upid_gca file
--lsf (optional): if specified, run on the cluster. By default, run locally.
--local-scheduler to run the pipeline with the local scheduler

Running the Central scheduler on the cluster:

1. Get an interactive node using bsub
    bsub -Is $SHELL

2. Start the central scheduler
    luigid

3. ssh to the interactive node

4. Run the pipeline
python genome_downloader.py GenomesDownloadEngine --project-name <project name>
--upid-file <upid_gca file path> --lsf True

luigi.cfg - See Luigi's Configuration
(http://luigi.readthedocs.io/en/stable/configuration.html)
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import luigi
import json
import subprocess


from config import gen_config as gc
from scripts.export.genomes import genome_fetch as gflib

# add parent directory to path
if __name__ == '__main__' and __package__ is None:
    os.sys.path.append(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ----------------------------------TASKS--------------------------------------


class UPAccessionLoader(luigi.Task):
    """
    Generate a dictionary of UPID-GCA pairs and the genome's domain by either
    parsing Uniprot's upid_gca file or using Uniprot's REST API.
    """
    project_dir = luigi.Parameter()
    upid_gca_file = luigi.Parameter()

    def run(self):
        """
        Call load_upid_gca_file to export all upid_gca accession pairs in
        json format.
        """
        id_pairs = None

        # load accessions from upid_gca file
        if self.upid_gca_file is None:
            id_pairs = gflib.load_upid_gca_pairs()
        # fetch accessions from Uniprots REST API
        else:
            id_pairs = gflib.load_upid_gca_file(self.upid_gca_file)

        outfile = self.output().open('w')
        json.dump(id_pairs, outfile)
        outfile.close()

    def output(self):
        """
        Check if the json file exists
        """
        out_file = os.path.join(self.project_dir, "upid_gca_dict.json")
        return luigi.LocalTarget(out_file)

# -----------------------------------------------------------------------------


class DownloadFile(luigi.Task):
    """
    Download a file for a specific accession.
    """
    ena_acc = luigi.Parameter()
    prot_dir = luigi.Parameter()

    def run(self):
        """
        Download ENA file.
        """
        # need to parametrise file format
        gflib.fetch_ena_file(self.ena_acc, "fasta", self.prot_dir)

    def output(self):
        """
        Check if file exists.
        """
        return luigi.LocalTarget(os.path.join(self.prot_dir,
                                              self.ena_acc + ".fa.gz"))

# -----------------------------------------------------------------------------


class DownloadGenome(luigi.Task):
    """
    Download all genome related accessions.
    """
    upid = luigi.Parameter()  # uniprot's ref proteome id
    gca_acc = luigi.Parameter()  # proteome corresponding GCA accession
    project_dir = luigi.Parameter()
    domain = luigi.Parameter()
    upid_dir = None
    acc_data = []
    db_entries = {}

    def setup_proteome_dir(self):
        """
        Setup project directory.
        """
        self.upid_dir = os.path.join(
            os.path.join(self.project_dir, self.domain), self.upid)
        if not os.path.exists(self.upid_dir):
            os.mkdir(self.upid_dir)
            os.chmod(self.upid_dir, 0777)

    def run(self):
        """
        Fetch genome accessions and call DownloadFile Task.
        """
        # generate directory for this proteome
        self.setup_proteome_dir()

        # fetch accessions
        accessions = gflib.fetch_genome_accessions(
            self.upid, str(self.gca_acc))

        # download genome accessions in proteome directory
        for acc in accessions:
            test = yield DownloadFile(ena_acc=acc, prot_dir=self.upid_dir)
            self.acc_data.append(test)

        self.db_entries[self.upid] = self.acc_data

# -----------------------------------------------------------------------------


class GenomesDownloadEngine(luigi.Task):
    """
    Initialize project environment parse Uniprot's the UPID_GCA file
    and call DownloadGenome task for each available reference proteome.
    """
    project_name = luigi.Parameter()  # a project name for the download
    upid_file = luigi.Parameter(default=None)  # uniprot's upid_gca file
    lsf = luigi.BoolParameter(default=False)  # if True run on lsf otherwise run locally
    proj_dir = ''

    def initialize_project(self):
        """
        Create main project directory and generate species domain subdirectories
        for proteome directories will be located.
        """
        location = ""

        if self.lsf == "True":
            location = gc.RFAM_GPFS_LOC
        else:
            location = gc.LOC_PATH

        self.proj_dir = os.path.join(location, self.project_name)

        # create project directory
        if not os.path.exists(self.proj_dir):
            os.mkdir(self.proj_dir)
            os.chmod(self.proj_dir, 0777)

        # create domain subdirs
        for domain in gc.DOMAINS:
            dom_dir = os.path.join(self.proj_dir, domain)
            if not os.path.exists(dom_dir):
                os.mkdir(dom_dir)
                os.chmod(dom_dir, 0777)

    def requires(self):
        """
        Initialize project directory and call UPAccessionLoader task to load
        upid-gca accession pairs into a json object.
        """
        self.initialize_project()

        return {
            "upid_gcas": UPAccessionLoader(project_dir=self.proj_dir,
                                           upid_gca_file=self.upid_file)
        }

    def run(self):
        """
        Execute the pipeline.
        """
        location = ""

        if self.lsf == "True":
            location = gc.RFAM_GPFS_LOC
        else:
            location = gc.LOC_PATH

        self.proj_dir = os.path.join(location, self.project_name)

        # load prot-gca pairs
        upid_gca_fp = self.input()["upid_gcas"].open('r')
        upid_gca_pairs = json.load(upid_gca_fp)
        upid_gca_fp.close()

        cmd = ''

        for upid in upid_gca_pairs.keys():

            # generate an lsf command
            if self.lsf == "True":
                cmd = gflib.lsf_cmd_generator(upid, upid_gca_pairs[upid]["GCA"],
                                              upid_gca_pairs[upid]["DOM"],
                                              os.path.realpath(__file__),
                                              self.proj_dir)

            # run the pipeline locally
            else:
                cmd = "python \"{this_file}\" DownloadGenome --upid {upid} " \
                      "--gca-acc {gca_acc} --project-dir {proj_dir} " \
                      "--domain {domain}".format(
                          this_file=os.path.realpath(__file__),
                          upid=upid,
                          gca_acc=upid_gca_pairs[upid]["GCA"],
                          proj_dir=self.proj_dir,
                          domain=upid_gca_pairs[upid]["DOM"])

            subprocess.call(cmd, shell=True)

            cmd = ''

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # defining main pipeline's main task
    luigi.run()
