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

"""
Download reference genome sequences.

Usage:

python genome_downloader.py GenomesDownloadEngine

To see more information about available options run
python genome_downloader.py GenomesDownloadEngine --help

Running the Central scheduler on the cluster:

1. Get an interactive node using bsub
    bsub -Is $SHELL

2. Start the central scheduler
    luigid

3. ssh to the interactive node

4. Run the pipeline
python genome_downloader.py GenomesDownloadEngine --project-name <project name>
--upid-file <upid_gca file path> --lsf

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


class CopyFileFromFTP(luigi.Task):
    """
    Download a file for a specific accession.
    """
    accession = luigi.Parameter()
    prot_dir = luigi.Parameter()

    def run(self):
        """
        Download ENA file.
        """
        # need to parametrise file format
        if self.accession[0:3] == "GCA":
            gflib.copy_gca_report_file_from_ftp(self.accession, self.prot_dir)
        else:
            gflib.copy_wgs_set_from_ftp(self.accession, self.prot_dir)

    def output(self):
        """
        Check if file exists.
        """

        filename = ''
        if self.genome_acc[0:3] == "GCA":
            filename = self.accession + "_sequence_report.txt"
        else:
            # only use the 5 first characters of the WGS string
            filename = self.accession[0:6] + ".fasta.gz"

        return luigi.LocalTarget(os.path.join(self.prot_dir,
                                              filename))

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
        Setup proteome directory.
        """

        # get project subdir index where the genome will be downloaded
        sub_dir_index = self.upid[8:]

        # construct the path for the new proteome directory
        self.upid_dir = os.path.join(
            os.path.join(self.project_dir, sub_dir_index), self.upid)

        # finally, create the directory if it does not exist
        if not os.path.exists(self.upid_dir):
            os.mkdir(self.upid_dir)
            os.chmod(self.upid_dir, 0777)

    def run(self):
        """
        Fetch genome accessions and call DownloadFile Task.
        """
        # generate directory for this proteome
        self.setup_proteome_dir()

        # fetch proteome accessions, this will also copy GCA file if available
        genome_accessions = gflib.get_genome_unique_accessions(self.upid, self.upid_dir)

        other_accessions = None
        wgs_set = None

        if genome_accessions["GCA"] != -1:
            # 1. check for assembly report file
            # This list is going to be empty
            other_accessions = genome_accessions["OTHER"]

            # fetch wgs set from ENA
            if len(other_accessions) == 0 or genome_accessions["WGS"] == 0:
                wgs_set = gflib.extract_wgs_acc_from_gca_xml(genome_accessions["GCA"])

        elif genome_accessions["WGS"] != -1 and genome_accessions["GCA"] == -1 or wgs_set is not None:
            # First copy WGS set in upid dir
            yield CopyFileFromFTP(wgs_set, self.upid_dir)

        # this should be done in all cases
        # download genome accessions in proteome directory
        if len(other_accessions) > 0:
            for acc in other_accessions:
                test = yield DownloadFile(ena_acc=acc, prot_dir=self.upid_dir)
                self.acc_data.append(test)

            self.db_entries[self.upid] = self.acc_data

        # merge and validate need to check final UPXXXXXXXXX.fa exists

# -----------------------------------------------------------------------------


class GenomesDownloadEngine(luigi.Task):
    """
    Initialize project environment parse Uniprot's the UPID_GCA file
    and call DownloadGenome task for each available reference proteome.
    """
    project_name = luigi.Parameter(default="genome-download",
        description="Files will be stored in LOC_PATH/project_name folder")
    upid_file = luigi.Parameter(default=None,
        description="UniProt upid_gca file. Use UniProt REST API by default")
    lsf = luigi.BoolParameter(default=True,
        description="If specified then run on lsf, otherwise run locally")

    proj_dir = ''

    def initialize_project(self):
        """
        Create main project directory and generate species domain subdirectories
        for proteome directories will be located.
        """
        location = ""

        if self.lsf is True:
            location = gc.RFAM_HPS_LOC
        else:
            location = gc.LOC_PATH

        self.proj_dir = os.path.join(location, self.project_name)

        # create project directory
        if not os.path.exists(self.proj_dir):
            os.mkdir(self.proj_dir)
            os.chmod(self.proj_dir, 0777)

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

        if self.lsf is True:
            location = gc.RFAM_HPS_LOC
        else:
            location = gc.LOC_PATH

        self.proj_dir = os.path.join(location, self.project_name)

        # load prot-gca pairs
        upid_gca_fp = self.input()["upid_gcas"].open('r')
        upid_gca_pairs = json.load(upid_gca_fp)
        upid_gca_fp.close()

        cmd = ''

        for upid in upid_gca_pairs.keys():

            sub_dir_index = upid[8:]
            sub_dir_loc = os.path.join(self.proj_dir, sub_dir_index)

            if not os.path.exists(sub_dir_loc):
                os.mkdir(sub_dir_loc)
                os.chmod(sub_dir_loc, 0777)

            # generate an lsf command
            if self.lsf is True:
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
