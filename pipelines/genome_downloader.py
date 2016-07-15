'''
Created on 7 Jul 2016

@author: ikalvari

Usage:

python genome_downloader.py GenomesDownloadEngine --project-name <project name>
--upid-file <upid_gca file path> --lsf <Bool> [--local-scheduler]

--project-name: A string indicating the project's name
--upid-file: The path to Uniprot's upid_gca file
--lsf: False (local), True (cluster)
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
'''

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


class ProjectInitializer(luigi.Task):  # redundant

    '''
    Forms the project directory structure
    '''

    priority = 10  # needs to be higher than DownloadGenomes

    location = luigi.Parameter()  # rfam's gpfs path on cluster
    project_name = luigi.Parameter(default='gen_downloads')
    project_dir = ''

    def requires(self):
        '''
        Return empty list to avoid pipeline kick off
        '''
        return []

    def run(self):

        dom_dir = None
        self.project_dir = os.path.join(self.location, self.project_name)

        # create project directory
        if not os.path.exists(self.project_dir):
            os.mkdir(self.project_dir)
            os.chmod(self.project_dir, 0777)

        # create domain subdirs
        for domain in gc.DOMAINS:
            dom_dir = os.path.join(self.project_dir, domain)
            if not os.path.exists(dom_dir):
                os.mkdir(dom_dir)
                os.chmod(dom_dir, 0777)

            dom_dir = None

# -----------------------------------------------------------------------------


class UPAccessionLoader(luigi.Task):

    '''
    Reads Uniprot's UPID_GCA file and generates a dictionary of UPID-GCA pairs
    and the genome's domain.
    '''

    # the 1st task to be executed
    priority = 9

    project_dir = luigi.Parameter()  # maybe implement an initialization Task
    upid_gca_file = luigi.Parameter()

    def run(self):

        id_pairs = gflib.load_upid_gca_file(self.upid_gca_file)

        outfile = self.output().open('w')
        json.dump(id_pairs, outfile)
        outfile.close()

    def output(self):
        '''
        Exports upid_gca accessions in json format
        '''
        out_file = os.path.join(self.project_dir, "upid_gca_dict.json")
        return luigi.LocalTarget(out_file)

# -----------------------------------------------------------------------------


class DownloadFile(luigi.Task):

    '''
    Downloads a file for a specific accession
    '''

    priority = 8

    ena_acc = luigi.Parameter()
    prot_dir = luigi.Parameter()

    def run(self):
        '''
        Download ENA file
        '''

        # need to parametrise file format
        gflib.fetch_ena_file(self.ena_acc, "fasta", self.prot_dir)

    def output(self):
        '''
        Check if file exists
        '''
        return luigi.LocalTarget(os.path.join(self.prot_dir,
                                              self.ena_acc + '.fa.gz'))

# -----------------------------------------------------------------------------


class DownloadGenome(luigi.Task):

    upid = luigi.Parameter()  # uniprot's ref proteome id
    gca_acc = luigi.Parameter()  # proteome corresponding GCA accession
    project_dir = luigi.Parameter()
    domain = luigi.Parameter()
    upid_dir = None
    acc_data = []
    db_entries = {}

    priority = 7

    def setup_proteome_dir(self):
        '''
            Setup project directory
        '''
        self.upid_dir = os.path.join(
            os.path.join(self.project_dir, self.domain), self.upid)
        if not os.path.exists(self.upid_dir):
            os.mkdir(self.upid_dir)
            os.chmod(self.upid_dir, 0777)

    def run(self):
        '''
        Fetch genome accessions and call DownloadFile Task
        '''

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
    '''
    This Task will initialize project environment parse Uniprot's the UPID_GCA
    file and call DownloadGenome Task for each available reference proteome
    '''

    priority = 6

    project_name = luigi.Parameter()  # a project name for the download
    upid_file = luigi.Parameter()  # uniprot's upid_gca file
    lsf = luigi.Parameter()  # if True run it on lsf otherwise locally
    proj_dir = ''

    def initialize_project(self):

        location = ""

        if self.lsf == "True":
            location = gc.RFAM_GPFS_LOC
        else:
            location = gc.LOC_PATH

        self.proj_dir = os.path.join(location, self.project_name)

        dom_dir = None

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

            dom_dir = None

    def requires(self):
        '''
        Initializes project directory and calls Task to create
        '''

        self.initialize_project()

        return {
            'upid_gcas': UPAccessionLoader(project_dir=self.proj_dir,
                                           upid_gca_file=self.upid_file)
        }

    def run(self):
        '''
        Executes the pipeline
        '''

        location = ""

        if self.lsf == "True":
            location = gc.RFAM_GPFS_LOC
        else:
            location = gc.LOC_PATH

        self.proj_dir = os.path.join(location, self.project_name)

        # load prot-gca pairs
        upid_gca_fp = self.input()['upid_gcas'].open('r')
        upid_gca_pairs = json.load(upid_gca_fp)
        upid_gca_fp.close()

        cmd = ''

        for upid in upid_gca_pairs.keys():

            # generate an lsf command
            if self.lsf == "True":
                cmd = gflib.lsf_cmd_generator(upid, upid_gca_pairs[upid]['GCA'],
                                              upid_gca_pairs[upid]['DOM'],
                                              os.path.realpath(__file__),
                                              self.proj_dir)

            # run the pipeline locally
            else:
                cmd = "python \"{this_file}\" DownloadGenome --upid {upid} " \
                    "--gca-acc {gca_acc} --project-dir {proj_dir} " \
                    "--domain {domain}".format(
                        this_file=os.path.realpath(__file__),
                        upid=upid,
                        gca_acc=upid_gca_pairs[upid]['GCA'],
                        proj_dir=self.proj_dir,
                        domain=upid_gca_pairs[upid]['DOM'])

            subprocess.call(cmd, shell=True)

            cmd = ''

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    luigi.run()
