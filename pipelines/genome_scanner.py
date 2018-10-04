"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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
Pipeline to search each genome individually. This pipeline also parallelises the step
of splitting and indexing the genome files in smaller chunks
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import luigi
import subprocess

from config import gen_config as gc
from utils import genome_search_utils as gsu

# add parent directory to path
if __name__ == '__main__' and __package__ is None:
    os.sys.path.append(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ----------------------------------TASKS--------------------------------------


class ScanGenome(luigi.Task):
    """
    Add comment here
    """
    updir = luigi.Parameter()
    upid = luigi.Parameter()
    tool = luigi.Parameter()
    lsf = luigi.Parameter()

    def run(self):
        """
        Add comment here
        """
        if self.lsf is True:
            gsu.single_genome_scan_from_download_directory(self.updir,
                                                       self.upid, tool=self.tool)

# -----------------------------------------------------------------------------


class GenomeSearchEngine(luigi.Task):
    """
    Launches genome search directly from the download directory which should
    have the following structure project_dir/xxx/UPYYYYYYxxx
    """
    project_dir = luigi.Parameter(description="Genome download project directory")
    genome_list = luigi.Parameter(default=None,
                                  description="A list of upids to search for ncRNAs")
    tool = luigi.Parameter(default=None,
                                  description="Infernal search tool (cmsearch/cmscan)")
    lsf = luigi.BoolParameter(default=False,
                              description="Run pipeline on lsf, otherwise run locally")

    def run(self):
        """
        Add comment here
        """
        # load upids
        upid_fp = open(self.genome_list, 'r')
        upids = [x.strip() for x in upid_fp]
        upid_fp.close()

        for upid in upids:
            # get updir location
            subdir = os.path.join(self.project_dir, upid[-3:])
            updir = os.path.join(subdir, upid)

            yield ScanGenome(updir, upid, self.tool.lower(), self.lsf)

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # defining pipeline's main task
    luigi.run()