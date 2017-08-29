#!/bin/bash

# Copyright [2009-2017] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


######################################################
## Launch XML export of Rfam data using LSF cluster ##
######################################################

set -e

usage = "Usage: lsf_rfam_xml_dumper.sh /path/to/output"

if [ "$#" -ne 1 ]
then
  echo $usage
  exit 1
fi

dir=$1

# activate virtual environment
source env/bin/activate
export PYTHONPATH=`pwd`

# check that output directories exist
mkdir -p $dir/families
mkdir -p $dir/motifs
mkdir -p $dir/clans
mkdir -p $dir/genomes
mkdir -p $dir/full_region

# launch jobs
echo "bsub python scripts/export/rfam_xml_dumper.py --type M --out $dir/motifs/"
echo "bsub python scripts/export/rfam_xml_dumper.py --type C --out $dir/clans/"
echo "bsub -M 16384 -R "rusage[mem=16384]" python scripts/export/rfam_xml_dumper.py --type G --out $dir/genomes/"
echo "bsub -M 16384 -R "rusage[mem=16384]" python scripts/export/rfam_xml_dumper.py --type F --out $dir/families/"

# -F (file size) is required to allow creation of large files
echo "bsub -M 16384 -R "rusage[mem=16384]" -F 1000000 python scripts/export/rfam_xml_dumper.py --type R --out $dir/full_region"
