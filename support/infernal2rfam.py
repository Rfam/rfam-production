"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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
Script to convert infernal output files to full_region
"""

# -------------------------------------------------------------------------

import sys
from utils import infernal_utils as iu

# -------------------------------------------------------------------------
if __name__ == '__main__':

    input_file = sys.argv[1]

    if "--tbl" in sys.argv or "--tblout" in sys.argv:
            iu.tblout_to_full_region(input_file, dest_dir=None)

    elif "-o" in sys.argv or "--out" in sys.argv:
        iu.infernal_to_full_region(input_file, dest_dir=None, filename=None)

    else:
        print "\nWrong input!\n"

        print "Usage infernal_file [-o|--tbl]\n"

        print "\n-o (--out): parse infernal output format"
        print "\n--tbl (--tblout): parse infernal tblout format"