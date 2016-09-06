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

import sys
from utils import infernal_utils as iu

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    infernal_output = sys.argv[1]
    out_dir = sys.argv[2]
    ss_notation = None

    if len(sys.argv) > 3:
        ss_notation = sys.argv[3]
        iu.infernal_output_parser(infernal_output, out_dir, ss_notation)

    else:
        iu.infernal_output_parser(infernal_output, out_dir, ss_notation="wuss")
