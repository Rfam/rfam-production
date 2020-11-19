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

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import json

# -----------------------------------------------------------------------------

def convert_genseq_to_rfamseq(genseq_dump):
    """
    Loads genseq json object and generates a rfamseq dump in txt format

    genseq_dump: This can be a directory or a genseq .json file

    return: void
    """

    if os.path.isfile(genseq_dump):

        rfamseq_entries = genseq_file_to_rfamseq(genseq_dump)

        for entry in rfamseq_entries:
            print '\t'.join(entry)

    elif os.path.isdir(genseq_dump):
        json_files = os.listdir(genseq_dump)

        for json_file in json_files:

            genseq_file_loc = os.path.join(genseq_dump, json_file)
            rfamseq_entries = genseq_file_to_rfamseq(genseq_file_loc)

            for entry in rfamseq_entries:
                print '\t'.join(entry)

    else:
        sys.exit("Wrong input! Please provide a json file or directory")

# -----------------------------------------------------------------------------

def genseq_file_to_rfamseq(genseq_json_dump):
    """
    Loads genseq json object and generates a rfamseq dump in txt format

    genseq_dump:  A genseq .json file generated from metadata export

    return: void
    """

    genseq_fp = open(genseq_json_dump, 'r')
    genseq_entries = json.load(genseq_fp)
    rfamseq_entries = []

    for genseq_dict in genseq_entries:

        # get fields
        fields = genseq_dict["fields"]
        genseq_fp.close()

        # initializing list with pk
        rfamseq_acc = genseq_dict["pk"]
        rfamseq_attributes = [rfamseq_acc]

        accession = ''
        version = ''
        if rfamseq_acc.find('.'):
            accession = rfamseq_acc.partition('.')[0]
            version = rfamseq_acc.partition('.')[2]
        else:
            accession = rfamseq_acc
            version = fields["seq_version"]

        rfamseq_attributes.append(accession)
        rfamseq_attributes.append(version)

        # ncbi_id
        rfamseq_attributes.append(str(fields["ncbi_id"]))

        # setting mol_type to other RNA for all new sequneces
        rfamseq_attributes.append("genomic DNA")

        # length
        rfamseq_attributes.append(str(fields["seq_length"]))

        # description
        rfamseq_attributes.append(fields["description"])

        # set previous accession and source
        previous_accession = ''
        source = "UNIPROT; ENA"
        rfamseq_attributes.append(previous_accession)
        rfamseq_attributes.append(source)

        rfamseq_entries.append(rfamseq_attributes)

    return rfamseq_entries

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    genseq_input = sys.argv[1]
    convert_genseq_to_rfamseq(genseq_input)
