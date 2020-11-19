#!/usr/bin/python

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

import os
import sys
import json
from scripts.export.genomes import fetch_gen_metadata as fgm

# -------------------------------------------------------------------------


def wet_wgs_genome_entry(upid, domain, gca_acc=''):

    entry_value_list = []

    uniprot_metadata = fgm.dump_uniprot_genome_metadata(upid, domain)["fields"]

    entry_value_list.append(upid)
    # append two empty strings for gca_acc
    if gca_acc != '':
        entry_value_list.append(gca_acc)
        entry_value_list.append(gca_acc.partition('.')[2])
    else:
        entry_value_list.append('\N')
        entry_value_list.append('\N')

    if uniprot_metadata["wgs_acc"] is not None:
        entry_value_list.append(uniprot_metadata["wgs_acc"])
        # wgs version - convert to int to remove any leading zeros and back to str
        # to concat with the rest of the string
        entry_value_list.append(str(int(uniprot_metadata["wgs_acc"][4:6])))
    else:
        # append an empty character for all wgs fields
        entry_value_list.append('\N')
        entry_value_list.append('\N')

    # assembly name
    if "assembly_name" not in uniprot_metadata:
        entry_value_list.append('\N')
    else:
        if uniprot_metadata["assembly_name"] is not None:
            entry_value_list.append(uniprot_metadata["assembly_name"])
        else:
            entry_value_list.append('\N')

    # assembly level
    if "assembly_level" not in uniprot_metadata:
        entry_value_list.append('\N')
    else:
        if uniprot_metadata["assembly_level"] is not None:
            entry_value_list.append(uniprot_metadata["assembly_level"])
        else:
            entry_value_list.append('\N')

    # study_ref value
    entry_value_list.append('\N')

    # genome description
    entry_value_list.append(uniprot_metadata["description"].replace('\n', ' '))

    # genome length fields
    # entry_value_list.append(str(ena_metadata["total_length"]))
    entry_value_list.append('\N')
    # entry_value_list.append(str(ena_metadata["ungapped_length"]))
    entry_value_list.append('\N')
    # entry_value_list.append(str(ena_metadata["circular"]))
    entry_value_list.append('\N')
    entry_value_list.append(str(uniprot_metadata["ncbi_id"]))
    entry_value_list.append(uniprot_metadata["scientific_name"])

    common_name = '\N'
    entry_value_list.append(common_name)
    entry_value_list.append(uniprot_metadata["kingdom"])

    # setting the fields num_rfam_regions and num_families to NULL
    # entry_value_list.append('\N')
    entry_value_list.append('\N')
    # entry_value_list.append('\N')
    entry_value_list.append('\N')

    entry_value_list.append(str(uniprot_metadata["is_reference"]))
    entry_value_list.append(str(uniprot_metadata["is_representative"]))

    entry_value_list.append(uniprot_metadata["created"])
    entry_value_list.append(uniprot_metadata["updated"])

    new_entry = '\t'.join(entry_value_list)

    return new_entry

# -------------------------------------------------------------------------


def get_gca_genome_entry(upid, assembly_acc, domain):

    entry_value_list = []
    new_entry = ''

    ena_dict = fgm.fetch_gca_data(upid, assembly_acc, domain)
    # add the upid
    entry_value_list.append(upid)

    # check for any metadata from ENA
    if "fields" in ena_dict:
        ena_metadata = ena_dict["fields"]
        uniprot_metadata = fgm.dump_uniprot_genome_metadata(upid, domain)["fields"]

        # build list
        entry_value_list.append(assembly_acc)
        entry_value_list.append(assembly_acc.partition('.')[2])

        if uniprot_metadata["wgs_acc"] is not None:
            entry_value_list.append(uniprot_metadata["wgs_acc"])
            # wgs version
            entry_value_list.append(str(int(uniprot_metadata["wgs_acc"][4:6])))
        else:
            # append an empty character for all wgs fields
            entry_value_list.append('\N')
            entry_value_list.append('\N')

        # assembly_name
        if ena_metadata["assembly_name"] is not None:
            entry_value_list.append(ena_metadata["assembly_name"])
        else:
            entry_value_list.append('\N')

        # assembly_level
        if ena_metadata["assembly_level"] is not None:
            entry_value_list.append(ena_metadata["assembly_level"])
        else:
            entry_value_list.append('\N')

        if ena_metadata["study_ref"] is not None:
            entry_value_list.append(ena_metadata["study_ref"])
        else:
            entry_value_list.append('\N')

        # genome description
        entry_value_list.append(ena_metadata["description"])

        # genome length fields
        entry_value_list.append(str(ena_metadata["total_length"]))
        entry_value_list.append(str(ena_metadata["ungapped_length"]))

        entry_value_list.append(str(ena_metadata["circular"]))
        entry_value_list.append(str(ena_metadata["ncbi_id"]))
        entry_value_list.append(ena_metadata["scientific_name"])

        common_name = ''
        if "common_name" in ena_metadata:
            common_name = ena_metadata["common_name"]

        entry_value_list.append(common_name)
        entry_value_list.append(ena_metadata["kingdom"])

        # setting the fields num_rfam_regions and num_families to NULL
        # entry_value_list.append('\N')
        entry_value_list.append('\N')
        # entry_value_list.append('\N')
        entry_value_list.append('\N')

        entry_value_list.append(str(uniprot_metadata["is_reference"]))
        entry_value_list.append(str(uniprot_metadata["is_representative"]))

        entry_value_list.append(uniprot_metadata["created"])
        entry_value_list.append(uniprot_metadata["updated"])

        new_entry = '\t'.join(entry_value_list)

    else:
        if assembly_acc != -1:
            new_entry = wet_wgs_genome_entry(upid, domain, gca_acc=assembly_acc)
        else:
            new_entry = wet_wgs_genome_entry(upid, domain, '')

    return new_entry

# -------------------------------------------------------------------------


if __name__ == '__main__':

    # a upid_gca file generated from Uniprot
    upid_gca_file = sys.argv[1]

    # load accessions
    fp = open(upid_gca_file, 'r')
    accessions = json.load(fp)
    fp.close()

    for upid in accessions.keys():

        if accessions[upid]["GCA"] != -1:
            assembly_acc = accessions[upid]["GCA"]
            domain = accessions[upid]["DOM"]
            new_entry = get_gca_genome_entry(upid, assembly_acc, domain)

        else:
            new_entry = wet_wgs_genome_entry(upid,
                                             accessions[upid]["DOM"], '')

        print new_entry