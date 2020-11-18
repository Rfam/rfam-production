#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Tool to generate taxonomy entries based on a list of tax ids
baseCommand: [python, '/Users/ikalvari/RfamWorkspace/rfam-bh19/rfam-production/scripts/export/generate_tax_data.py']
inputs:
  taxid_list:
    type: File
    inputBinding:
      position: 1
outputs:
  taxonomy_data:
    type: stdout

stdout: taxonomy.txt

$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
