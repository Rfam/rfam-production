#/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Tool to generate taxonomy entries based on a list of tax ids

requirements:
  DockerRequirement:
    dockerPull: rfam/rfam-production:cwl-infernal

inputs:
  fasta_file:
    type: File
    inputBinding:
      position: 1
  tax_id:
    type: int
    inputBinding:
      position: 2
  source: 
    type: string
    inputBinding:
      position: 3

baseCommand: [python, '/rfam/rfam-production/support/fasta2rfamseq.py']

outputs:
  rfamseq_entries:
    type: stdout

stdout: test.rfamseq

$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
