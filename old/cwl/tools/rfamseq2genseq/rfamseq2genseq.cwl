cwlVersion: v1.0
class: CommandLineTool
label: Script to map an rfamseq file to a genseq file

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: rfam/rfam-production:cwl-infernal

inputs:
  rfamseq:
    type: File
    inputBinding:
      position: 2

baseCommand: [python, /rfam/rfam-production/cwl/tools/rfamseq2genseq/rfamseq2genseq_single.py]

outputs:
  genseq:
    type: File
    outputBinding:
      glob: $(inputs.rfamseq.nameroot + '.genseq')




$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
