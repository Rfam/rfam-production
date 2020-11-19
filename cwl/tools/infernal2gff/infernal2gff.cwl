cwlVersion: v1.0
class: CommandLineTool
label: index a sequence file for use by esl-sfetch
doc: "https://github.com/EddyRivasLab/easel"

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMax: 1
    ramMin: 1024  # just a default, could be lowered

hints:
  DockerRequirement:
    dockerPull: rfam/rfam-production:cwl-infernal

inputs:
  infernal_file:
    type: File
    inputBinding:
      position: 2

  max_e_val:
    label: 'maximum E-value to include'
    type: int?
    inputBinding:
      prefix: -E

  min_bit_score:
    label: 'minimum bit score to include'
    type: int?
    inputBinding:
      prefix: -T

  cmscan:
    label: 'set if infernal file is in cmscan format'
    type: boolean?
    default: false
    inputBinding:
      prefix: '--cmscan'

  fmt2:
    label: 'set if infernal file used --fmt 2 option'
    type: boolean?
    default: false
    inputBinding:
      prefix: '--fmt2'

  all_info:
    label: 'output all info in `attributes` column [default: E-value]'
    type: boolean?
    default: false
    inputBinding:
      prefix: '--all'

  no_info:
    label: 'output no info in `attributes` column [default: E-value]'
    type: boolean?
    default: false
    inputBinding:
      prefix: '--none'

baseCommand: [ perl, /src/infernal2gff.pl ]

stdout: out.gff  # to aid cwltool's cache feature

outputs:
  gff:
    type: File
    outputBinding:
      glob: out.gff
      outputEval: ${self[0].basename=inputs.infernal_file.basename + '.gff'; return self;}


$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
