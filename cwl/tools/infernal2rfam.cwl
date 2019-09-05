cwlVersion: v1.0
class: CommandLineTool
label: Script to map an infernal file to a rfam file

hints:
  DockerRequirement:
    dockerPull: miguelboland/rfam-production:latest

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.infernal_file)


inputs:
  infernal_file:
    type: File
    inputBinding:
      position: 1

  infernal_tbl:
    label: 'Read infernal tbl input'
    type: boolean
    default: false
    inputBinding:
      prefix: --tbl
      position: 2

  infernal_out:
    label: 'parse infernal output format'
    type: boolean
    default: false
    inputBinding:
      prefix: --out
      position: 3

baseCommand: [python, /rfam/rfam-production/support/infernal2rfam.py]

outputs:
  rfam_file:
    type: File
    outputBinding:
      glob: '*.txt'
      outputEval: ${self[0].basename=inputs.infernal_file.basename + '_rfam.txt'; return self;}

$namespaces:
  edam: http://edamontology.org/

$schemas:
 - http://edamontology.org/EDAM_1.16.owl

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"