#/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Tool to generate taxonomy entries based on a list of tax ids

# use javascript for input name parsing
requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMax: 1
    ramMin: 1024  # just a default, could be lowered

hints:
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
    type: File
    outputBinding:
      glob: ${self[0].basename=inputs.fasta_file.basename + '.fa'; return self;}.rfamseq
