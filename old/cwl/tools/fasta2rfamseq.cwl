#/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Tool to generate taxonomy entries based on a list of tax ids

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
    type: stdout

stdout: test.rfamseq
