cwlVersion: v1.0
class: CommandLineTool
label: Chunk a sequence file using esl-ssplit
doc: "https://github.com/EddyRivasLab/bio-easel"

requirements:
  ResourceRequirement:
    coresMax: 1
    ramMin: 1024  # just a default, could be lowered
  InitialWorkDirRequirement:
    listing: [ $(inputs.sequences) ]
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: ikalvari/rfam-cloud:dev
  SoftwareRequirement:
    packages:
      easel: {}
        # specs: [ https://identifiers.org/rrid/RRID:TBD ]
        # version: [ "???" ]

inputs:
  sequences:
    type: File
#    format: edam:format_1929  # FASTA
    inputBinding:
      position: 3
      valueFrom: $(self.basename)
  num_chunks:
    type: int
    inputBinding:
      position: 4

baseCommand: [ esl-ssplit.pl, '-n', '-r']

outputs:
  chunks:
    type: File[]
    format: edam:format_1929  # FASTA
    outputBinding:
      glob: $(inputs.sequences.basename + '.*')

$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
