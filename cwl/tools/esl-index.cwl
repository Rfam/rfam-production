cwlVersion: v1.0
class: CommandLineTool
label: index a sequence file for use by esl-sfetch
doc: "https://github.com/EddyRivasLab/easel"

requirements:
  ResourceRequirement:
    coresMax: 1
    ramMin: 1024  # just a default, could be lowered
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: [ $(inputs.sequences) ]
hints:
  DockerRequirement:
    dockerPull: rfam/rfam-production:cwl-infernal
  SoftwareRequirement:
    packages:
      easel: {}

inputs:
  sequences:
    type: File
    inputBinding:
      position: 1
      valueFrom: $(self.basename)
    format: edam:format_1929  # FASTA

baseCommand: [ esl-sfetch, --index ]

outputs:
  sequence_index_file:
    type: File
    outputBinding:
      glob: $(inputs.sequences.basename + '.ssi')

$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
