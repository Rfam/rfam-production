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
    dockerPull: miguelboland/easel:latest
  SoftwareRequirement:
    packages:
      easel: {}

inputs:
  sequences:
    type: File
    format: edam:format_1929  # FASTA
    inputBinding:
      position: 3
      valueFrom: $(self.basename)
  num_chunks:
    type: int
    inputBinding:
      position: 4

baseCommand: [ esl-ssplit.pl, '-n', '-r']

successCodes: [0, 2] # esl-ssplit will fail if chunk number == 1, which is fine for us.

outputs:
  sequence_chunks:
    type: File[]
    format: edam:format_1929  # FASTA
    outputBinding:
      glob: $(inputs.sequences.basename + '*')
      # If number of files is larger than 1, do not include the input file which was matched by glob
      outputEval: |
        ${
          if (self.length === 1){
            return self;
          } else {
            return self.filter(function(el){return el.basename !== inputs.sequences.basename});
          }
        }

$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
