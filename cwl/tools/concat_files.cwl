#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}  # to propagate the file format
  ResourceRequirement:
    coresMax: 1
    ramMin: 100  # just a default, could be lowered

inputs:
  files:
    type: File[]
    streamable: true
    inputBinding:
      position: 1
  outfile_name:
    type: string?
    default: result

baseCommand: cat

stdout: result  # to aid cwltool's cache feature

outputs:
  result:
    type: File
    outputBinding:
      glob: result
      outputEval: |
        ${ self[0].format = inputs.files[0].format;
           self[0].basename = inputs.outfile_name;
           return self;
         }


s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"