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
    dockerPull: miguelboland/easel:latest
  SoftwareRequirement:
    packages:
      easel: {}
        # specs: [ https://identifiers.org/rrid/RRID:TBD ]
        # version: [ "???" ]

inputs:
  no_sum:
    label: report per-sequence info line, not just a summary
    type: boolean?
    default: false
    inputBinding:
      prefix: -a
  residue:
    label: count and report residue composition
    type: boolean?
    default: false
    inputBinding:
      prefix: -c
  in_fmt:
    label: specify that input file is in format <s>
    type: string?
    inputBinding:
      prefix: --informat
  is_rna:
    label: specify that input contains RNA sequence
    type: boolean?
    default: false
    inputBinding:
      prefix: --rna
  is_dna:
    label: specify that input contains DNA sequence
    type: boolean?
    default: false
    inputBinding:
      prefix: --dna
  amino:
    label: specify that input contains protein sequence
    type: boolean?
    default: false
    inputBinding:
      prefix: --amino
  comptbl:
    label: 'alternative output: a table of residue compositions per seq'
    type: boolean?
    default: false
    inputBinding:
      prefix: --comptbl
  sequences:
    type: File
    inputBinding:
      valueFrom: $(self.basename)
    format: edam:format_1929  # FASTA

stdout: seqstat.txt

baseCommand: [esl-seqstat]

outputs:
  nt_count:
    type: int
    outputBinding:
      glob: seqstat.txt
      loadContents: true
      outputEval: |
        ${
          var re = /Total.+:\s+(\d+)/;
          var val = re.exec(self[0].contents)[1];
          return parseInt(val);
        }




$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
