cwlVersion: v1.0
class: CommandLineTool
label: search sequence(s) against a covariance model database
doc: "http://eddylab.org/infernal/Userguide.pdf"

requirements:
  ResourceRequirement:
    ramMin: 4096
    coresMin: 4
hints:
  DockerRequirement:
    dockerPull: rfam/rfam-production:easel
  SoftwareRequirement:
    packages:
      infernal:
        specs: [ "https://identifiers.org/rrid/RRID:SCR_011809" ]
        version: [ "1.1.2" ]
  #- $import: infernal-docker.yml

inputs:
  covariance_model_database:
    type: File
    inputBinding:
      position: 1

  query_sequences:
    type: File
    streamable: true
    inputBinding:
      position: 1
    format:
      - edam:format_1929  # FASTA
      # - edam:format_1927  # EMBL
      # - edam:format_1936  # Genbank entry format
      # - edam:format_1961  # Stockholm
      # - edam:format_3281  # A2M
      # - edam:format_1982  # Clustal
      # - edam:format_1997  # PHYLIP
      # ddbj ?
      # pfam ?
      # afa ?

  search_space_size:
    type: int
    label: search space size in *Mb* to <x> for E-value calculations
    inputBinding:
      prefix: -Z

  omit_alignment_section:
    label: Omit the alignment section from the main output.
    type: boolean?
    inputBinding:
      prefix: --noali
    doc: This can greatly reduce the output volume.

  only_hmm:
    label: Only use the filter profile HMM for searches, do not use the CM
    type: boolean?
    default: false
    inputBinding:
      prefix: --hmmonly
    doc: |
      Only filter stages F1 through F3 will be executed, using strict P-value
      thresholds (0.02 for F1, 0.001 for F2 and 0.00001 for F3). Additionally
      a bias composition filter is used after the F1 stage (with P=0.02
      survival threshold).  Any hit that survives all stages and has an HMM
      E-value or bit score above the reporting threshold will be output.

  cut_ga:
    label: use CM's GA gathering cutoffs as reporting thresholds
    type: boolean?
    default: false
    inputBinding:
      prefix: --cut_ga

  rfam:
    label: set heuristic filters at Rfam-level (fast)
    type: boolean?
    default: false
    inputBinding:
      prefix: --rfam

  nohmmonly:
    label: Disable hmm-only thresholds
    type: boolean?
    default: false
    inputBinding:
      prefix: --nohmmonly

  notextw:
    label: Unlimit the length of each line in the main output
    type: boolean?
    default: false
    inputBinding:
      prefix: --notextw
  acc:
    label: Output RFam accession rather than Rfam ID
    type: boolean?
    default: false
    inputBinding:
      prefix: --acc

baseCommand: cmsearch

arguments:
  - valueFrom: matches.tbl
    prefix: --tblout
  - valueFrom: $(runtime.cores)
    prefix: --cpu

outputs:
  matches:
    label: target hits table, format 2
    doc: "http://eddylab.org/infernal/Userguide.pdf#page=60"
    type: File
    outputBinding:
      glob: matches.tbl
      outputEval: ${self[0].basename=inputs.query_sequences.basename + '_matches.tbl'; return self;}

$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
