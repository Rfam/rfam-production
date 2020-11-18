#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  sequences:
    type: File
  covariance_model_database:
    type: File

outputs:
  genome_scanner_matches:
    type: File[]
    outputSource: genome_scanner/cmsearch_matches
  rfam_files:
    type: File[]
    outputSource: infernal2rfam/rfam_file
  merged_infernal:
    type: File
    outputSource: merge_infernal_files/result
  merged_gff:
    type: File
    outputSource: merge_gff/result

steps:
  genome_scanner:
    in:
      sequences: sequences
      covariance_model_database: covariance_model_database
    out:
      - cmsearch_matches
    run: ./genome_scanner.cwl

  infernal2rfam:
    scatter: infernal_file
    in:
      infernal_file: genome_scanner/cmsearch_matches
      infernal_tbl:
        default: True
    out:
      - rfam_file
    run: ../tools/infernal2rfam.cwl

  merge_infernal_files:
    in:
      files: infernal2rfam/rfam_file
      outfile_name:
        default: 'merged.txt'
    out:
      - result
    run: ../tools/concat_files.cwl

  infernal2gff:
    scatter: infernal_file
    in:
      infernal_file: genome_scanner/cmsearch_matches
      all:
        default: True
    out:
      - gff
    run: ../tools/infernal2gff/infernal2gff.cwl

  merge_gff:
    in:
      files: infernal2gff/gff
      outfile_name:
        default: 'merged.gff'
    out:
      - result
    run: ../tools/concat_files.cwl

$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
