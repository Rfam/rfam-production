#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  sequences:
    type: File
  covariance_model_database:
    type: File

outputs:
  chunked_sequence:
    type: File[]
    outputSource: chunk/sequence_chunks
  cmsearch_matches:
    type: File[]
    outputSource: cmsearch/matches

steps:
  chunk:
    in:
      sequences: sequences
      num_chunks:
        valueFrom: $(Math.ceil(inputs.sequences.size / 2000000))
    out:
      - sequence_chunks
    run: ../tools/esl-chunk.cwl

  seqstat:
    in:
      sequences: sequences
    out:
      - nt_count
    run: ../tools/esl-seqstat.cwl

  cmsearch:
    scatter: query_sequences
    in:
      covariance_model_database: covariance_model_database
      query_sequences: chunk/sequence_chunks
      cut_ga:
        default: True
      acc:
        default: True
      notextw:
        default: True
      nohmmonly:
        default: True
      search_space_size:
        source: seqstat/nt_count
        # size * 2 as both strands are being searched and divided by 1M as the size needs to be in Mb
        valueFrom: $(Math.ceil(self * 2 /Math.pow(10,6)))
    out:
      - matches
    run: ../tools/cmsearch.cwl

$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/docs/schema_org_rdfa.html

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
