#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

include { GENERATE_CM } from './workflows/cm'
include { GENERATE_TREE } from './workflows/tree'
include { GENERATE_SEED } from './workflows/seed'
include { GENERATE_FASTA_FILES } from './workflows/fasta_files'

process FETCH_FAMILIES {
  input:
  val(_flag)

  output:
  file('families')

  """
  mysql -s \
    --host ${params.db.live.host} \
    --port ${params.db.live.port} \
    --user ${params.db.live.user} \
    --password ${params.db.live.password} \
    <<< "select rfam_acc from family" > families
  """
}

process LOAD_CM_AND_SEED {
  input:
  tuple path(cm_file), path(rfam_seed_file)

  """
  load_cm_seed_in_db.py \
    --host ${params.db.live.host} \
    --port ${params.db.live.port} \
    --user ${params.db.live.user} \
    --password ${params.db.live.password} \
    '${rfam_seed_file}' '$cm_file'
  """
}

workflow {
  main:
    FETCH_FAMILIES | set { family_file }

    family_file | splitText | set { families }

    families | GENERATE_RFAM_SEED | set { rfam_seed }
    families | GENERATE_TREE

    GENERATE_CM(family_file, families, rfam_seed) | set { cm }
    GENERATE_FASTA_FILES(families, rfam_seed)

    LOAD_CM_AND_SEED(cm, rfam_seed)
}
