#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENERATE_CM } from './workflows/cm'
include { GENERATE_TREE } from './workflows/tree'
include { GENERATE_SEED } from './workflows/seed'
include { GENERATE_FASTA_FILES } from './workflows/fasta_files'
include { FETCH_FAMILIES } from './workflows/fetch_families'
// include { LOAD_CM_AND_SEED } from './workflows/load_cm_seed_in_db'

workflow {
  main:
    FETCH_FAMILIES | set { family_file }

    family_file | splitText | set { families }

    families | GENERATE_RFAM_SEED | set { rfam_seed }
    families | GENERATE_TREE

    GENERATE_CM(family_file, families, rfam_seed) | set { cm }
    GENERATE_FASTA_FILES(families, rfam_seed)

    // LOAD_CM_AND_SEED(family_file, cm, rfam_seed)
}
