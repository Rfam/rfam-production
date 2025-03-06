#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

// include { GENERATE_CM } from './workflows/cm'
include { GENERATE_TREE } from './workflows/tree'
include { GENERATE_SEED } from './workflows/seed'
include { GENERATE_FASTA_FILES } from './workflows/fasta_files'
include { FETCH_FAMILIES } from './workflows/fetch_families'
// include { GENERATE_ALIGNMENTS } from './workflows/generate_alignments'
// include { LOAD_CM_AND_SEED } from './workflows/load_cm_seed_in_db'

workflow {
  main:
    FETCH_FAMILIES | set { family_file }

    family_file | splitText | map { it.trim() } | set { families }

    families | GENERATE_SEED | set { seed_alignments }
    families | GENERATE_TREE

    seed_alignments | GENERATE_CM | set { cm }
    seed_alignments | GENERATE_FASTA_FILES | set { fasta }
    // GENERATE_ALIGNMENTS(fasta, rfam_seed)

    // LOAD_CM_AND_SEED(family_file, cm, rfam_seed)
}
