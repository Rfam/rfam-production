#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

include { GENERATE_CM } from './workflows/cm'
include { GENERATE_TREE } from './workflows/tree'
include { GENERATE_SEED } from './workflows/seed'
include { GENERATE_FASTA_FILES } from './workflows/fasta'
include { FETCH_FAMILIES } from './workflows/fetch_families'
include { GENERATE_FULL_ALIGNMENTS } from './workflows/full_alignments'
include { GENERATE_FULL_REGION } from './workflows/full_region'
include { GENERATE_CLANIN } from './workflows/clanin'
include { GENERATE_PDB } from './workflows/pdb'
// include { LOAD_CM_AND_SEED } from './workflows/load_cm_seed_in_db'

// TODO: Setup publishing properly

// Easy:
// TODO: Rfam.clanin
// TODO: Rfam.full_region
// TODO: Rfam.pdb
// TODO: Rfam2go

// Harder:
// TODO: Upload ENA mapping
// TODO: 3D seed - Need to get list of families with 3D

workflow {
  main:
    FETCH_FAMILIES | set { family_file }

    GENERATE_CLANIN | set { clanin }
    GENERATE_FULL_REGION | set { full_region }
    GENERATE_PDB | set { pdb }

    family_file | splitText | map { it.trim() } | set { families }

    families | GENERATE_SEED | set { seed_alignments }
    families | GENERATE_TREE

    seed_alignments | GENERATE_CM
    seed_alignments | GENERATE_FASTA_FILES
    seed_alignments | GENERATE_FULL_ALIGNMENTS

    // LOAD_CM_AND_SEED(GENERATE_FASTA_FILES.out.all_fasta, GENERATE_CM.out.all_cms, seed_alignments)
  publish:
    // Files published to 'ftp' can be copied directly to the final location as they are already complete
    GENERATE_SEED.out.seed_gz >> 'ftp'
    GENERATE_CM.out.cm_gzip >> 'ftp'
    clanin >> 'ftp'
    full_region >> 'ftp'
    pdb >> 'ftp'
    rfam2go >> 'ftp/rfam2go'
    GENERATE_FASTA_FILES.out.fasta >> 'ftp/fasta_files'
    GENERATE_FASTA_FILES.out.all_fasta >> 'ftp/fasta_files'
    GENERATE_FULL_ALIGNMENTS.out.full_alignments >> 'ftp/full_alignments'

    /// This files require further steps before publishing, generally building
    // a tarball, which seems to be an issue for me in nextflow.
    GENERATE_TREE.out.seed_trees >> 'Rfam.seed_tree'
    GENERATE_CM.out.all_cms >> 'Rfam'
}

output {
  mode 'copy'
}
