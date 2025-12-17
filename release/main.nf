#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
nextflow.preview.output = true

// -- those uncommented below have been completed --
include { FETCH_FAMILIES } from './workflows/fetch_families'

include { GENERATE_3D_SEED } from './workflows/3d_seed'
include { GENERATE_FASTA_FILES } from './workflows/fasta'

include { GENERATE_CLANIN } from './workflows/clanin'
include { GENERATE_CM } from './workflows/cm'
include { GENERATE_FULL_ALIGNMENTS } from './workflows/full_alignments'
include { GENERATE_FULL_REGION } from './workflows/full_region'
include { GENERATE_PDB } from './workflows/pdb'
include { GENERATE_RFAM2GO } from './workflows/rfam2go'
include { GENERATE_SEED } from './workflows/seed'
include { GENERATE_TREE } from './workflows/tree'
include { RUN_VIEW_PROCESS } from './workflows/view_process'
include { LOAD_CM_AND_SEED } from './workflows/load_cm_seed_in_db'

include { clan_competition } from './workflows/clan_competition'
include { update_stockholm_s3 } from './workflows/update_stockholm_s3'

// nice to do, not crucial - needs fixing
// include { prepare_rfam_live } from './workflows/prepare_rfam_live'

// missing params.ena.password, user and hostname
//include { UPLOAD_ENA_MAPPING } from './workflows/ena_mapping'

// requires troubleshooting, not compulsory for release
// include { apicuron } from './workflows/apicuron'

// final step (must be run last)
// include { stage_rfam_live } from './workflows/stage_rfam_live'

// finally
// include { text_search } from './workflows/update_text_search_dev'

workflow {
  main:
    //RUN_VIEW_PROCESS()
    
    FETCH_FAMILIES | set { family_file }
    family_file | splitText | map { it.trim() } | set { families }

    //families | GENERATE_TREE
    //families | GENERATE_FULL_ALIGNMENTS
    // //families | GENERATE_SEED | set { seed_alignments }
    // //seed_alignments | GENERATE_CM
    families | GENERATE_SEED

    GENERATE_SEED.out.seeds | GENERATE_CM
    //GENERATE_SEED.out.seeds | GENERATE_3D_SEED
    //GENERATE_SEED.out.seeds | GENERATE_FASTA_FILES
    // Load CM and seed into database
    LOAD_CM_AND_SEED(
      GENERATE_CM.out.cm_gzip,
      GENERATE_SEED.out.seed_gz
    )



    //GENERATE_PDB | set { pdb }
    //GENERATE_FULL_REGION | set { full_region }
    //GENERATE_RFAM2GO | set { rfam2go }

    // GENERATE_CLANIN | set { clanin }
    
    //clan_competition(Channel.of('start'))
    //update_stockholm_s3(Channel.of(true))


    //pending
    //UPLOAD_ENA_MAPPING()
    //apicuron(Channel.of('start'))
    //prepare_rfam_live(Channel.of('start'))
    
    // stage_rfam_live(Channel.of('start'))

    // Dumps Rfam database XML files, validates them, checks error logs are empty, 
    // creates a release note, and symlinks the data to a dev directory for text search indexing
    // text_search(Channel.of('start'))
}



    //UPLOAD_ENA_MAPPING()
    //RUN_VIEW_PROCESS()
    // clan_competition(Channel.of('start'))
    //FETCH_FAMILIES | set { family_file }

    //GENERATE_CLANIN | set { clanin }
    //GENERATE_FULL_REGION | set { full_region }
    //GENERATE_PDB | set { pdb }
    //GENERATE_RFAM2GO | set { rfam2go }

    //family_file | splitText | map { it.trim() } | set { families }

    //families | GENERATE_TREE
    //families | GENERATE_FULL_ALIGNMENTS
    //families | GENERATE_SEED | set { seed_alignments }

    //seed_alignments | GENERATE_CM
    //seed_alignments | GENERATE_FASTA_FILES
    //seed_alignments | GENERATE_3D_SEED

    // needs fixing:
    // LOAD_CM_AND_SEED(GENERATE_FASTA_FILES.out.all_fasta, GENERATE_CM.out.all_cms, seed_alignments)
  
  //publish:
    // Files published to 'ftp' can be copied directly to the final location as they are already complete
   // GENERATE_SEED.out.seed_gz >> 'ftp'
   // GENERATE_CM.out.cm_gzip >> 'ftp'
   // GENERATE_3D_SEED.out.seeds >> 'ftp'
   // clanin >> 'ftp'
   // full_region >> 'ftp'
   // pdb >> 'ftp'

    // These end up as subdirectories within the ftp directory, but are
    // otherwise complete and ready to be put in place
   // rfam2go >> 'rfam2go'
    //GENERATE_FASTA_FILES.out.fasta >> 'fasta_files'
    //GENERATE_FASTA_FILES.out.all_fasta >> 'fasta_files'
    //GENERATE_FULL_ALIGNMENTS.out.full_alignments >> 'full_alignments'

    // This files require further steps before publishing, generally building
    // a tarball, which seems to be an issue for me in nextflow.
   // GENERATE_TREE.out.seed_trees >> 'seed_tree'
   // GENERATE_CM.out.all_cms >> 'Rfam'
//}

// output {
//   ftp {
//     mode 'copy'
//     path 'ftp'
//   }
// 
//   fasta_files {
//     mode 'copy'
//     path 'ftp/fasta_files'
//   }
// 
//   Rfam {
//     mode 'copy'
//   }
// 
//   rfam2go {
//     mode 'copy'
//     path 'ftp/rfam2go'
//   }
// 
//   seed_tree {
//     mode 'copy'
//     path 'Rfam.seed_tree'
//   }
// }
