#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { annotated_files } from './annotated_files'
include { mapping_and_updates } from 'pdb_mapping/pdb_mapping'
include { view_process } from './view_process'
include { clan_competition } from './clan_competition'
include { prepare_rfam_live } from './prepare_rfam_live'
include { generate_ftp_files } from './generate_ftp_files'
include { generate_rfam2go } from './generate_rfam2go'
include { stage_rfam_live } from './stage_rfam_live'
include { update_text_search_dev } from './update_text_search_dev'


workflow {
  import_data \
  | ifEmpty('no import') \
  | analyze \
  | ifEmpty('no analysis') \
  | precompute \
  | ifEmpty('no precompute') \
  | genes \
  | ifEmpty('no genes') \
  | export
}