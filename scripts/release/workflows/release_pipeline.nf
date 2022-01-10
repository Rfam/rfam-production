#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { generate_annotated_files } from './annotated_files'
include { mapping_and_updates } from 'pdb_mapping/pdb_mapping'
include { view_process } from './view_process'
include { clan_competition } from './clan_competition'
include { prepare_rfam_live } from './prepare_rfam_live'
include { generate_ftp_files } from './generate_ftp_files'
include { generate_rfam2go } from './rfam2go'
include { stage_rfam_live } from './stage_rfam_live'
include { text_search } from './update_text_search_dev'


workflow {
  generate_annotated_files \
  | ifEmpty('no annotated files') \
  | mapping_and_updates \
  | ifEmpty('no pdb_mapping') \
  | view_process \
  | ifEmpty('no view processes') \
  | clan_competition \
  | ifEmpty('no clan compete') \
  | prepare_rfam_live
  | ifEmpty('no prepare rfam_live') \
  | generate_ftp_files \
  | ifEmpty('no ftp files') \
  | generate_rfam2go \
  | ifEmpty('no rfam2go') \
  | stage_rfam_live \
  | ifEmpty('no stage rfam_live') \
  | text_search
}
