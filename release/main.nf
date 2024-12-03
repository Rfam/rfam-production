#!/usr/bin/env nextflow

include { cms } from './workflows/cms'
include { complete_clans } from './workflows/complete_clans'
include { database_export } from './workflows/database_export'
include { fasta_export } from './workflows/fasta_files'
include { rfam2go } from './workflows/rfam2go'
include { seeds } from './workflows/seeds'
include { templates } from './workflows/templates'
include { trees } from './workflows/trees'
include { view_processes } from './workflows/view_processes'

process fetch_families {
  output:
  path('families.tsv')

  """
  mysql \
    --host "${params.db.host}" \
    --port "${params.db.port}" \
    --user "${params.db.user}" \
    "-p${params.db.password}" \
    --database "${params.database.name}" <<< 'select rfam_acc from family order by rfam_acc' \
    --skip-column-names > families.txt
  """
}

process load_to_database {
  input:
  path(families)
  path(cm)
  path(seed)

  output:
  val('done')

  """
  gzip -c $cm > ${accession}.cm.gz
  gzip -c ${seed} > ${accession}.seed.gz
  load_cm_seed_in_db.py ${accession}.seed.gz ${accession}.cm.gz
  """
}

workflow {
  main:
    Channel.fromPath(params.families_with_3d) | set { fam_3d }
    fetch_families | set { accession_file }
    accession_file | splitCsv(sep:"\t") |  set { families }

    seeds(families, fam_3d)
    cms(accession_file, families)
    trees(families)

    load_to_database(accession_file, seeds.out.rfam_seed, cms.out.gzip_cm) \
    | view_processes \
    | complete_clans \
    | prepare_database \
    | (database_export & templates & rfam2go) \
    | mix \
    | collect \
    | map { "done" }
  publish:
    fasta_export.out.sequences >> 'fasta_files'
    fasta_export.out.merged >> 'fasta_files'
}

output {
  directory "$launchDir/ftp"
  mode "copy"
}
