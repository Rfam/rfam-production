process full_region {
  publishDir "$params.release_ftp", mode: "copy"

  input:
  path(query)

  output:
  path('Rfam.full_region.gz')

  """
  mysql \
    --host "${params.db.host}" \
    --port "${params.db.port}" \
    --user "${params.db.user}" \
    "-p${params.db.password}" \
    --database "${params.database.name}" < $query > Rfam.full_region
  gzip Rfam.full_region
  """
}

process pdb_file {
  publishDir "$params.release_ftp", mode: "copy"

  input:
  path(query)

  output:
  path('Rfam.pdb.gz')

  """
  mysql \
    --host "${params.db.host}" \
    --port "${params.db.port}" \
    --user "${params.db.user}" \
    "-p${params.db.password}" \
    --database "${params.database.name}" < $query > Rfam.pdb
  gzip Rfam.pdb
  """
}

process clanin {
  publishDir "$params.release_ftp", mode: "copy"

  output:
  path('Rfam.clanin')

  """
  clanin_file_generator.py Rfam.clanin
  """
}

process upload_ena_mapping {
  input:
  path(query)

  output:
  path('ena.txt')

  """
  mysql \
    --host "${params.db.host}" \
    --port "${params.db.port}" \
    --user "${params.db.user}" \
    "-p${params.db.password}" \
    --database "${params.database.name}" < $query > ena.txt

  upload_ena.py ena.txt
  """
}

workflow database_export {
  take:
    start
  main:
    Channel.fromPath('sql/ftp_rfam_full_reqion.sql') \
    | full_region \
    | set { full_region_gz }

    Channel.fromPath('sql/ftp_rfam_pdb.sql') \
    | pdb_file \
    | set { pdb_gz }

    Channel.fromPath('sql/ena_mapping.sql') \
    | upload_ena
    | set { ena }

    clanin | set { clanin_file }

    ena \
    | mix(clanin_file, full_region_gz) \
    | collect \
    | map { "database export done" } \
    | ifEmpty("No export") \
    | set { done }
  emit:
    done
}
