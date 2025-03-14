process BUILD {
  input:
  path(query)

  output:
  path("Rfam_${params.date}_rfam2embl_crossrefs.txt")

  script:
  """
  mysql -s \
    --host=${params.db.live.host} \
    --port=${params.db.live.port} \
    --user=${params.db.live.user} \
    --database=${params.db.live.database} \
    --password=${params.db.live.password} \
    < $query > Rfam_${params.date}_rfam2embl_crossrefs.txt
  """
}

process UPLOAD {
  when params.ena_mapping.upload

  input:
  path("Rfam_${params.date}_rfam2embl_crossrefs.txt")

  """
  ftp -n <<EOF
  open ${params.ena.hostname}
  user ${params.ena.user} ${params.ena.password}
  cwd /xref
  stor Rfam_${params.date}_rfam2embl_crossrefs.txt
  EOF
  """
}

workflow UPLOAD_ENA_MAPPING {
  emit:
    ena_mapping
  main:
    channel.fromPath("${moduleDir}/sql/ena_mapping.sql") \
    | BUILD \
    | set { ena_mapping }

    ena_mapping | UPLOAD
}
