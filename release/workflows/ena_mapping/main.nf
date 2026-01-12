process BUILD {
  input:
  path(query)

  output:
  path("Rfam_${params.date}_rfam2embl_crossrefs.txt")

  script:
  """
  mysql -s \
    --host=${params.db.host} \
    --port=${params.db.port} \
    --user=${params.db.user} \
    --database=${params.db.name} \
    --password=${params.db.password} \
    < ${query} > Rfam_${params.date}_rfam2embl_crossrefs.txt
  """
}

process UPLOAD {
  when:
  params.ena_mapping.upload

  input:
  path(crossref_file)
  //path("Rfam_${params.date}_rfam2embl_crossrefs.txt")

  output:
  path(crossref_file), optional: true

  script:
  """
  ftp -n <<'EOF'
  open ${params.ena.hostname}
  user ${params.ena.user} ${params.ena.password}
  cwd /xref
  stor ${crossref_file}
  EOF
  """
}

workflow UPLOAD_ENA_MAPPING {
  main:
    Channel.fromPath("${moduleDir}/sql/ena_mapping.sql") \
    | BUILD \
    | set { ena_mapping }

    ena_mapping | UPLOAD

  emit:
    ena_mapping
}
