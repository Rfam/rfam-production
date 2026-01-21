process BUILD {
  input:
  path(query)

  output:
  path("Rfam_${params.release_info.date}_rfam2embl_crossrefs.txt")

  script:
  """
  mysql -s \
    --host=${params.db.host} \
    --port=${params.db.port} \
    --user=${params.db.user} \
    --database=${params.db.name} \
    --password=${params.db.password} \
    < ${query} > Rfam_${params.release_info.date}_rfam2embl_crossrefs.txt
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
  # Upload file using curl
  curl -T ${crossref_file} \
       -u ${params.ena.user}:${params.ena.password} \
       ftp://${params.ena.hostname}/xref/
  
  # Verify upload
  echo "Verifying upload..."
  curl -u ${params.ena.user}:${params.ena.password} \
       --head ftp://${params.ena.hostname}/xref/${crossref_file}
  
  if [ \$? -eq 0 ]; then
    echo "Upload successful!"
  else
    echo "Upload verification failed"
    exit 1
  fi
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
