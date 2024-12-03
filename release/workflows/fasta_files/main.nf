process export_family {
  tag { "$family" }
  maxForks 10

  input:
  tuple val(family), path(query)

  output:
  path("${family}.fa.gz")

  """
  mysql \
    --host "${params.db.host}" \
    --port "${params.db.port}" \
    --user "${params.db.user}" \
    "-p${params.db.password}" \
    --database "${params.database.name}"
    --skip-column-names \
    -e "set @family=${family}; source $query" < $query > ids
  esl-sfetch -Cf ids ${params.rfamseq.complete} > ${family}.fa
  gzip ${family}.fa
  """
}

process merge_fasta {
  input:
  path(families)

  output:
  path("Rfam.fa.gz")

  """
  zcat $families > Rfam.fa
  gzip Rfam.fa
  """
}

workflow fasta_export {
  take:
    families
  emit:
    sequences
    merged
  main:
    Channel.fromPath('sql/sequence-ids.sql') | set { query }

    families | join(query) | export_family | set { sequences }
    sequences | collect | merge_fasta | set { merged }
}
