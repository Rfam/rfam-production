// It is a slightly odd design to split getting ids from building the files. But
// for some reason sometimes the cluster kills the combined process. So this way
// there is less to redo. Plus we don't have to limit the number of jobs for
// GENERATE_FASTA.
process FETCH_IDS {
  tag "${acc}"
  maxForks 100

  input:
  tuple val(acc), val(sql)

  output:
  tuple val(acc), path("${acc}.ids")

  script:
  """
  mysql \
    --host='${params.db.live.host}' \
    --port='${params.db.live.port}' \
    --user='${params.db.live.user}' \
    -p'${params.db.live.password}' \
    --database='${params.db.live.database}' \
    --skip-column-names <<< "$sql" > ${acc}.ids
  """
}

process GENERATE_FASTA {
  tag "${acc}"
  memory { params.fasta.largeFamilies.contains(acc) ? 10.GB : 2.GB }
  time 12.hours

  // For some reason this will fail with no error message. This seems to mean
  // the process was killed by something outside our control. To better deal
  // with that we just retry and hope for the best.
  errorStrategy 'retry'
  maxRetries 5

  input:
  tuple val(acc), path(ids), path(rfam_seed)

  output:
  val("${acc}.fa.gz")

  """
  esl-reformat fasta '${rfam_seed}'                                                     > '${acc}.unsorted.fa'
  esl-sfetch -Cf '${params.rfamseq.directory}/${params.rfamseq.combined_fasta}' ${ids} >> '${acc}.unsorted.fa'
  seqkit rmdup '${acc}.unsorted.fa' | seqkit seq --upper-case > '${acc}.fa'
  gzip '${acc}.fa'
  """
}

process COMBINE_FASTA {
  input:
  path("family*.fa.gz")

  output:
  path('Rfam.fa.gz')

  """
  find . -name 'family*.fa.gz' | xargs -I {} zcat {} > Rfam.fa.unsorted
  seqkit rmdup Rfam.fa.unsorted | seqkit seq --upper-case > Rfam.fa
  gzip Rfam.fa
  """
}

workflow GENERATE_FASTA_FILES {
  take:
    accessions
    seed_alignments
  emit:
    fasta
    all_fasta
  main:
    query_template = channel.fromPath('sql/ids.sql')

    accessions \
    | combine(query_template) \
    | map { acc, template ->
      def binding = [acc: acc]
      def engine = new groovy.text.GStringTemplateEngine()
      def result = engine.createTemplate(template).make(binding)
      [acc, result.toString()]
    } \
    | FETCH_IDS \
    | join(seed_alignments) \
    | GENERATE_FASTA \
    | set { fasta }

    fasta | collect | COMBINE_FASTA | set { all_fasta }
}
