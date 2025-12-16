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
    --host='${params.db.host}' \
    --port='${params.db.port}' \
    --user='${params.db.user}' \
    --password='${params.db.password}' \
    --database='${params.db.name}' \
    --skip-column-names <<< "$sql" > ${acc}.ids
  """
}

process GENERATE_FASTA {
  tag "${acc}"
  memory { params.fasta.largeFamilies.contains(acc) ? '10 GB' : '2 GB' }
  time '12h'

  // For some reason this will fail with no error message. This seems to mean
  // the process was killed by something outside our control. To better deal
  // with that we just retry and hope for the best.
  errorStrategy 'retry'
  maxRetries 5

  input:
  tuple val(acc), path(ids), path(rfam_seed)

  output:
  path("${acc}.fa.gz"), emit: fasta
  //val("${acc}.fa.gz")

  script:
  """
  esl-reformat fasta '${rfam_seed}' > '${acc}.unsorted.fa'
  esl-sfetch -Cf '${params.rfamseq.directory}/${params.rfamseq.combined_fasta}' ${ids} >> '${acc}.unsorted.fa'
  seqkit rmdup '${acc}.unsorted.fa' | seqkit seq --upper-case > '${acc}.fa'
  gzip '${acc}.fa'
  """
}

process COMBINE_FASTA {
  input:
  path(fastas)
  //path("*.fa.gz")  // or rename files in GENERATE_FASTA to match pattern
  //path("family*.fa.gz")

  output:
  path('Rfam.fa.gz'), emit: combined
  //path('Rfam.fa.gz')

  script:
  """
  #find . -name 'family*.fa.gz' | xargs -I {} zcat {} > Rfam.fa.unsorted
  #seqkit rmdup Rfam.fa.unsorted | seqkit seq --upper-case > Rfam.fa
  #gzip Rfam.fa
  if [ -n "${fastas}" ]; then
    cat ${fastas} | gunzip -c > Rfam.fa.unsorted
  else
    touch Rfam.fa.unsorted
  fi
  seqkit rmdup Rfam.fa.unsorted | seqkit seq --upper-case > Rfam.fa
  gzip -f Rfam.fa
  """
}

workflow GENERATE_FASTA_FILES {
  take:
    seeds  // This is coming in as a list/collection

  main:
    query_template = channel.fromPath('sql/ids.sql')
    
    seeds \
      | flatten \
      | map { seed -> 
          def acc = seed.baseName.replaceAll(/\.seed$/, '')
          [acc, seed] 
      } \
      | combine(query_template) \
      | map { acc, seed, template ->
        def binding = [acc: acc]
        def engine = new groovy.text.GStringTemplateEngine()
        def result = engine.createTemplate(template).make(binding)
        [acc, result.toString()]
      } \
      | FETCH_IDS \
      | join(
          seeds \
          | flatten \
          | map { seed -> 
              def acc = seed.baseName.replaceAll(/\.seed$/, '')
              [acc, seed]
          }
      ) \
      | GENERATE_FASTA \
      | collect \
      | COMBINE_FASTA

  emit:
    combined = COMBINE_FASTA.out.combined

}

//workflow GENERATE_FASTA_FILES {
//  take:
//    accessions
//    seed_alignments
//  emit:
//    fasta = fasta
//    all_fasta = all_fasta
//  main:
//    query_template = channel.fromPath('sql/ids.sql')
//
//    accessions \
//    | combine(query_template) \
//    | map { acc, template ->
//      def binding = [acc: acc]
//      def engine = new groovy.text.GStringTemplateEngine()
//      def result = engine.createTemplate(template).make(binding)
//      [acc, result.toString()]
//    } \
//    | FETCH_IDS \
//    | join(seed_alignments) \
//    | GENERATE_FASTA \
//    | set { fasta }
//
//    fasta | collect | COMBINE_FASTA | set { all_fasta }
//}
