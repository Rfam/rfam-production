process fetch_genomes {
  input:
  val(_flag)

  output:
  path('genomes/*.jsonl')

  """
  rfam_export --rfam-version ${params.rfam_version} --rfamseq-version ${params.rfamseq_version} list-genomes all-genomes.txt
  mkdir genomes
  shuf all-genomes.txt > genomes-shuf.txt
  split \
    -n l/${params.genome_chunks} \
    --additional-suffix='.txt' \
    genomes-shuf.txt genomes/
  """
}

process extract_hits {
  tag { "$genomes.baseName" }
  maxForks 10
  memory 10.GB
  publishDir "${params.text_search}/full_region", mode: "copy"
  errorStrategy 'ignore'

  input:
  path(genomes)

  output:
  path("full_regions/*.xml")

  """
  rfam_export --rfam-version ${params.rfam_version} --rfamseq-version ${params.rfamseq_version} full-regions $genomes full_regions
  """
}

workflow export {
  take:
    start
  main:
    start | fetch_genomes | flatten | extract_hits | collect | map { "export done" } | set { done }
  emit:
    done
}

workflow {
  export(Channel.of('start'))
}
