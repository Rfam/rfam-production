process fetch_ncbi_locations {
  time '1h'

  input:
  path(urls)

  output:
  path('ncbi.db')

  """
  mkdir summaries
  wget --quiet --input-file ncbi-urls.txt -P summaries
  rfamseq ncbi parse-assembly-summary summaries/*.txt ncbi.db
  """
}

process fetch_viral_additions {
  time '1h'

  output:
  path("virus-proteomes.txt")

  """
  wget "${params.additional.viruses.source}" -O pir-virus.txt
  grep '^>' pir-virus.txt | awk '{ print \$1 }' | tr -d '>' > virus-proteomes.txt
  """
}

process download_uniprot_summary {
  time '1h'
  publishDir 'genomes/uniprot', mode: 'copy'

  output:
  path('summary.json')

  """
  curl '$params.proteome_json' > summary.json
  """
}

process process_uniprot_proteomes {
  time '1h'

  input:
  tuple path(reference), path("additional.txt")

  output:
  path("unique.jsonl")

  """
  rfamseq uniprot parse-proteomes $reference reference.jsonl
  rfamseq uniprot fetch-proteomes additional.txt additional.jsonl
  rfamseq uniprot deduplicate reference.jsonl additional.jsonl unique.jsonl
  """
}

process chunk_genomes {
  input:
  path("genomes.jsonl")

  output:
  path("parts/*.jsonl")

  """
  mkdir parts
  shuf genomes.jsonl > genomes-shuf.jsonl
  split \
    -n l/${params.proteome_chunks} \
    --additional-suffix='.jsonl' \
    genomes-shuf.jsonl parts/
  """
}

process download {
  tag { "$proteome_file.baseName" }
  queue 'datamover'
  maxForks params.download_forks

  publishDir "genomes/fasta/${proteome_file.baseName}", mode: "copy", pattern: "UP*.fa"
  publishDir "genomes/metadata/${proteome_file.baseName}", mode: "copy", pattern: "UP*.jsonl"
  publishDir "genomes/failures/", mode: "copy", pattern: "*-failed.jsonl"

  memory { 6.GB * task.attempt }
  errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'finish' }
  maxRetries 4

  input:
  tuple path(proteome_file), path(info)

  output:
  tuple val("${proteome_file.baseName}"), path("UP*.fa"), emit: sequences
  tuple val("${proteome_file.baseName}"), path("UP*.jsonl"), emit: metadata
  tuple val("${proteome_file.baseName}"), path("${proteome_file.baseName}-failed.jsonl"), emit: failed

  """
  rfamseq uniprot download-genomes \
    --failed-filename "${proteome_file.baseName}-failed.jsonl"  \
    ${params.version} $info $proteome_file .
  """
}

process validate_chunk {
  tag { "$short_name" }
  time '12h'
  memory { 10.GB * task.attempt }
  errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'finish' }
  maxRetries 3

  input:
  tuple val(short_name), path('genomes*.fa')

  output:
  path("${short_name}.fa")

  """
  set -euo pipefail

  find . -name 'genomes*.fa' | xargs -I {} seqkit rmdup -s {} | seqkit rmdup > unique.fa
  seqkit shuffle --two-pass unique.fa > ${short_name}.fa
  esl-sfetch --index ${short_name}.fa
  """
}

process merge_chunks {
  time '2d'
  publishDir 'genomes/rfamseq', mode: 'copy'
  memory 10.GB

  input:
  path('genomes*.fa')

  output:
  tuple path("rfamseq_${params.version}.fa"), path("rfamseq_${params.version}.fa.ssi"), emit: sequences
  path("rfamseq_${params.version}.seqstat"), emit: seqstat

  """
  set -euo pipefail

  find . -name 'genomes*.fa' | xargs cat > rfamseq_${params.version}.fa
  esl-seqstat -a rfamseq_${params.version}.fa > rfamseq_${params.version}.seqstat
  esl-sfetch --index rfamseq_${params.version}.fa
  """
}

process split_seqstat {
  memory 10.GB

  input:
  path(seqstat)

  output:
  path("rfamseq-chunks/*"), emit: rfamseq_chunks
  path("rev-chunks/*"), emit: rev_chunks

  """
  shuf ${seqstat} > shuffled.seqstat
  rfamseq seqstat chunk-ids shuffled.seqstat ${params.rfam_seq.chunk_size} rfamseq-chunks/

  rfamseq seqstat take-fraction shuffled.seqstat ${params.rfam_seq.rev_fraction} rev-selected.seqstat
  rfamseq seqstat chunk-ids rev-selected.seqstat ${params.rfam_seq.chunk_size} rev-chunks/
  """
}

process chunk_rfamseq {
  tag { "${ids.baseName}" }
  publishDir 'genomes/rfamseq', mode: 'copy'

  input:
  tuple path(merged), path(ssi), path(ids)

  output:
  path("rfamseq_${ids.baseName}.fa.gz")

  script:
  """
  esl-sfetch -f $merged $ids > rfamseq_${ids.baseName}.fa
  gzip rfamseq_${ids.baseName}.fa
  """
}

process chunk_rev_rfamseq {
  tag { "${ids.baseName}" }
  publishDir 'genomes/rfamseq', mode: 'copy'

  input:
  tuple path(merged), path(ssi), path(ids)

  output:
  path("rev-rfamseq_${ids.baseName}.fa.gz")

  script:
  """
  esl-sfetch -f $merged $ids > rfamseq_${ids.baseName}.fa
  esl-shuffle -r rfamseq_${ids.baseName}.fa > rev-rfamseq_${ids.baseName}.fa
  gzip rev-rfamseq_${ids.baseName}.fa
  """
}

workflow genome_download {
  main:
    Channel.fromPath('ncbi-urls.txt') | fetch_ncbi_locations | set { ncbi_info }

    fetch_viral_additions | set { additional_uniprot }

    download_uniprot_summary \
    | combine(additional_uniprot) \

    | process_uniprot_proteomes \
    | chunk_genomes \
    | flatten \
    | combine(ncbi_info) \
    | download

    download.out.sequences \
    | validate_chunk \
    | collect \
    | merge_chunks

    merge_chunks.out.seqstat | split_seqstat
    split_seqstat.out.rfam_chunks | flatten | set { rfam_chunks }
    split_seqstat.out.rev_chunks | flatten | set { rev_chunks }

    merge_chunks.out.sequences \
    | combine(rfam_chunks) \
    | chunk_rfamseq

    merge_chunks.out.sequences \
    | combine(rev_chunks) \
    | chunk_rev_rfamseq
}

workflow {
  genome_download()
}
