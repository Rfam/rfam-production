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
  rfamseq uniprot fetch-proteomes additional.txt addtional.jsonl
  rfamseq uniprot deduplicate reference.jsonl additional.jsonl unique.jsonl
  """
}

process chunk_genomes {
  input:
  path("genomes.json")

  output:
  path("parts/*.jsonl")

  """
  mkdir parts
  find . -name 'proteomes*.jsonl' > summary.jsonl
  shuf summary.jsonl > summary-shuf.jsonl
  split \
    -n l/${params.proteome_chunks} \
    --additional-suffix='.jsonl' \
    summary-shuf.jsonl parts/
  """
}

process download {
  tag { "$proteome_file.baseName" }
  queue 'datamover'
  maxForks params.download_forks

  publishDir "genomes/fasta/${proteome_file.baseName}", mode: "copy", pattern: "UP*.fa"
  publishDir "genomes/metadata/${proteome_file.baseName}", mode: "copy", pattern: "UP*.jsonl"
  publishDir "genomes/failures/${proteome_file.baseName}", mode: "copy", pattern: "*failed.jsonl"

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

  input:
  path('genomes*.fa')

  output:
  tuple path("rfamseq_${params.version}.fa"), path("rfamseq_${params.version}.fa.ssi")

  """
  set -euo pipefail

  find . -name 'genomes*.fa' | xargs cat > rfamseq_${params.version}.fa
  esl-seqstat -a rfamseq_${params.version}.fa > rfamseq_${params.version}.seqstat
  esl-sfetch --index rfamseq_${params.version}.fa
  """
}

process build_rfamseq {
  container ''
  publishDir 'genomes/rfamseq', mode: 'copy'

  input:
  tuple path(merged), path(ssi)

  output:
  path('*.fa.gz')

  script:
  name = "r${params.rfam_seq.main_chunks}_rfamseq${params.version}"
  """
  mkdir rfamseq
  pushd rfamseq
  /homes/rfamprod/Bio-Easel/scripts/esl-ssplit.pl -v --oroot ${name}.fa -n -r -z ../${merged} ${params.rfam_seq.main_chunks}
  for((i = 1; i <= ${params.rfam_seq.main_chunks}; i++)); do mv ${name}.fa.\$i ${name}_\$i.fa; done
  find . -name '*.fa' | xargs -P 4 -I {} gzip {}
  popd
  """
}

process build_rev {
  publishDir 'genomes/rfamseq', mode: 'copy'

  input:
  tuple path(merged), path(ssi)

  output:
  path 'rev-rfamseq*', emit: to_rev

  """
  mkdir to-rev
  export SEED=\$RANDOM
  echo \$SEED
  seqkit sample --rand-seed \$SEED -p ${params.rfam_seq.rev_fraction} merged.fa > sampled.fa
  seqkit split2 --out-dir to-rev -p ${params.rfam_seq.rev_chunks} sampled.fa
  pushd to-rev
  for((i = 1; i <= ${params.rfam_seq.rev_chunks}; i++)); do
    pretty_count="\$(printf '%03i' \$i)"
    esl-shuffle -r sampled.part_\$pretty_count.fa > rev-rfamseq${params.version}_\$i.fa
  done
  esl-seqstat ../sampled.fa > rev-rfamseq${params.version}-all.seqstat
  find . -name '*.fa' | xargs -P 4 -I {} gzip {}
  popd
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
    | merge_chunks \
    | (build_rfamseq & build_rev)
}

workflow {
  genome_download()
}
