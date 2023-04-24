process fetch_ncbi_locations {
  queue 'short'
  time '1h'
  errorStrategy 'finish'

  input:
  path(urls)

  output:
  path('ncbi.db')

  """
  mkdir summaries
  wget --input-file ncbi-urls.txt -P summaries
  find summaries -type f | xargs cat > merged
  grep '^#' merged | tail -1 | sed 's|# ||' > info
  grep -v '^#' merged | sort -u >> info
  rfamseq parse-assembly-summary \
    summaries/assembly_summary_refseq.txt \
    summaries/assembly_summary_genbank.txt \
    summaries/assembly_summary_refseq_historical.txt \
    summaries/assembly_summary_genbank_historical.txt \
    ncbi.db
  """
}

process download_all_proteomes {
  queue 'short'
  time '1h'
  publishDir 'genomes/uniprot', mode: 'copy'

  output:
  path('summary.xml')

  """
  curl '$params.proteome_xml' > summary.xml
  """
}

process find_genomes {
  queue 'short'
  time '1h'

  input:
  tuple path(summary), path(to_skip)

  output:
  path("parts/*.jsonl")

  """
  rfamseq proteomes2genomes --ignore $to_skip $summary summary.jsonl

  mkdir parts
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
  maxForks 30
  publishDir "genomes/fasta/${proteome_file.baseName}", mode: "copy"
  memory { 6.GB * task.attempt }
  errorStrategy { task.exitStatus in 125..140 ? 'retry' : 'finish' }
  maxRetries 4

  input:
  tuple path(proteome_file), path(info)

  output:
  tuple val("${proteome_file.baseName}"), path("UP*.fa"), path("UP*.jsonl")

  """
  rfamseq download ${params.version} $info $proteome_file .
  """
}

process validate_chunk {
  tag { "$short_name" }
  queue 'standard'
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
  queue 'standard'
  time '24h'
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
  queue 'standard'
  container ''
  publishDir 'genomes/rfamseq', mode: 'copy'
  errorStrategy 'finish'

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
  queue 'standard'
  publishDir 'genomes/rfamseq', mode: 'copy'
  errorStrategy 'finish'

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
    Channel.fromPath(params.ignore_upi) | set { to_ignore }
    Channel.fromPath('ncbi-urls.txt') | fetch_ncbi_locations | set { ncbi_info }

    download_all_proteomes \
    | combine(to_ignore)
    | find_genomes \
    | flatten \
    | combine(ncbi_info) \
    | download \
    | map { s, gs, _ -> [s, gs] } \
    | validate_chunk \
    | collect \
    | merge_chunks \
    | (build_rfamseq & build_rev)
}

workflow {
  genome_download()
}
