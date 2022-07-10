process fetch_ncbi_locations {
  queue 'short'
  errorStrategy 'finish'

  input:
  path(urls)

  output:
  path('ncbi.db')

  """
  mkdir summaries
  pushd summaries
  xargs -P 4 -a "../$urls" wget --no-verbose
  popd
  find summaries -type f | xargs cat > merged
  grep '^#' merged | tail -1 | sed 's|# ||' > info
  grep -v '^#' merged >> info
  parse_assembly_info.py info ncbi.db
  """
}

process find_genomes {
  queue 'short'

  input:
  path(to_skip)

  output:
  path("ncbi/*.jsonl"), emit: ncbi
  path("ena/*.jsonl"), emit: ena

  """
  curl '$params.proteome_xml' > summary.xml
  proteomes_to_genomes.py --ignorable $to_skip summary.xml .

  mkdir ncbi
  shuf ncbi.jsonl > ncbi-shuf.jsonl
  split -n l/550 ncbi-shuf.jsonl --additional-suffix='.jsonl' ncbi/

  mkdir ena
  shuf ena.jsonl > ena-shuf.jsonl
  split -n l/10 ena-shuf.jsonl --additional-suffix='.jsonl' ena/
  """
}

process download_ncbi {
  tag { "$gca_file.baseName" }
  maxForks 30
  queue 'short'
  publishDir "genomes/ncbi/${gca_file.baseName}"
  memory { 6.GB * task.attempt }
  errorStrategy { task.exitStatus in 129..140 ? 'retry' : 'finish' }
  maxRetries 4

  input:
  tuple path(gca_file), path(info)

  output:
  tuple val("ncbi-${gca_file.baseName}"), path("UP*.fa")

  """
  set -euo pipefail

  mkdir complete
  ncbi_urls.py $info $gca_file complete urls ena-only.jsonl
  xargs -a urls -L 2 -P 4 wget --no-verbose -O || true

  # It turns out that not all files which are specified by NCBI will actually
  # exist. This can be dealt with by falling back to ENA based lookup
  find complete -name '*.fa.gz' -empty > missing
  if [[ -s missing ]]; then
    xargs -a missing rm
    xargs -a missing -I {} basename {} \
    | cut -d. -f1 \
    | xargs -I {} jq -c 'select(.upi == "{}") | .accession = null | .kind = "ena"' $gca_file >> ena-only.jsonl
  fi

  gzip -d complete/*.fa.gz
  select_ids.py --ignore-file ena-only.jsonl $gca_file complete/
  if [[ -e ena-only.jsonl ]]; then
    download_ena.py ena-only.jsonl
  fi
  """
}

process download_ena {
  tag { "$ena_file.baseName" }
  maxForks 10
  publishDir 'genomes/ena'
  memory '5GB'
  queue 'short'
  errorStrategy 'finish'

  input:
  path(ena_file)

  output:
  tuple val("ena-${ena_file.baseName}"), path('UP*.fa')

  """
  download_ena.py $ena_file
  """
}

process validate_chunk {
  queue 'short'
  tag { "$short_name" }

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
  queue 'short'

  input:
  path('genomes*.fa')

  output:
  tuple path('merged.fa'), path('merged.fa.ssi')

  """
  set -euo pipefail

  find . -name 'genomes*.fa' | xargs cat > merged.fa
  esl-seqstat -a merged.fa > merged.seqstat
  esl-sfetch --index merged.fa
  """
}

process build_rfamseq {
  queue 'short'
  container ''
  publishDir 'genomes/rfamseq'
  errorStrategy 'finish'

  input:
  tuple path(merged), path(ssi)

  output:
  path('rfamseq/*')

  script:
  name = "r${params.rfam_seq.main_chunks}_rfamseq${params.version}"
  """
  mkdir rfamseq
  pushd rfamseq
  /homes/rfamprod/Bio-Easel/scripts/esl-ssplit.pl -v --oroot ${name}.fa -n -r -z ../${merged} ${params.rfam_seq.main_chunks}
  for((i = 1; i < ${params.rfam_seq.main_chunks}; i++)); do mv ${name}.fa.\$i ${name}_\$i.fa; done
  gzip *.fa
  popd
  """
}

process build_rev {
  queue 'short'
  publishDir 'genomes/rfamseq'
  errorStrategy 'finish'

  input:
  tuple path(merged), path(ssi)

  output:
  path 'to-rev/*', emit: to_rev

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
  gzip *.fa
  popd
  """
}

workflow genome_download {
  main:
    Channel.fromPath('ncbi-urls.txt') | fetch_ncbi_locations | set { ncbi_info }

    Channel.fromPath(params.ignore_upi) | find_genomes

    find_genomes.out.ncbi \
    | flatten \
    | combine(ncbi_info) \
    | download_ncbi \
    | set { ncbi }

    find_genomes.out.ena
    | flatten \
    | download_ena \
    | set { ena }

    ncbi.mix(ena) \
    | validate_chunk \
    | collect \
    | merge_chunks \
    | (build_rfamseq & build_rev)
}

workflow {
  genome_download()
}
