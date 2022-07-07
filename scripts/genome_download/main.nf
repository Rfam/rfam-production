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
  memory '5GB'
  memory { 4.GB * task.attempt }
  errorStrategy { task.exitStatus in 129..140 ? 'retry' : 'ignore' }
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
  tag { "$short_name" }

  input:
  tuple val(short_name), path('genomes*.fa')

  output:
  path("${short_name}.fa")

  """
  find . -name 'genomes*.fa' | xargs -I {} seqkit rmdup -s {} > unique.fa
  seqkit shuffle --two-pass unique.fa > ${short_name}.fa
  esl-sfetch --index ${short_name}.fa
  """
}

process merge_and_split {
  publishDir 'genomes'
  container ''

  input:
  path('genomes*.fa')

  output:
  path('rfamseq')

  """
  set -euo pipefail

  mkdir rfamseq to-rev
  find . -name 'genomes*.fa' | xargs cat > merged.fa
  esl-sfetch --index merged.fa

  pushd rfamseq
  esl-ssplit.pl -v --oroot r${params.rfam_seq.main_chunks}_rfamseq${params.version}.fa -n -r -z ../merged.fa ${params.rfam_seq.main_chunks}
  popd

  pushd to-rev
  esl-ssplit.pl -v --oroot rev-rfamseq${params.version}.fa -n -r -z ../merged.fa ${params.rfam_seq.rev_fraction}
  popd

  find to-rev -name '*.fa' -print '%f\\n" | xargs -I {} esl-shuffle -r -o rfamseq/{} to-rev/{}
  find rfamseq -name 'rev-*.fa' | xargs cat | esl-seqstat - > rfamseq/rev-rfamseq${params.version}-all.seqstat
  find rfamseq -name '*.fa' | xargs gzip
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
    | merge_and_split
}

workflow {
  genome_download()
}
