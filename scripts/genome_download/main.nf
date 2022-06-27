process fetch_ncbi_locations {
  errorStrategy 'finish'

  output:
  path('ncbi.db')

  """
  curl '$params.genbank_assembly_info' > initial
  grep '^#' initial | tail -1 | sed 's|# ||' > info
  grep -v '^#' initial >> info
  {
    curl '$params.genbank_old_assembly_info'
    curl '$params.refseq_assembly_info'
    curl '$params.refseq_old_assembly_info'
  } | grep -v '^#' >> info
  parse_assembly_info.py info ncbi.db
  """
}

process find_genomes {
  input:
  path(to_skip)

  output:
  path("ncbi/*.jsonl"), emit: ncbi
  path("ena/*.jsonl"), emit: ena

  """
  curl '$params.proteome_xml' > summary.xml
  proteomes_to_genomes.py --ignorable $to_skip summary.xml .

  mkdir ncbi
  split -n l/500 ncbi.jsonl --additional-suffix='.jsonl' ncbi/

  mkdir ena
  split -n l/10 ena.jsonl --additional-suffix='.jsonl' ena/
  """
}

process download_ncbi {
  tag { "$gca_file.baseName" }
  maxForks 10
  publishDir "genomes/ncbi/${gca_file.baseName}"
  memory '5GB'
  memory { 4.GB * task.attempt }
  errorStrategy { task.exitStatus in 130..140 ? 'retry' : 'finish' }
  maxRetries 4

  input:
  tuple path(gca_file), path(info)

  output:
  path("*.fa")

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
    | xargs -I {} jq 'select(.upi == "{}") | .accession = null | .kind = "ena"' $gca_file >> ena-only.jsonl
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
  errorStrategy 'finish'

  input:
  path(ena_file)

  output:
  path('*.fa')

  """
  download_ena.py $ena_file
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
  find . -name 'genomes*.fa' | xargs -I {} cat {} > merged.fa

  pushd rfamseq
  esl-randomize-sqfile.pl -O r${params.rfam_seq.main_chunks}_rfamseq${params.version}.fa -L -N ${params.rfam_seq.main_chunks} ../merged.fa 1.0
  popd

  pushd to-rev
  esl-randomize-sqfile.pl -O rev-rfamseq${params.version}.fa -L -N ${params.rfam_seq.rev_chunks} ../merged.fa ${params.rfam_seq.rev_fraction}
  popd

  find to-rev -name '*.fa' -print '%f\\n" | xargs -I {} esl-shuffle -r -o rfamseq/{} to-rev/{}
  find rfam-seq -name 'rev-*.fa' | xargs cat | esl-seqstat - > rfamseq/rev-rfamseq${params.version}-all.seqstat
  find rfam-seq -name '*.fa' | xargs gzip
  """
}

workflow genome_download {
  main:
    fetch_ncbi_locations | set { ncbi_info }

    Channel.fromPath(params.ignore_upi) | find_genomes

    find_genomes.out.ncbi \
    | flatten \
    | combine(ncbi_info) \
    | download_ncbi
    | set { ncbi }

    find_genomes.out.ena \
    | flatten \
    | download_ena \
    | set { ena }

    ncbi.mix(ena) \
    | collect \
    | merge_and_split
}

workflow {
  genome_download()
}
