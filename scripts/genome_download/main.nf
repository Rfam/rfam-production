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
  maxForks 20
  publishDir "genomes/ncbi/${gca_file.baseName}"

  input:
  tuple path(gca_file), path(info)

  output:
  path("*.fa")

  """
  set -euo pipefail

  mkdir complete
  ncbi_urls.py $info $gca_file complete urls ena-only.jsonl
  xargs -a urls -L 2 -P 4 wget -O
  gzip -d complete/*.fa.gz
  select_ids.py $gca_file complete/
  download_ena.py ena_only.jsonl
  """
}

process download_ena {
  tag { "$ena_file.baseName" }
  maxForks 10
  publishDir 'genomes/ena'

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
