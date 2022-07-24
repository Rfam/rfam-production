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

process fetch_uniprot_taxonomy {
  memory '2GB'
  queue 'short'
  errorStrategy 'finish'

  output:
  path('taxonomy.tsv')

  """
  curl 'https://rest.uniprot.org/taxonomy/stream?compressed=true&fields=id%2Ccommon_name%2Cscientific_name%2Clineage&format=tsv&query=%2A' > taxonomy.tsv.gz
  gzip -d taxonomy.tsv.gz
  """
}

process download_all_proteomes {
  queue 'short'

  output:
  path('summary.xml')

  """
  curl '$params.proteome_xml' > summary.xml
  """
}

process find_genomes {
  queue 'short'

  input:
  tuple path(summary), path(to_skip)

  output:
  path("summary.jsonl"), emit: summary
  path("parts/*.jsonl"), emit: chunks

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

// process find_genomes {
//   queue 'short'
//
//   input:
//   tuple path(summary), path(to_skip)
//
//   output:
//   path("summary.jsonl"), emit: summary
//   path("ncbi/*.jsonl"), emit: ncbi
//   path("ena/*.jsonl"), emit: ena
//
//   """
//   proteomes_to_genomes.py --ignorable $to_skip $summary .
//
//   mkdir ncbi
//   shuf ncbi.jsonl > ncbi-shuf.jsonl
//   split -n l/${params.ncbi_chunks} ncbi-shuf.jsonl --additional-suffix='.jsonl' ncbi/
//
//   mkdir ena
//   shuf ena.jsonl > ena-shuf.jsonl
//   split -n l/${params.ena_chunks} ena-shuf.jsonl --additional-suffix='.jsonl' ena/
//
//   cat ncbi.jsonl ena.jsonl > summary.jsonl
//   """
// }

process lookup_taxonomy_info {
  queue 'short'
  publishDir 'genomes/metadata', mode: 'copy'

  input:
  tuple path(summary), path(taxonomy)

  output:
  path('taxonomy.csv')

  """
  uniprot_taxonomy.py $taxonomy $summary taxonomy.csv
  """
}

process download {
  tag { "$proteome_file.baseName" }
  maxForks 30
  queue 'short'
  publishDir "genomes/fasta/${proteome_file.baseName}", mode: "copy"
  memory { 6.GB * task.attempt }
  errorStrategy 'ignore'
  /* errorStrategy { task.exitStatus in 129..140 ? 'retry' : 'finish' } */
  maxRetries 4

  input:
  tuple path(proteome_file), path(info)

  output:
  tuple val("${proteome_file.baseName}"), path("UP*.fa"), path("UP*.jsonl")

  """
  rfamseq download $info $proteome_file .
  """
}

// process download_ncbi {
//   tag { "$gca_file.baseName" }
//   maxForks 30
//   queue 'short'
//   publishDir "genomes/ncbi/${gca_file.baseName}"
//   memory { 6.GB * task.attempt }
//   errorStrategy { task.exitStatus in 129..140 ? 'retry' : 'finish' }
//   maxRetries 4
//
//   input:
//   tuple path(gca_file), path(info)
//
//   output:
//   tuple val("ncbi-${gca_file.baseName}"), path("UP*.fa"), path("metadata.db")
//
//   """
//   set -euo pipefail
//
//   mkdir complete
//   ncbi_urls.py $info $gca_file complete urls ena-only.jsonl
//   xargs -a urls -L 2 -P 4 wget --no-verbose -O || true
//
//   # It turns out that not all files which are specified by NCBI will actually
//   # exist. This can be dealt with by falling back to ENA based lookup for any
//   # empty downloads
//   find complete -name '*.fa.gz' -empty > missing
//   if [[ -s missing ]]; then
//     xargs -a missing rm
//     xargs -a missing -I {} basename {} \
//     | cut -d. -f1 \
//     | xargs -I {} jq -c 'select(.upi == "{}") | .accession = null | .kind = "ena"' $gca_file >> ena-only.jsonl
//   fi
//
//   gzip -d complete/*.fa.gz
//   select_ids.py --ignore-file ena-only.jsonl $gca_file complete/ metadata.db .
//   if [[ -e ena-only.jsonl ]]; then
//     download_ena.py ena-only.jsonl metadata.db .
//   fi
//   """
// }
//
// process download_ena {
//   tag { "$ena_file.baseName" }
//   maxForks 10
//   publishDir 'genomes/ena', mode: 'copy'
//   memory '5GB'
//   queue 'short'
//   errorStrategy 'retry'
//   maxRetries 4
//
//   input:
//   path(ena_file)
//
//   output:
//   tuple val("ena-${ena_file.baseName}"), path('UP*.fa'), path("metadata.db")
//
//   """
//   download_ena.py $ena_file metadata.db .
//   """
// }

process validate_chunk {
  queue 'short'
  tag { "$short_name" }

  input:
  tuple val(short_name), path('genomes*.fa')

  output:
  path("${short_name}.fa"), emit: merged
  tuple val(short_name), path("${short_name}.fa"), emit: for_metadata

  """
  set -euo pipefail

  find . -name 'genomes*.fa' | xargs -I {} seqkit rmdup -s {} | seqkit rmdup > unique.fa
  seqkit shuffle --two-pass unique.fa > ${short_name}.fa
  esl-sfetch --index ${short_name}.fa
  """
}

process generate_rfamseq_metadata {
  queue 'short'
  publishDir 'genomes/metadata', mode: 'copy'
  tag { "$short_name" }

  input:
  tuple val(short_name), path(fasta), path(info)

  output:
  path("rfamseq.${short_name}.tsv")

  """
  set -euo pipefail

  esl-seqstat -a --dna "$fasta" \
  | grep '^=' \
  | seqstat2rfamseq.py $info - rfamseq.${short_name}.tsv
  """
}

process generate_genseq_metadata {
  queue 'short'
  publishDir 'genomes/metadata', mode: 'copy'
  tag { "$short_name" }

  input:
  tuple val(short_name), path(fasta), path(info)

  output:
  path("rfamseq.${short_name}.tsv")

  """
  set -euo pipefail

  grep '^>' '$fasta' \
  | fasta2genseq.py $info $params.version - geseq.${short_name}.tsv
  """
}

process merge_chunks {
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
  queue 'short'
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
  queue 'short'
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

    download_all_proteomes | set { summary_xml }

    summary_xml | combine(to_ignore) | find_genomes

    find_genomes.out.chunks | flatten | combine(ncbi_info) | download

    /* find_genomes.out.ncbi \ */
    /* | flatten \ */
    /* | combine(ncbi_info) \ */
    /* | download_ncbi \ */
    /* | set { ncbi } */

    /* find_genomes.out.ena */
    /* | flatten \ */
    /* | download_ena \ */
    /* | set { ena } */

    /* find_genomes.out.summary \ */
    /* | combine(fetch_uniprot_taxonomy()) \ */
    /* | lookup_taxonomy_info */

    /* ncbi.mix(ena) | set { chunks } */

    /* chunks | generate_rfamseq_metadata */
    /* chunks | generate_genseq_metadata */

    /* chunks \ */
    /* | map { s, gs, _ -> [s, gs] } \ */
    /* | validate_chunk */

    /* validate_chunk.out.genome | collect | merge_genome */

    /* validate_chunk.out.merged */
    /* | collect \ */
    /* | merge_chunks \ */
    /* | (build_rfamseq & build_rev) */
}

workflow {
  genome_download()
}
