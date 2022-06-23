process fetch_ncbi_locations {
  errorStrategy 'finish'

  output:
  path('ncbi.pickle')

  """
  curl '$params.genbank_assembly_info' > initial
  grep '^#' initial | tail -1 | sed 's|# ||' > info
  grep -v '^#' initial >> info
  {
    curl '$params.genbank_old_assembly_info' 
    curl '$params.refseq_assembly_info' 
    curl '$params.refseq_old_assembly_info' 
  } | grep -v '^#' >> info
  parse_assembly_info.py info ncbi.pickle
  """
}

process find_genomes {
  input:
  path(to_skip)

  output:
  path("ncbi.csv"), emit: ncbi
  path("ena.csv"), emit: ena

  """
  curl '$params.proteome_xml' > summary.xml
  proteomes_to_genomes.py --ignorable $to_skip summary.xml
  """
}

process download_gca {
  tag { "$proteome" }
  maxForks 20
  errorStrategy 'ignore'
  publishDir 'genomes'

  input:
  tuple val(proteome), val(kind), val(gca), path(info)

  output:
  path("${proteome}.fa")

  """
  ncbi_url.py $info $gca | xargs -I {} wget -O - {} | gzip -d > ${proteome}.fa
  """
}

process download_ena {
  tag { "$proteome-$accession" }
  maxForks 10
  errorStrategy 'ignore'

  input:
  tuple val(proteome), val(kind), val(accession)

  output:
  tuple val(proteome), path("${proteome}-${accession}.fa")

  """
  curl 'https://www.ebi.ac.uk/ena/browser/api/fasta/${accession}?download=true' > ${proteome}-${accession}.fa
  """
}

process merge_ena_fasta {
  tag { "$proteome" }
  errorStrategy 'ignore'
  publishDir 'genomes'

  input:
  tuple val(proteome), path('parts*.fa')

  output:
  path("${proteome}.fa")
  
  """
  cat parts*.fa > ${proteome}.fa
  """
}

workflow genome_download {
  main:
    fetch_ncbi_locations | set { ncbi }

    Channel.fromPath(params.ignore_upi) \
    | find_genomes \
    | splitCsv \
    | branch {
      ncbi: it[1] == 'ncbi'
      ena: it[1] == 'other'
    } \
    | set { to_download }

    to_download.ncbi | combine(ncbi) | download_gca

    to_download.ena \
    | download_ena \
    | groupTuple \
    | merge_ena_fasta
}

workflow {
  genome_download()
}
