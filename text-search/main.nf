// Finds all genomes and chunks them into params.export.genome_chunks files
process fetch_genomes {
  label 'db'

  input:
  val(_flag)

  output:
  path('genomes/*.txt')

  """
  rfam_export list-genomes all-genomes.txt
  mkdir genomes
  shuf all-genomes.txt > genomes-shuf.txt
  split \
    -n l/${params.export.genome_chunks} \
    --additional-suffix='.txt' \
    genomes-shuf.txt genomes/
  """
}

// Export all families in one process. This does them all at once since they are
// so small and there aren't very may of them.
process export_families {
  label 'db'
  publishDir "${params.text_search}/families", mode: "copy"
  errorStrategy 'ignore'

  input:
  val(_flag)

  output:
  path("${accession}.xml")

  """
  rfam_export families families/
  find families -type f -name '*.xml' | xargs -I {} xmllint --noout --schema http://www.ebi.ac.uk/ebisearch/XML4dbDumps.xsd {}
  """
}

// Export all seed and full hits in the input genomes file.
process export_hits {
  tag { "$genomes.baseName" }
  label 'db'
  memory 10.GB
  publishDir "${params.text_search}/full_region", mode: "copy"
  errorStrategy 'ignore'

  input:
  path(genomes)

  output:
  path("full_regions/*.xml")

  """
  rfam_export full-regions $genomes full_regions
  find full_regions -type f -name '*.xml' | xargs -I {} xmllint --noout --schema http://www.ebi.ac.uk/ebisearch/XML4dbDumps.xsd {}
  """
}

// Export all genome data in the
process export_genomes {
  tag { "$genomes.baseName" }
  label 'db'
  memory 10.GB
  publishDir "${params.text_search}/genomes", mode: "copy"
  errorStrategy 'ignore'

  input:
  path(genomes)

  output:
  path("genomes/*.xml")

  """
  rfam_export genomes $genomes genomes
  find genomes -type f -name '*.xml' | xargs -I {} xmllint --noout --schema http://www.ebi.ac.uk/ebisearch/XML4dbDumps.xsd {}
  """
}

// Export all motifs
process export_motifs {
  label 'db'
  publishDir "${params.text_search}/motifs", mode: "copy"
  errorStrategy 'ignore'

  input:
  val(_flag)

  output:
  path("motifs/*.xml")

  """
  rfam_export export-motifs motifs/
  find motifs -type f -name '*.xml' | xargs -I {} xmllint --noout --schema http://www.ebi.ac.uk/ebisearch/XML4dbDumps.xsd {}
  """
}

workflow export {
  take:
    start

  main:
    // Export all families
    start \
    | export_families \
    | map { "family done" } \
    | set { family_done }

    // Export all motifs
    start \
    | export_motifs \
    | map { "motif done" } \
    set { motifs_done }

    // Compute the genomes to export
    start \
    | fetch_genomes \
    | flatten \
    | set { genomes }

    // Export all genomes
    genomes \
    | export_genomes \
    | collect \
    | map { "genomes done" } \
    | set { genome_done }

    // Export all full regions
    genomes \
    | export_hits \
    | collect \
    | map { "hits done" } \
    | set { hits_done }

    // Collect all results and produce a single output result
    family_done \
    | mix(hits_done, genome_done, motifs_done, hits_done) \
    | collect \
    | map { "export done" } \
    | set { done }

  emit:
    done
}

workflow {
  export(Channel.of('start'))
}
