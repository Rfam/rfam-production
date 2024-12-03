#!/usr/bin/env nextflow

process fill_template {
  maxForks 1
  publishDir "${params.release_ftp}", mode: "copy"

  input:
  tuple val(name), path(template_path), val(_flag)

  output:
  tuple path("${template_path.baseName}"), val('done')

  """
  templates.py ${name} ${template_path} ${template_path.baseName}
  """
}

workflow templates {
  take:
    flag
  emit:
    done
  main:
    Channel.fromPath('templates/*.template') \
    | map { [it.baseName, it] } \
    | combine(flag)
    | fill_template \
    | collect \
    | set { done }
}
