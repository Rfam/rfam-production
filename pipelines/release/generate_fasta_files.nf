process fetch_families {

    input:
    val(_flag)

    output:
    file('families')

    """
    mysql -s `python ${params.rfamprod}/scripts/view/mysql_options.py $params.db` <<< "select rfam_acc from family" > families
    """
}

process generate_fasta {

    input:
    val(acc)

    output:
    val('done')

    """
    python $params.rfamprod/scripts/export/fasta_file_generator.py --seq-db $params.rfamseqfa --rfam-seed $params.release_ftp/seed/Rfam.seed --acc $acc --outdir $params.release_ftp/fasta_files
    """
}

workflow generate_fasta_files {
  take: start
  emit: done
  main:
    start
    | fetch_families \
    | splitCsv(sep: "\t") \
    | generate_fasta \
    | set { done }
}

workflow {
  generate_fasta_files(Channel.of('start'))
}