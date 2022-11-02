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
    queue 'short'

    input:
    val(acc)

    output:
    val('files_done')

    """
    python $params.rfamprod/scripts/export/fasta_file_generator.py --seq-db $params.rfamseqfa --rfam-seed $params.release_ftp/seed/Rfam.seed --outdir $params.release_ftp/fasta_files --acc $acc
    """
}

process combine_fasta {

    input:
    val('files_done')

    output:
    val('done')

    """
    cd $params.release_ftp/fasta_files
    mkdir unzipped
    cp *.fa.gz unzipped/
    cd unzipped
    gunzip *.gz
    cat *.fa > Rfam.fa
    cp Rfam.fa ../
    """
}

workflow generate_fasta_files {
  take: start
  emit: done
  main:
    start
    | fetch_families \
    | splitText \
    | generate_fasta \
    | set { done }
}

workflow {
  generate_fasta_files(Channel.of('start'))
}