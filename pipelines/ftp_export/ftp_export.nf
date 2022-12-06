include { fromQuery } from 'plugin/nf-sqldb'

process generate_seed_files {

    input:
    val(acc)

    output:
    path('*.seed')

    """
    writeAnnotatedSeed.pl $acc
    mv $acc > $acc.seed
    """

}

workflow {

  ch = channel.fromQuery('select rfam_acc from family', db: 'rfamlive')
  ch \
  | view
  | splitText \
  | generate_seed_files \
  | collectFile(name: 'Rfam.seed', newLine: true)


}
