Rfamseq Update
==============

This is a nextflow pipeline meant to update the rfamseq database.

To run do:

```bash
nextflow run main.nf
```

This is containerized so setting `process.container` and `docker.enabled` are a
good idea.

What this does
--------------

The general task is to fetch the list of [reference proteomes](https://www.uniprot.org/help/reference_proteome) from UniProt.
This list is then parsed to extract the genome and all [components](https://legacy.uniprot.org/help/proteome_component) that were used in by UniProt.
Generally these genomes are GCA or GCF accessions and can be fetched easily from NCBI/ENA.
However, there is a complication as Uniprot may or may not use all chromosomes, scaffolds, etc that are in any given genome.
Because of this, each genome needs to be filtered both to remove any unneeded sequences (eg patch regions in human are generally not used) and then any missing components need to be fetched.

Additionally, the components can be [Whole Genome Sets](https://www.ncbi.nlm.nih.gov/genbank/wgs/) (WGS).
When this occurs the WGS set needs to be resolved into individual sequences, which can be at any assembly level from small 50mers to large contigs.
This tries to correctly resolve WGS sets into the largest possible sequence as well.
