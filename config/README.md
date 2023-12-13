# Configs

`gen_config` contains info for the genome-related scripts. This may be unneccessary if we update the process for downloading genomes and creating rfamseq in future releases (e.g. 15.0).  
`rfam_local` and `rfam_config` were tidied suring the codon migration, but will need further clean up for the move to SLURM, a lot of the variables are LSF specifics e.g. bsub groups, queue names.
`rfam_search` is used by the XML dumper scripts. 
