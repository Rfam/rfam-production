#!/bin/bash

#SBATCH --time 4-0:0:0
#SBATCH --mem=2G
#SBATCH -p standard

module load nextflow-22.10.1-gcc-11.2.0-ju5saqw
nextflow run -ansi-log false main.nf
