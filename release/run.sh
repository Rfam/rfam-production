#!/bin/bash

#SBATCH --time 4-0:0:0
#SBATCH --mem=6G
#SBATCH -p production
#SBATCH -o main-vp.out
#SBATCH -e main-vp.err


export NXF_OPTS="-Dnxf.pool.type=sync -Dnxf.pool.maxThreads=10000"

module load nextflow
module load r/4.4.0

nextflow run -ansi-log false main.nf
#nextflow run -ansi-log false view_process.nf -resume
