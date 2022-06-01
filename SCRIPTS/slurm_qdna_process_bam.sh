#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --job-name=QDNA_process_bam
#SBATCH --mem=
#SBATCH --output=QDNA_process_bam.out
cd /rds/project/rds-csoP2nj6Y6Y/ctDNA/dl

module load R/4.2.0-icelake

Rscript 2_qdna_process_bam.R