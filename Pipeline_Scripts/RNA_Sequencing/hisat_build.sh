#!/usr/bin/env bash

#SBATCH --job-name=execute_hisat-build
#SBATCH --nodes=1
#SBATCH --ntasks=24 # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=10:00:00  # modify this to reflect how long to let the job go.
#SBATCH --output=../03_output/logfiles/log_hisat-build_%J.txt

hisat2-build -p 24 GCF_003704035.1_HU_Pman_2.1.3_genomic.fna pman