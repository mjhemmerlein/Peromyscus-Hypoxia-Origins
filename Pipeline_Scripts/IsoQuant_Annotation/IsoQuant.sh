#!/usr/bin/env bash

#SBATCH --job-name=execute_isoquant
#SBATCH --nodes=1
#SBATCH --ntasks=24 # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=unlimited   # modify this to reflect how long to let the job go.
#SBATCH --output=/home/mjhemm/projects/IsoQuant/03_output/logfiles/log_isoquant_%J.txt

ulimit -S -n 131072

isoquant.py --reference GCF_003704035.1_HU_Pman_2.1.3_genomic.fna \
--genedb GCF_003704035.1_HU_Pman_2.1.3_genomic.gtf \
--complete_genedb \
--fastq Sample_1.hifi_reads.fastq.gz \
		Sample_2.hifi_reads.fastq.gz \
		Sample_3.hifi_reads.fastq.gz \
		Sample_4.hifi_reads.fastq.gz \
		Sample_5.hifi_reads.fastq.gz \
		Sample_6.hifi_reads.fastq.gz \
		Sample_7.hifi_reads.fastq.gz \
		Sample_8.hifi_reads.fastq.gz \
		Sample_9.hifi_reads.fastq.gz \
		Sample_10.hifi_reads.fastq.gz \
		Sample_11.hifi_reads.fastq.gz \
		Sample_12.hifi_reads.fastq.gz \
--data_type pacbio_ccs -o ../03_output