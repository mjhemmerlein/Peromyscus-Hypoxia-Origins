#!/usr/bin/env bash
#SBATCH --job-name=execute_featureCounts
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=unlimited
#SBATCH --output=../logfiles/log_featureCounts_ExtMM_%J.txt
## Tabulate counts from .bam files using genome annotation file

featureCounts -O -M -p -F GTF -t transcript -T 24 -a OUT.extended_annotation.gtf \
    -o ../FeatureCounts_Ext_Multimap/EP_Pman_ExtMM_readcounts.txt \
    *.bam