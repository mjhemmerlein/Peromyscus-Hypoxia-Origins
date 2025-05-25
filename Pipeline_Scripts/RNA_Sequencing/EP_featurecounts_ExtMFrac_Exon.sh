#!/usr/bin/env bash
#SBATCH --job-name=execute_featureCounts
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=unlimited
#SBATCH --output=../logfiles/log_featureCounts_MMFrac_%J.txt
## Tabulate counts from .bam files using genome annotation file

featureCounts -O -M --fraction -p -F GTF -t exon -T 24 -a OUT.extended_annotation.gtf \
    -o ../FeatureCounts_Ext_MultimapFrac_Exon/EP_Pman_ExtMMFrac_readcounts_Exon.txt \
    *.bam