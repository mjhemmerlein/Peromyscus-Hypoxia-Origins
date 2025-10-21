#!/usr/bin/env bash
#SBATCH --job-name=execute_hisat-align
#SBATCH --nodes=1
#SBATCH --ntasks=24 # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=unlimited
#SBATCH --output=../logfiles/log_hisat-align_%J.txt

ls *.fastq.gz | while read file; do
        name=$(echo "${file}" | cut -d "_" -f 1-5)

        echo ${name}

hisat2 -p 24 -q --mp 2,0 \
-x ../../01_input/pman \
--summary-file ../reports/${name}.txt \
-1 ${name}_1_fastp.1.fastq.gz \
-2 ${name}_2_fastp.2.fastq.gz | \
samtools view -Sbo ../hisat2_align_paired/${name}_fastp_pman_Halign_liberal.bam

done

echo "Done!"