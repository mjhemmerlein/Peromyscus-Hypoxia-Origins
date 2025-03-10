#!/usr/bin/env bash

#SBATCH --job-name=execute_fastp
#SBATCH --nodes=1
#SBATCH --ntasks=24 # modify this number to reflect how many cores you want to use (up to 24)
#SBATCH --time=10:00:00   # modify this to reflect how long to let the job go.
#SBATCH --output=../03_output/logfiles/log_fastp_%J.txt

ls *.fastq.gz | while read file; do
        name=$(echo "${file}" | cut -d "_" -f 1-4)

        echo ${name}

        fastp\
        --in1 ${name}_R1_001.fastq.gz \
        --out1 ../03_output/fastp_output/${name}_1_fastp.1.fastq.gz \
        --in2 ${name}_R2_001.fastq.gz \
   		--out2 ../03_output/fastp_output/${name}_2_fastp.2.fastq.gz \
    	--html ../03_output/reports/${name}.fastp.html \
		--cut_front \
		--cut_front_window_size=1 \
		--cut_front_mean_quality=3 \
		--cut_tail \
		--cut_tail_window_size=1  \
		--cut_tail_mean_quality=3 \
		--cut_right \
		--cut_right_window_size=5 \
		--cut_right_mean_quality=20 \
		--length_required=36

done

echo "Done!"