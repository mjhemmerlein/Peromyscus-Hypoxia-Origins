#!/usr/bin/env bash
#SBATCH --job-name=execute_SendJob
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=unlimited
#SBATCH --output=/home/mjhemm/projects/Hypoxia_Dream/Peromyscus-Hypoxia-Origins/logfiles/log_SendJob_%J.txt

Rscript RNA_Seq_RScripts/EP_Dream.R