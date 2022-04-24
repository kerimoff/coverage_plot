#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --job-name="plot"
#SBATCH --partition=amd

module load any/jdk/1.8.0_265
module load nextflow
module load any/singularity/3.7.3
module load squashfs/4.4
module load tabix

nextflow run main.nf -profile tartu_hpc -resume \
  --studyFile input_test_Alasoo_2018_all_quants.tsv \
  --outdir test_results/results_complete_Alasoo_2018 \
  --debug_mode


