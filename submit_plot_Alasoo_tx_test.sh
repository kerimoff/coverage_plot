#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --job-name="plot"
#SBATCH --partition=testing

module load any/jdk/1.8.0_265
module load nextflow
module load any/singularity/3.7.3
module load squashfs/4.4
module load tabix

nextflow run main.nf -profile tartu_hpc,test -resume \
  --studyFile /gpfs/space/home/kerimov/GitHub/coverage_plots/data/input_test_tx.tsv \
  --outdir /gpfs/space/home/kerimov/GitHub/coverage_plots/results/results_new_wiggleplotr_container 