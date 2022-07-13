#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --job-name="data_plot"
#SBATCH --partition=amd

module load any/jdk/1.8.0_265
module load nextflow
module load any/singularity/3.7.3
module load squashfs/4.4
module load tabix

nextflow run plot_from_data.nf -profile tartu_hpc,test -resume \
  --from_files "/gpfs/space/home/kerimov/GitHub/coverage_plots/data/plot_data_tar_rds/*" \
  --outdir /gpfs/space/home/kerimov/GitHub/coverage_plots/results/results_plot_from_test_data_new


