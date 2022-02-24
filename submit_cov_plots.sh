#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --job-name="cov_plot"
#SBATCH --partition=amd

Rscript plot_cov_loop.R \
  -n Alasoo_2018 \
  -f /gpfs/space/home/kerimov/qtlmap_lc_hisat/results/Alasoo_2018_17Feb/susie/Alasoo_2018_leafcutter_macrophage_naive.purity_filtered.txt.gz \
  -s /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/Alasoo_2018.tsv\
  -p /gpfs/space/home/kerimov/qcnorm_fast/results/Alasoo_2018_16Feb/Alasoo_2018_run/Alasoo_2018/normalised/leafcutter/leafcutter_metadata.txt.gz \
  -q macrophage_naive \
  -o /gpfs/space/home/kerimov/coverage_plots/plots_new_pval/macrophage_naive \
  -v /gpfs/space/home/kerimov/coverage_plots/data/Alasoo_2018.filtered.vcf.gz \
  -b /gpfs/space/home/kerimov/rnaseq/results/Alasoo_2018_16Feb/bigwig \
  -m /gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v1.0/MANE_transcript_gene_map.txt \
  -g /gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v1.0/Homo_sapiens.GRCh38.105.gtf \
  -u /gpfs/space/home/kerimov/qcnorm_fast/results/Alasoo_2018_16Feb/Alasoo_2018_run/Alasoo_2018/normalised/leafcutter/qtl_group_split_norm/Alasoo_2018.macrophage_naive.tsv \
  -a /gpfs/space/home/kerimov/qtlmap_lc_hisat/results/Alasoo_2018_17Feb/sumstats/Alasoo_2018_leafcutter_macrophage_naive.all.tsv.gz
  

# Rscript plot_cov_loop.R \
#   -n Alasoo_2018 \
#   -f /gpfs/space/home/kerimov/qtlmap_lc_hisat/results/Alasoo_2018_17Feb/susie/Alasoo_2018_leafcutter_macrophage_IFNg.purity_filtered.txt.gz \
#   -s /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/Alasoo_2018.tsv\
#   -p /gpfs/space/home/kerimov/qcnorm_fast/results/Alasoo_2018_16Feb/Alasoo_2018_run/Alasoo_2018/normalised/leafcutter/leafcutter_metadata.txt.gz \
#   -q macrophage_IFNg \
#   -o /gpfs/space/home/kerimov/coverage_plots/plots_new/macrophage_IFNg \
#   -v /gpfs/space/home/kerimov/coverage_plots/data/Alasoo_2018.filtered.vcf.gz \
#   -b /gpfs/space/home/kerimov/rnaseq/results/Alasoo_2018_16Feb/bigwig \
#   -m /gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v1.0/MANE_transcript_gene_map.txt \
#   -g /gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v1.0/Homo_sapiens.GRCh38.105.gtf \
#   -u /gpfs/space/home/kerimov/qcnorm_fast/results/Alasoo_2018_16Feb/Alasoo_2018_run/Alasoo_2018/normalised/leafcutter/qtl_group_split_norm/Alasoo_2018.macrophage_IFNg.tsv


# Rscript plot_cov_loop.R \
#   -n Alasoo_2018 \
#   -f /gpfs/space/home/kerimov/qtlmap_lc_hisat/results/Alasoo_2018_17Feb/susie/Alasoo_2018_leafcutter_macrophage_IFNg+Salmonella.purity_filtered.txt.gz \
#   -s /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/Alasoo_2018.tsv\
#   -p /gpfs/space/home/kerimov/qcnorm_fast/results/Alasoo_2018_16Feb/Alasoo_2018_run/Alasoo_2018/normalised/leafcutter/leafcutter_metadata.txt.gz \
#   -q macrophage_IFNg+Salmonella \
#   -o /gpfs/space/home/kerimov/coverage_plots/plots_new/macrophage_IFNg+Salmonella \
#   -v /gpfs/space/home/kerimov/coverage_plots/data/Alasoo_2018.filtered.vcf.gz \
#   -b /gpfs/space/home/kerimov/rnaseq/results/Alasoo_2018_16Feb/bigwig \
#   -m /gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v1.0/MANE_transcript_gene_map.txt \
#   -g /gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v1.0/Homo_sapiens.GRCh38.105.gtf \
#   -u /gpfs/space/home/kerimov/qcnorm_fast/results/Alasoo_2018_16Feb/Alasoo_2018_run/Alasoo_2018/normalised/leafcutter/qtl_group_split_norm/Alasoo_2018.macrophage_IFNg+Salmonella.tsv


# Rscript plot_cov_loop.R \
#   -n Alasoo_2018 \
#   -f /gpfs/space/home/kerimov/qtlmap_lc_hisat/results/Alasoo_2018_17Feb/susie/Alasoo_2018_leafcutter_macrophage_Salmonella.purity_filtered.txt.gz \
#   -s /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/Alasoo_2018.tsv\
#   -p /gpfs/space/home/kerimov/qcnorm_fast/results/Alasoo_2018_16Feb/Alasoo_2018_run/Alasoo_2018/normalised/leafcutter/leafcutter_metadata.txt.gz \
#   -q macrophage_Salmonella \
#   -o /gpfs/space/home/kerimov/coverage_plots/plots_new/macrophage_Salmonella \
#   -v /gpfs/space/home/kerimov/coverage_plots/data/Alasoo_2018.filtered.vcf.gz \
#   -b /gpfs/space/home/kerimov/rnaseq/results/Alasoo_2018_16Feb/bigwig \
#   -m /gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v1.0/MANE_transcript_gene_map.txt \
#   -g /gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v1.0/Homo_sapiens.GRCh38.105.gtf \
#   -u /gpfs/space/home/kerimov/qcnorm_fast/results/Alasoo_2018_16Feb/Alasoo_2018_run/Alasoo_2018/normalised/leafcutter/qtl_group_split_norm/Alasoo_2018.macrophage_Salmonella.tsv