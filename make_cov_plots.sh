
Rscript plot_cov_loop.R \
  -n Alasoo_2018 \
  -f /gpfs/space/home/kerimov/qtlmap_lc_hisat/results/Alasoo_2018_17Feb/susie/Alasoo_2018_leafcutter_macrophage_naive.purity_filtered.txt.gz \
  -s /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/Alasoo_2018.tsv\
  -p /gpfs/space/home/kerimov/qcnorm_fast/results/Alasoo_2018_16Feb/Alasoo_2018_run/Alasoo_2018/normalised/leafcutter/leafcutter_metadata.txt.gz \
  -q macrophage_naive \
  -o /gpfs/space/home/kerimov/coverage_plots/plots \
  -v /gpfs/space/home/kerimov/coverage_plots/Alasoo_2018.filtered.vcf.gz \
  -b /gpfs/space/home/kerimov/rnaseq/results/Alasoo_2018_16Feb/bigwig \
  -m /gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v1.0/MANE_transcript_gene_map.txt \
  -g /gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v1.0/Homo_sapiens.GRCh38.105.gtf \
  -u /gpfs/space/home/kerimov/qcnorm_fast/results/Alasoo_2018_16Feb/Alasoo_2018_run/Alasoo_2018/normalised/leafcutter/qtl_group_split_norm/Alasoo_2018.macrophage_naive.tsv