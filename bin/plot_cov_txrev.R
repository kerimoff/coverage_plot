message(" ## Loading libraries: optparse")
suppressPackageStartupMessages(library("optparse"))

#Parse command-line options
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  make_option(c("-f", "--finemap_susie"), type="character", default=NULL,
              help="Purity filtered susie output. Tab separated file", metavar = "type"),
  make_option(c("-s", "--sample_meta"), type="character", default=NULL,
              help="Sample metadata file. Tab separated file", metavar = "type"),
  make_option(c("-p", "--phenotype_meta"), type="character", default=NULL,
              help="Phenotype metadata file. Tab separated file", metavar = "type"),
    make_option(c("-q", "--qtl_group"), type="character", default=NULL,
              help="The selected qtl_group in the study", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default="./normalised_results/",
              help="Path to the output directory. [default \"%default\"]", metavar = "type"),
  make_option(c("-n", "--name_of_study"), type="character", default=NULL,
              help="Name of the study. Optional", metavar = "type"),
  make_option(c("-v", "--vcf_file"), type="character", default=NULL,
              help="TPM quantile TSV file with phenotype_id column", metavar = "type"),
  make_option(c("-b", "--bigwig_path"), type="character", default=NULL,
              help="Path to the bigwig files", metavar = "type"),
  make_option(c("-m", "--mane_transcript_gene_map"), type="character", default=NULL,
              help="Path to the MANE transcripts of genes map file", metavar = "type"),
  make_option(c("-g", "--gtf_file"), type="character", default=NULL,
              help="Path to the GTF file to get exons of transcripts", metavar = "type"),
  make_option(c("-t", "--txrev_gtf_file"), type="character", default=NULL,
              help="Path to the GTF file to get exons of txrev transcripts", metavar = "type"),
  make_option(c("-u", "--usage_matrix_norm"), type="character", default=NULL,
              help="Path to the normalised usage matrix", metavar = "type"),
  make_option(c("-a", "--all_summ_stats"), type="character", default=NULL,
              help="Path to the full gene nominal summary statistics", metavar = "type"),
  make_option(c("-e", "--exon_summ_stats"), type="character", default=NULL,
              help="Path to the full exon nominal summary statistics", metavar = "type"),
  make_option(c("-w", "--wiggle_plotr_path"), type="character", default=NULL,
              help="Local path to wiggleplotr package", metavar = "type"),
  make_option(c("-i", "--individual_boxplots"), type="logical", 
              action = "store_true", default=FALSE,
              help="Flag to generate individual boxplots", metavar = "type"),
  make_option(c("-d", "--debug_mode"), type="logical", 
              action = "store_true", default=FALSE,
              help="If to run the script in debug mode", metavar = "type")
)

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, cqn, data.table")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ggplot2"))
library(seqminer)

make_transcript_exon_granges <- function(gff, transcript_ids) {
  exon_list <- list()
  for (transcript_id in transcript_ids) {
    transcript_exons_temp <- gff[(base::gsub("\\..*","",SummarizedExperiment::elementMetadata(gff)[,"transcript_id"]) == transcript_id)]
    gene_id = transcript_exons_temp$gene_name[1]
    exon_list[[paste0("GENE:", gene_id, ":", transcript_id)]] <- transcript_exons_temp
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

make_transcript_exon_granges_ccds <- function(gff, transcript_ids) {
  exon_list <- list()
  for (transcript_id in transcript_ids) {
    transcript_exons_temp <- gff[(base::gsub("\\..*","",SummarizedExperiment::elementMetadata(gff)[,"transcript_id"]) == transcript_id & !is.na(SummarizedExperiment::elementMetadata(gff)[,"ccds_id"]))]
    gene_id = transcript_exons_temp$gene_name[1]
    exon_list[[paste0("GENE:", gene_id, ":", transcript_id)]] <- transcript_exons_temp
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

make_txrev_exon_granges <- function(gff, txrev_pheno_ids) {
  exon_list <- list()
  for (transcript_id in txrev_pheno_ids) {
    exon_list[[transcript_id]] <- gff[unlist(SummarizedExperiment::elementMetadata(gff)[,"Parent"]) == transcript_id]
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

generate_beta_plot <- function(transcript_struct_df_loc, 
                               nom_exon_cc_sumstats_filt_loc, 
                               limits, 
                               vert_lines = NULL,
                               vert_line_alpha = 0.3){
  exon_exp_rescaled_exons <- transcript_struct_df_loc %>% 
    dplyr::filter(transcript_id == paste0("GENE:",ss_oi$gene_name)) %>% 
    dplyr::mutate(exon_row_num = dplyr::row_number()) %>%
    dplyr::mutate(exon_rescaled_center = round((end+start)/2)) %>%
    dplyr::mutate(beta_wrap_label = "BETA with 95% CI") %>%
    dplyr::left_join(nom_exon_cc_sumstats_filt_loc %>% dplyr::select(exon_row_num, molecular_trait_id, beta, interval, p_fdr), by = "exon_row_num") %>% 
    dplyr::mutate(is_sign_exon_qtl = if_else(p_fdr <= 0.01, 1, 0.2)) 
  
  beta_plot <- ggplot(exon_exp_rescaled_exons, aes(x = exon_rescaled_center, y = beta, ymin = beta - interval, ymax = beta + interval, alpha = is_sign_exon_qtl)) + 
    ggplot2::geom_blank() +
    ggplot2::geom_point(color = "darkblue", alpha=exon_exp_rescaled_exons$is_sign_exon_qtl) + 
    ggplot2::geom_errorbar(width = max(limits)/200, alpha=exon_exp_rescaled_exons$is_sign_exon_qtl) + 
    ggplot2::geom_hline(yintercept = 0, alpha = 0.5, color = "lightgrey") +
    ggplot2::facet_grid(beta_wrap_label~.) +
    ggplot2::theme_light() +
    ggplot2::scale_x_continuous(expand = c(0,0)) + 
    ggplot2::ylab("Effect size") +
    coord_cartesian(xlim = limits) +
    theme(plot.margin=unit(c(0,1,0,1),"line"), 
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(colour = "grey10"),
          strip.background = element_rect(fill = "grey85")) 
  
  if (!is.null(vert_lines)) {
    beta_plot <- beta_plot + ggplot2::geom_vline(xintercept = vert_lines, alpha = vert_line_alpha, color = "lightgrey")
  }
  return(beta_plot)
}


#Debugging
if (FALSE) {
  opt = list()
  opt$n = "Alasoo_2018"
  opt$f = "/Users/kerimov/Work/temp_files/txrevise_stuff/Alasoo_2018_3Apr/susie/Alasoo_2018_txrev_macrophage_IFNg.purity_filtered.txt.gz"
  opt$s = "/Users/kerimov/Work/GitHub/SampleArcheology/studies/cleaned/Alasoo_2018.tsv"
  opt$p = "/Users/kerimov/Work/temp_files/txrevise_stuff/Alasoo_2018_3Apr/txrevise_Ensembl_105_phenotype_metadata.tsv.gz"
  opt$q = "macrophage_IFNg"
  opt$o = "/Users/kerimov/Work/temp_files/txrevise_stuff/Alasoo_2018_3Apr/plots_data_generate"
  opt$v = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018/vcf/Alasoo_2018_filtered_noX.vcf.gz"
  opt$b = "/Users/kerimov/Work/temp_files/txrevise_stuff/Alasoo_2018_3Apr/bigwig"
  opt$m = "/Users/kerimov/Work/temp_files/leafcutter_data/MANE_transcript_gene_map.txt"
  opt$g = "/Users/kerimov/Work/temp_files/leafcutter_data/Homo_sapiens.GRCh38.105.gtf"
  opt$t = "/Users/kerimov/Work/temp_files/txrevise_stuff/Homo_sapiens.GRCh38.105.txrevise.CAGE25/txrevise.all.gff3"
  opt$u = "/Users/kerimov/Work/temp_files/txrevise_stuff/Alasoo_2018_3Apr/qcnorm/Alasoo_2018.macrophage_IFNg.tsv.gz"
  opt$a = "/Users/kerimov/Work/temp_files/txrevise_stuff/Alasoo_2018_3Apr/sumstats/Alasoo_2018_txrev_macrophage_IFNg.all.tsv.gz"
  opt$e = "/Users/kerimov/Work/temp_files/txrevise_stuff/Alasoo_2018_3Apr/sumstats/Alasoo_2018_exon_macrophage_IFNg.all.tsv.gz"
  opt$w = "/Users/kerimov/Work/R_env/wiggleplotr/"
  index = 1
}

susie_file_path = opt$f
sample_meta_path = opt$s
phenotype_meta_path = opt$p
output_dir = opt$o
qtl_group_in = opt$q
study_name = opt$n
vcf_file_path = opt$v
bigwig_files_path = opt$b
mane_transcript_gene_map_file = opt$m
gtf_file_path = opt$g
txrev_gtf_file_path = opt$t
norm_usage_matrix_path = opt$u
nominal_sumstats_path = opt$a
nominal_exon_sumstats_path = opt$e
wiggleplotr_path = opt$w
individual_boxplots = opt$individual_boxplots
debug_mode = opt$debug_mode

message("######### Options: ######### ")
message("######### Working Directory  : ", getwd())
message("######### qtl_group          : ", qtl_group_in)
message("######### susie_file_path    : ", susie_file_path)
message("######### sample_meta_path   : ", sample_meta_path)
message("######### phenotype_meta_path: ", phenotype_meta_path)
message("######### output_dir         : ", output_dir)
message("######### opt_study_name     : ", study_name)
message("######### vcf_file_path      : ", vcf_file_path)
message("######### bigwig_files_path  : ", bigwig_files_path)
message("######### mane_map_file_path : ", mane_transcript_gene_map_file)
message("######### gtf_file_path      : ", gtf_file_path)
message("######### txrev_gtf_file_path: ", txrev_gtf_file_path)
message("######### norm_usage_matrix  : ", norm_usage_matrix_path)
message("######### nominal_sumstats   : ", nominal_sumstats_path)
message("######### exon_sumstats      : ", nominal_exon_sumstats_path)
message("######### wiggleplotr_path   : ", wiggleplotr_path)
message("######### individual_boxplots: ", individual_boxplots)
message("######### debug_mode         : ", debug_mode)
if (!is.null(wiggleplotr_path)) {
  devtools::load_all(wiggleplotr_path)
}

################## Global variable definitions ################
conf.level = 0.95
ci.value <- -qnorm( ( 1 - conf.level ) / 2 )

sumstat_colnames <- c("molecular_trait_id", "chromosome", "position", 
                      "ref", "alt", "variant", "ma_samples", "maf", 
                      "pvalue", "beta", "se", "type", "ac", "an", 
                      "r2", "molecular_trait_object_id", "gene_id", 
                      "median_tpm", "rsid")

###############################################################

message(" ## Reading Txref GTF file")
gtf_ref_txrev <- rtracklayer::import(txrev_gtf_file_path, 
                               colnames = c("ID", "gene_id", "Parent"),
                               feature.type = c("exon"))

message(" ## Reading GTF file for MANE")
gtf_ref <- rtracklayer::import(gtf_file_path, 
                               colnames = c("type", "gene_id", "gene_name", "gene_biotype", 
                                            "transcript_id", "transcript_name","transcript_biotype", "exon_number", "exon_id", "ccds_id"),
                               feature.type = c("exon"))

message(" ## Reading sample_metadata file")
sample_metadata <- readr::read_tsv(sample_meta_path)

message(" ## Reading mane_transcript_gene_map file")
mane_transcript_gene_map <- readr::read_tsv(mane_transcript_gene_map_file)

message(" ## Reading susie_purity_filtered file")
susie_purity_filtered <- readr::read_tsv(file = susie_file_path, col_types = "cccicccccdddddddd")

message(" ## Reading normalised usage matrix")
norm_exp_df <- readr::read_tsv(norm_usage_matrix_path)

message(" ## Reading leafcutter metadata file")
phenotype_metadata <- readr::read_tsv(phenotype_meta_path, col_types = "cccccddiccid") 

if (is.null(study_name)) { 
  assertthat::has_name(sample_metadata, "study" )
  study_name <- sample_metadata$study[1] 
}

if(assertthat::assert_that(all(!is.na(phenotype_metadata$gene_id) && all(!is.na(phenotype_metadata$gene_name))), 
                           msg = "All gene_id's and gene_name's in phenotype_metadata should be non-NA")) {
  message("All the phenotypes in phenotype_metadata has properly assigned to a certain gene!")
}

susie_high_pip_with_gene <- susie_purity_filtered %>% 
  dplyr::left_join(phenotype_metadata %>% dplyr::select(-chromosome), by = c("molecular_trait_id" = "phenotype_id")) 

highest_pip_vars_per_cs <- susie_high_pip_with_gene %>% 
  dplyr::group_by(cs_id) %>% 
  dplyr::arrange(-pip) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup()

if (debug_mode) {
  message(" ## Slicing 10 highest pip credible set variants for debug_mode")
  highest_pip_vars_per_cs <- highest_pip_vars_per_cs %>% 
    dplyr::arrange(-pip) %>% 
    dplyr::slice_head(n = 10)
}

message(" ## Will plot ", nrow(highest_pip_vars_per_cs), " highest pip per credible set signals.")
message(" ## Starting to plot")
for (index in 1:nrow(highest_pip_vars_per_cs)) {
  ss_oi = highest_pip_vars_per_cs[index,]
  message("index: ", index, ", txrev_phenotype_id: ", ss_oi$molecular_trait_id, ", variant: ", ss_oi$variant, ", gene_id: ", ss_oi$gene_id)
  
  # get all the intons in leafcutter cluster
  txrev_quant_pheno_ids <- phenotype_metadata %>% dplyr::filter(quant_id %in% ss_oi$quant_id)
  txrev_exons <- make_txrev_exon_granges(gff = gtf_ref_txrev, txrev_pheno_ids = txrev_quant_pheno_ids$phenotype_id)
  txrev_exons_cdss <- make_txrev_exon_granges(gff = gtf_ref_txrev, ss_oi$molecular_trait_id)
  
  exons_to_plot = txrev_exons
  exon_cdss_to_plot = txrev_exons_cdss
  
  variant_regions_vcf <- ss_oi %>% 
    dplyr::select(variant, chromosome, position) %>% 
    dplyr::mutate(region = paste0(chromosome,":", position, "-", position))
  snps <- seqminer::tabix.read.table(vcf_file_path, variant_regions_vcf$region)
  snps_filt <- snps %>% 
    dplyr::filter(ID %in% variant_regions_vcf$variant) %>% 
    dplyr::arrange(CHROM, POS)
  
  var_genotype <- snps_filt %>% 
    dplyr::select(-c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")) %>% 
    base::t() %>% 
    BiocGenerics::as.data.frame() %>% 
    dplyr::rename("GT_DS" = "V1") %>% 
    dplyr::mutate(GT = gsub(pattern = "\\:.*", replacement = "", x = GT_DS)) %>% 
    dplyr::mutate(REF = gsub(pattern = "\\|.*", replacement = "", x = GT)) %>% 
    dplyr::mutate(ALT = gsub(pattern = ".*\\|", replacement = "", x = GT)) %>% 
    dplyr::mutate(DS = as.numeric(REF) + as.numeric(ALT)) %>% 
    dplyr::mutate(genotype_id = BiocGenerics::rownames(.)) %>% 
    dplyr::mutate(genotype_id = gsub(pattern = "\\.", replacement = "-", x = genotype_id)) %>%
    left_join(sample_metadata %>% dplyr::select(sample_id, genotype_id, qtl_group), by = "genotype_id")
  
  bigwig_files <- list.files(bigwig_files_path, full.names = T)
  
  track_data_study <- data.frame(file_name = (gsub(pattern = paste0(bigwig_files_path, "/"), replacement = "", x = bigwig_files)),
                                  bigWig = bigwig_files)
  track_data_study <- track_data_study %>% dplyr::mutate(sample_id = gsub(pattern = ".bigwig", replacement = "", x = file_name))
  
  track_data_study <-track_data_study %>% 
    dplyr::left_join(var_genotype %>% dplyr::select(sample_id, qtl_group, DS), by = c("sample_id")) %>% 
    dplyr::mutate(scaling_factor = 1, track_id = qtl_group) %>% 
    dplyr::mutate(colour_group = as.character(DS)) %>% 
    dplyr::select(sample_id, scaling_factor, bigWig, track_id, colour_group, qtl_group)
  
  # Generate the output path 
  signal_name <- paste0(gsub(pattern = ":", replacement = "_", x = ss_oi$molecular_trait_id), "&", ss_oi$variant, "&", ss_oi$gene_id)
  path_plt = file.path(output_dir, signal_name)
  if (!dir.exists(path_plt)){
    dir.create(path_plt, recursive = TRUE)
  }
  
  if(!is.null(nominal_exon_sumstats_path)) {
    message(" ## Reading exon summary statistics")
    nom_exon_cc_sumstats <- seqminer::tabix.read.table(nominal_exon_sumstats_path, variant_regions_vcf$region) 
    colnames(nom_exon_cc_sumstats) <- sumstat_colnames
  
    # Extract the QTLs of exons according to gene and variant of interest
    nom_exon_cc_sumstats_filt <- nom_exon_cc_sumstats %>% 
      dplyr::filter(variant == ss_oi$variant, molecular_trait_object_id == ss_oi$gene_id) %>% 
      dplyr::filter(rsid == rsid[1]) %>% # if variant has more than one rsid keep only the first unique rsid 
      dplyr::mutate(exon_end = as.numeric(gsub(pattern = ".*\\_", replacement = "", x = molecular_trait_id))) %>% 
      dplyr::mutate(exon_start = gsub(pattern = "_[^_]+$", replacement = "", x = molecular_trait_id)) %>% 
      dplyr::mutate(exon_start = as.numeric(gsub(pattern = ".*\\_", replacement = "", x = exon_start))) %>% 
      dplyr::mutate(exon_center = round((exon_start + exon_end) / 2)) %>% 
      dplyr::mutate(exon_length = abs(exon_start - exon_end) + 1) %>% 
      dplyr::mutate(exon_row_num = dplyr::row_number()) %>% 
      dplyr::mutate(se_top = beta + se) %>% 
      dplyr::mutate(se_bottom = beta - se) %>% 
      dplyr::mutate(interval = ci.value * se) %>% 
      dplyr::mutate(p_fdr = p.adjust(pvalue, method = "fdr"))
    
    if (nrow(nom_exon_cc_sumstats_filt) > 0) {
      nom_exon_granges <- list()
      nom_exon_granges[[paste0("GENE:",ss_oi$gene_name)]] = GenomicRanges::GRanges(
        seqnames = nom_exon_cc_sumstats_filt$chromosome,
        ranges = IRanges::IRanges(start = nom_exon_cc_sumstats_filt$exon_start, end = nom_exon_cc_sumstats_filt$exon_end),
        strand = ifelse(test = ss_oi$strand == 1, yes = "+", no = "-"),
        mcols = data.frame(exon_id = nom_exon_cc_sumstats_filt$molecular_trait_id, 
                           gene_id = nom_exon_cc_sumstats_filt$gene_id))
      
      exons_to_plot <- append(exons_to_plot, nom_exon_granges)
    }
  }
  
  if (!is.null(mane_transcript_gene_map_file)) {
    MANE_transcript_oi <- mane_transcript_gene_map %>% 
      dplyr::filter(gene_id %in% ss_oi$gene_id) %>% 
      dplyr::pull(transcript_id)

    if (length(MANE_transcript_oi) > 0) {
      mane_transcript_exons <-  make_transcript_exon_granges(gff = gtf_ref, transcript_ids = MANE_transcript_oi)
      mane_transcript_exons_cdss <-  make_transcript_exon_granges_ccds(gff = gtf_ref, transcript_ids = MANE_transcript_oi)

      # mane_transcript_exons_df <- mane_transcript_exons[[1]] %>% BiocGenerics::as.data.frame()
      exons_to_plot <- append(exons_to_plot, mane_transcript_exons)
      exon_cdss_to_plot <- append(exon_cdss_to_plot, mane_transcript_exons_cdss)
    }
  }
  
  plot_rel_height = ifelse(length(txrev_exons)+1 <= 5, 3, length(txrev_exons)) 
  plot_rel_height = ifelse(plot_rel_height > 20, 20, plot_rel_height)

  exon_plot_data <- wiggleplotr::generateTxStructurePlotData(exons = exons_to_plot,
                                                             cdss = exon_cdss_to_plot)
  
  exon_plot <- wiggleplotr::plotTranscriptStructure(exons_df = exon_plot_data$transcript_struct_df, limits = exon_plot_data$limits)
  
  coverage_plot_data = wiggleplotr::generateCoveragePlotData(exons = exons_to_plot, 
                                                             cdss = exon_cdss_to_plot, 
                                                             plot_fraction = 0.2,
                                                             track_data = track_data_study %>% dplyr::filter(qtl_group==qtl_group_in))
  
  coverage_plot_data$coverage_df <- coverage_plot_data$coverage_df %>% dplyr::filter(!is.na(coverage))
  coverage_plot = wiggleplotr::makeCoveragePlot(coverage_df = coverage_plot_data$coverage_df, 
                                                limits = coverage_plot_data$limits, 
                                                alpha = 1, 
                                                fill_palette = wiggleplotr::getGenotypePalette(), 
                                                coverage_type = "line", 
                                                show_genotype_legend = TRUE)

  # permute the rows so that it becomes anonymous
  coverage_plot_data$coverage_df <- coverage_plot_data$coverage_df[sample(nrow(coverage_plot_data$coverage_df)),]
  
  if (nrow(nom_exon_cc_sumstats_filt) > 0) {
    beta_plot <- generate_beta_plot(transcript_struct_df_loc = exon_plot_data$transcript_struct_df, 
                                    nom_exon_cc_sumstats_filt_loc = nom_exon_cc_sumstats_filt,
                                    limits = exon_plot_data$limits)
    merged_plot <- cowplot::plot_grid(coverage_plot, beta_plot, exon_plot , align = "v", axis = "lr", rel_heights = c(3, 3, plot_rel_height), ncol = 1)
  } else {
    merged_plot <- cowplot::plot_grid(coverage_plot, exon_plot , align = "v", axis = "lr", rel_heights = c(3, plot_rel_height), ncol = 1)
  }
  
  filename_plt = paste0("cov_plot_", signal_name,".pdf")
  ggplot2::ggsave(path = path_plt, filename = filename_plt, plot = merged_plot, device = "pdf", width = 10, height = 8)
  
  # BOXPLOTS START HERE
  norm_exp_df_oi <- norm_exp_df %>% dplyr::filter(phenotype_id %in% txrev_quant_pheno_ids$phenotype_id)
  norm_exp_df_oi <- tibble::column_to_rownames(.data = norm_exp_df_oi,var = "phenotype_id")
  norm_exp_df_oi <- norm_exp_df_oi %>% base::t() %>% 
    GenomicRanges::as.data.frame() %>% 
    tibble::rownames_to_column(var = "sample_id")
  norm_exp_df_oi <- norm_exp_df_oi %>% 
    tidyr::pivot_longer(cols = -sample_id, names_to="intron_id", values_to = "norm_exp")
  
  track_data_study_box <- track_data_study %>% 
    dplyr::mutate(genotype_text = as.factor(colour_group)) %>% 
    dplyr::filter(qtl_group %in% qtl_group_in) %>%
    dplyr::mutate(condition_name = qtl_group) %>% 
    dplyr::mutate(gene_name = ss_oi$gene_name) %>% 
    dplyr::mutate(snp_id = ss_oi$variant) 
  
  track_data_study_box <- norm_exp_df_oi %>%  
    dplyr::left_join(track_data_study_box, by = "sample_id") %>% 
    dplyr::mutate(is_significant = intron_id == ss_oi$molecular_trait_id)
  
  nom_cc_sumstats <- seqminer::tabix.read.table(nominal_sumstats_path, variant_regions_vcf$region) 
  colnames(nom_cc_sumstats) <- sumstat_colnames
  nom_cc_sumstats <- nom_cc_sumstats %>% 
    dplyr::filter(variant %in% variant_regions_vcf$variant)
  
  # Keep only 1 rsid per variant per molecular_trait_id
  nom_cc_sumstats <- nom_cc_sumstats %>% 
    dplyr::group_by(molecular_trait_id, variant) %>% 
    dplyr::filter(rsid == rsid[1]) %>% 
    dplyr::ungroup()
  
  nom_cc_sumstats_filt <- nom_cc_sumstats %>% 
    dplyr::select(molecular_trait_id, pvalue, beta, se, maf) %>% 
    dplyr::rename(intron_id = molecular_trait_id)
  
  track_data_study_box_wrap <- track_data_study_box %>% 
    dplyr::left_join(nom_cc_sumstats_filt, by = "intron_id") %>% 
    dplyr::mutate(stats_text = paste0("Pval: ", pvalue, "			BETA: ", beta, 
                                      "\nSE: ", se)) %>% 
    dplyr::mutate(intron_id_with_stats = paste0(intron_id, "\n", stats_text))
  
  message(" ## Plotting box plots")

  boxplot_facet <- ggplot2::ggplot(track_data_study_box_wrap, 
                                   ggplot2::aes(x = genotype_text, 
                                                y = norm_exp, 
                                                color = is_significant, 
                                                group = genotype_text)) + 
    ggplot2::facet_wrap(~ intron_id_with_stats, dir = "v") + 
    ggplot2::geom_boxplot(outlier.shape = NA) + 
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
    ggplot2::ylab(paste0("Usage of txrev events ", ss_oi$group_id, "\n(", ss_oi$gene_name, " gene region)")) +
    ggplot2::xlab(paste0(track_data_study_box_wrap$snp_id[1], "    (MAF:", track_data_study_box_wrap$maf[1], ")")) + 
    ggplot2::theme_light() + 
    ggplot2::theme(strip.text.x = ggplot2::element_text(colour = "grey10"), strip.background = ggplot2::element_rect(fill = "grey85"))
  
  filename_plt_box_facet = paste0("box_facet_plot_", signal_name,".pdf")
  ggplot2::ggsave(path = path_plt, filename = filename_plt_box_facet, plot = boxplot_facet, device = "pdf", width = 10, height = 11)
  
  track_data_study_box_wrap_for_RDS <- track_data_study_box_wrap %>%
    dplyr::select(genotype_text, norm_exp, is_significant, intron_id, pvalue, beta, se, snp_id, maf)
  
  track_data_study_box_wrap_for_RDS <- track_data_study_box_wrap_for_RDS[sample(nrow(track_data_study_box_wrap_for_RDS)),]
  
  limit_max <- max(coverage_plot_data$limits, exon_plot_data$limits)
  
  tx_str_df <- exon_plot_data$transcript_struct_df %>% dplyr::mutate(limit_max = limit_max)
  Rds_list <- list(coverage_plot_df = coverage_plot_data$coverage_df, ss_oi = ss_oi)
  Rds_list[["tx_str_df"]] <- tx_str_df
  Rds_list[["box_plot_wrap_df"]] <- track_data_study_box_wrap_for_RDS
  Rds_plot_file_name <- paste0(path_plt, "/plot_data_", signal_name,".Rds")
  saveRDS(object = Rds_list, compress = "gzip", file = Rds_plot_file_name)
  
  tar_path = paste0(path_plt, "/plot_data_tsv/")
  if (!dir.exists(tar_path)){
    dir.create(tar_path, recursive = TRUE)
  }

  gzfile = gzfile(paste0(tar_path, "/coverage_df_", signal_name, ".tsv.gz"), "w")
  write.table(x = coverage_plot_data$coverage_df, file = gzfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  close(gzfile)
  
  gzfile = gzfile(paste0(tar_path, "/tx_str_", signal_name, ".tsv.gz"), "w")
  write.table(x = tx_str_df, file = gzfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  close(gzfile)

  gzfile = gzfile(paste0(tar_path, "/box_plot_df_", signal_name, ".tsv.gz"), "w")
  write.table(x = track_data_study_box_wrap_for_RDS, file = gzfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  close(gzfile)

  gzfile = gzfile(paste0(tar_path, "/ss_oi_df_", signal_name, ".tsv.gz"), "w")
  write.table(x = ss_oi, file = gzfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  close(gzfile)
  # write_tsv(x = coverage_plot_data$coverage_df, file = paste0(tar_path, "/coverage_df_", signal_name, ".tsv") )
  # write_tsv(x = tx_str_df, file = paste0(tar_path, "/tx_str_", signal_name, ".tsv") )
  # write_tsv(x = track_data_study_box_wrap_for_RDS, file = paste0(tar_path, "/box_plot_df_", signal_name, ".tsv") )
  # write_tsv(x = ss_oi, file = paste0(tar_path, "/ss_oi_df_", signal_name, ".tsv") )
  
  # signal_name <- gsub(pattern = "&", replacement = "\\&", x = signal_name)
  
  # filename_all_plt_data_tar = paste0(path_plt, "/plot_data_", signal_name,".tar.gz")
  # setwd(path_plt)
  # tar(tarfile = filename_all_plt_data_tar, files = "plot_data_tsv",
  #     compression = "gzip")
  # unlink("plot_data_tsv", recursive = TRUE)
  # setwd("../..")

  if (!individual_boxplots) {
    next
  }
  
  message(" ## Plotting individual boxplots")
  for (intron_id_sel in track_data_study_box$intron_id %>% BiocGenerics::unique()) {
    nom_cc_sumstats_filt <- nom_cc_sumstats %>% 
      dplyr::filter(molecular_trait_id == intron_id_sel, variant %in% ss_oi$variant) %>% 
      dplyr::slice(1)
    
    labels <- c(paste0("P-value: ", nom_cc_sumstats_filt$pvalue, "			BETA: ", nom_cc_sumstats_filt$beta, 
                       "\nStd. Err: ", nom_cc_sumstats_filt$se, "				MAF: ", nom_cc_sumstats_filt$maf))
    
    track_data_study_box_intron <- track_data_study_box %>% dplyr::filter(intron_id == intron_id_sel)
    
    box_plot <- ggplot2::ggplot(track_data_study_box_intron, 
                                ggplot2::aes(x = genotype_text, 
                                             y = norm_exp, 
                                             color = condition_name, 
                                             group = genotype_text)) + 
      ggplot2::geom_boxplot(outlier.shape = NA) + 
      ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
      ggplot2::ylab(paste0(track_data_study_box_intron$intron_id[1], " usage")) +
      ggplot2::xlab(track_data_study_box_intron$snp_id[1]) + 
      ggplot2::theme_light() + 
      ggplot2::labs(subtitle = labels) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(colour = "grey10"), strip.background = ggplot2::element_rect(fill = "grey85"))
    
    filename_plt_box = paste0("box_plot_", gsub(pattern = ":", replacement = "_", x = intron_id_sel), "&", ss_oi$variant, "&", ss_oi$gene_id,".pdf")
    ggplot2::ggsave(path = path_plt, filename = filename_plt_box, plot = box_plot, device = "pdf", width = 7, height = 5)
  }
}