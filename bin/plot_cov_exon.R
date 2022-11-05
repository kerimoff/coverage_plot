message(" ## Loading libraries: optparse")
suppressPackageStartupMessages(library("optparse"))

#Parse command-line options
option_list <- list(
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
  make_option(c("-u", "--usage_matrix_norm"), type="character", default=NULL,
              help="Path to the normalised usage matrix", metavar = "type"),
  make_option(c("--div_scaling_factors"), type="character", default=NULL,
              help="Path to the scaling_factors file", metavar = "type"),
  make_option(c("-a", "--all_summ_stats"), type="character", default=NULL,
              help="Path to the full gene nominal summary statistics", metavar = "type"),
  make_option(c("-e", "--exon_summ_stats"), type="character", default=NULL,
              help="Path to the full exon nominal summary statistics", metavar = "type"),
  make_option(c("-w", "--wiggle_plotr_path"), type="character", default=NULL,
              help="Local path to wiggleplotr package", metavar = "type"),
  make_option(c("-i", "--individual_boxplots"), type="logical", 
              action = "store_true", default=FALSE,
              help="Flag to generate individual boxplots", metavar = "type")
)

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, cqn, data.table")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("seqminer"))

make_transcript_exon_granges <- function(gff, transcript_ids, add_gene = FALSE) {
  exon_list <- list()
  for (transcript_id in transcript_ids) {
    transcript_exons_temp <- gff[(base::gsub("\\..*","",SummarizedExperiment::elementMetadata(gff)[,"transcript_id"]) == transcript_id)]
    gene_id = transcript_exons_temp$gene_name[1]
    if (add_gene) {
      exon_list[[paste0("GENE:", gene_id, ":", transcript_id)]] <- transcript_exons_temp
    } else {
      exon_list[[transcript_id]] <- transcript_exons_temp
    }
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

make_transcript_exon_granges_ccds <- function(gff, transcript_ids, add_gene = FALSE) {
  exon_list <- list()
  for (transcript_id in transcript_ids) {
    transcript_exons_temp <- gff[(base::gsub("\\..*","",SummarizedExperiment::elementMetadata(gff)[,"transcript_id"]) == transcript_id & !is.na(SummarizedExperiment::elementMetadata(gff)[,"ccds_id"]))]
    gene_id = transcript_exons_temp$gene_name[1]
    if (add_gene) {
      exon_list[[paste0("GENE:", gene_id, ":", transcript_id)]] <- transcript_exons_temp
    } else {
      exon_list[[transcript_id]] <- transcript_exons_temp
    }
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

generate_beta_plot <- function(transcript_struct_df_loc, 
                               nom_exon_cc_sumstats_filt_loc, 
                               limits, 
                               ss_oi,
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

prepareTranscriptStructureForPlotting <- function(exon_ranges, cds_ranges, transcript_annotations){
  #Combine exon_ranges and cds_ranges into a single data.frame that also contains transcript rank
  
  #Convert exon ranges into data.frame and add transcript rank
  exons_df = purrr::map_df(exon_ranges, data.frame, .id = "transcript_id")
  exons_df = dplyr::mutate(exons_df, transcript_rank = as.numeric(factor(exons_df$transcript_id)), type = "")
  transcript_rank = unique(exons_df[,c("transcript_id", "transcript_rank", "type")])
  
  #Convert CDS ranges into a data.frame
  cds_df = purrr::map_df(cds_ranges, data.frame, .id = "transcript_id")
  cds_df = dplyr::left_join(cds_df, transcript_rank, by = "transcript_id") #Add matching transcript rank
  
  #Join exons and cdss together
  exons_df = dplyr::mutate(exons_df, feature_type = "exon")
  cds_df = dplyr::mutate(cds_df, feature_type = "cds")
  transcript_struct = rbind(exons_df, cds_df)
  
  #Add transcript label to transcript structure
  transcript_struct = dplyr::left_join(transcript_struct, transcript_annotations, by = "transcript_id")
  return(transcript_struct)
}

#Debugging
if (FALSE) {
  opt = list()
  opt$n = "Alasoo_2018"
  opt$f = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/single_batch_debug_mode.tsv"
  opt$s = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/Alasoo_2018.tsv"
  opt$p = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/exon_counts_Ensembl_105_phenotype_metadata.tsv.gz"
  opt$q = "macrophage_IFNg"
  opt$o = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/plots_data_generate_EXON"
  opt$v = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/Alasoo_2018.filtered.vcf.gz"
  opt$b = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/bigwig/"
  opt$m = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/MANE_transcript_gene_map.txt"
  opt$g = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/Homo_sapiens.GRCh38.105.gtf"
  opt$div_scaling_factors = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/Alasoo_2018.macrophage_IFNg.scaling_factors.tsv.gz"
  opt$u = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/Alasoo_2018.macrophage_IFNg.exon_counts_TPM_norm.tsv.gz"
  opt$e = "/Users/kerimov/Work/temp_files/debug/6e27637dec31f4eaa6230b3f23a9a7/debug/nom_exon_cc_sumstats_all.tsv"
  opt$w = "/Users/kerimov/Work/R_env/wiggleplotr/"
  opt$individual_boxplots = TRUE
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
norm_usage_matrix_path = opt$u
scaling_factors_path = opt$div_scaling_factors
nominal_sumstats_path = opt$a
nominal_exon_sumstats_path = opt$e
wiggleplotr_path = opt$w
individual_boxplots = opt$individual_boxplots


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
message("######### scaling_fct_path   : ", scaling_factors_path)
message("######### norm_usage_matrix  : ", norm_usage_matrix_path)
message("######### nominal_sumstats   : ", nominal_sumstats_path)
message("######### exon_sumstats      : ", nominal_exon_sumstats_path)
message("######### wiggleplotr_path   : ", wiggleplotr_path)
message("######### individual_boxplots: ", individual_boxplots)

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

time_here <- function(prev_time, message_text = "Time in this point: "){
  message(message_text, Sys.time() - prev_time)
  return(Sys.time())
}

###############################################################
start_time <- Sys.time()

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
highest_pip_vars_per_cs <- readr::read_tsv(file = susie_file_path, col_types = "cccicccccddddddddcccddccccd")

message(" ## Reading normalised usage matrix")
norm_exp_df <- readr::read_tsv(norm_usage_matrix_path)

message(" ## Reading metadata file")
phenotype_metadata <- readr::read_tsv(phenotype_meta_path, col_types = "cccccddicciddi") 

message(" ## Reading scaling_factors file")
scaling_factor_data <- readr::read_tsv(scaling_factors_path, col_types = "cd") 

start_time <- time_here(prev_time = start_time, message_text = " >> Read input TSVs in: ")
if (is.null(study_name)) { 
  assertthat::has_name(sample_metadata, "study" )
  study_name <- sample_metadata$study[1] 
}

if(assertthat::assert_that(all(!is.na(phenotype_metadata$gene_id) && all(!is.na(phenotype_metadata$gene_name))), 
                           msg = "All gene_id's and gene_name's in phenotype_metadata should be non-NA")) {
  message("All the phenotypes in phenotype_metadata has properly assigned to a certain gene!")
}

# susie_high_pip_with_gene <- susie_purity_filtered %>% 
#   dplyr::left_join(phenotype_metadata %>% dplyr::select(-chromosome), by = c("molecular_trait_id" = "phenotype_id")) 

# # Get highest PIP variant for each credible set
# highest_pip_vars_per_cs <- susie_high_pip_with_gene %>% 
#   dplyr::group_by(cs_id) %>% 
#   dplyr::arrange(-pip) %>% 
#   dplyr::slice(1) %>% 
#   dplyr::ungroup()

# assertthat::assert_that(!is.null(nominal_exon_sumstats_path), msg = "ERROR: nominal_exon_sumstats_path is null!")

variant_regions_vcf <- highest_pip_vars_per_cs %>% 
  dplyr::select(variant, chromosome, position) %>% 
  dplyr::mutate(region = paste0(chromosome,":", position, "-", position))

message(" ## Reading exon summary statistics")
nom_exon_cc_sumstats_all <- seqminer::tabix.read.table(nominal_exon_sumstats_path, variant_regions_vcf$region) 
colnames(nom_exon_cc_sumstats_all) <- sumstat_colnames
if (is.null(nom_exon_cc_sumstats_all) || nrow(nom_exon_cc_sumstats_all) == 0) {
  message("Exiting! Weirdly there are no exon summary statistics for this variant: ", variant_regions_vcf$region)
  if (!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
}

if (nrow(highest_pip_vars_per_cs) == 10) {
  dir.create("debug")
  readr::write_tsv(nom_exon_cc_sumstats_all, "debug/nom_exon_cc_sumstats_all.tsv")
}
message(" ## Reading exon summary statistics complete")

message(" ## Reading all variants from VCF_file")
snps_all <- seqminer::tabix.read.table(vcf_file_path, variant_regions_vcf$region)
if (study_name == "Steinberg_2020") {
  names(snps_all) <- gsub(pattern = ".", replacement = ":", x = names(snps_all), fixed = T)
}
if (study_name == "Quach_2016") {
  names(snps_all) <- gsub(pattern = ".", replacement = "@", x = names(snps_all), fixed = T)
}
if (study_name %in% c("Schmiedel_2018", "Bossini-Castillo_2019", "ROSMAP", "iPSCORE")) {
  names(snps_all) <- gsub(pattern = "X", replacement = "", x = names(snps_all), fixed = T)  
}
if (study_name %in% c("iPSCORE")) {
  names(snps_all) <- gsub(pattern = ".", replacement = "-", x = names(snps_all), fixed = T)  
}
message(" ## Reading all variants from VCF_file complete")

start_time <- time_here(prev_time = start_time, message_text = " >> seqminer tabix took: ")

message(" ## Will plot batch of ", nrow(highest_pip_vars_per_cs), " highest pip per credible set signals.")
message(" ## Starting to plot")


for (index in 1:nrow(highest_pip_vars_per_cs)) {
  ss_oi = highest_pip_vars_per_cs[index,]
  message("index: ", index, ", exon_phenotype_id: ", ss_oi$molecular_trait_id, ", variant: ", ss_oi$variant)
  
  message(" ## Preparing Track Data")
  
  # Get all exons of the gene (with group_id)
  exon_quant_pheno_ids <- phenotype_metadata %>% dplyr::filter(group_id %in% ss_oi$group_id)

  variant_regions_vcf <- ss_oi %>% 
    dplyr::select(variant, chromosome, position) %>% 
    dplyr::mutate(region = paste0(chromosome,":", position, "-", position))

  snps_filt <- snps_all %>% 
    dplyr::filter(ID %in% variant_regions_vcf$variant) %>% 
    dplyr::arrange(CHROM, POS)
  
  # var_genotype <- snps_filt %>% 
  #   dplyr::select(-c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")) %>% 
  #   base::t() %>% 
  #   BiocGenerics::as.data.frame() %>% 
  #   dplyr::rename("GT_DS" = "V1") %>% 
  #   dplyr::mutate(GT = gsub(pattern = "\\:.*", replacement = "", x = GT_DS)) %>% 
  #   dplyr::mutate(REF = gsub(pattern = "\\|.*", replacement = "", x = GT)) %>% 
  #   dplyr::mutate(ALT = gsub(pattern = ".*\\|", replacement = "", x = GT)) %>% 
  #   dplyr::mutate(DS = as.numeric(REF) + as.numeric(ALT)) %>% 
  #   dplyr::mutate(genotype_id = BiocGenerics::rownames(.)) %>% 
  #   dplyr::mutate(genotype_id = gsub(pattern = "\\.", replacement = "-", x = genotype_id)) %>%
  #   left_join(sample_metadata %>% dplyr::select(sample_id, genotype_id, qtl_group), by = "genotype_id")
  
if (study_name == "Lepik_2017") {
    var_genotype <- snps_filt %>% 
      dplyr::select(-c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")) %>% 
      base::t() %>% 
      BiocGenerics::as.data.frame() %>% 
      dplyr::rename("GT_DS" = "V1") %>% 
      dplyr::mutate(GT = gsub(pattern = "\\:.*", replacement = "", x = GT_DS)) %>% 
      dplyr::mutate(REF = gsub(pattern = "\\/.*", replacement = "", x = GT)) %>% 
      dplyr::mutate(ALT = gsub(pattern = ".*\\/", replacement = "", x = GT)) %>% 
      dplyr::mutate(DS = as.numeric(REF) + as.numeric(ALT)) %>% 
      dplyr::mutate(genotype_id = BiocGenerics::rownames(.)) %>% 
      dplyr::mutate(genotype_id = gsub(pattern = "\\.", replacement = "-", x = genotype_id))
  } else  {
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
      dplyr::mutate(genotype_id = gsub(pattern = "\\.", replacement = "-", x = genotype_id))
  }

  sample_meta_clean = sample_metadata %>% 
    dplyr::filter(rna_qc_passed, genotype_qc_passed) %>%  
    dplyr::select(sample_id, genotype_id, qtl_group)
  
  var_genotype <- sample_meta_clean %>% 
    left_join(var_genotype, by = "genotype_id")

  bigwig_files <- list.files(bigwig_files_path, full.names = T)
  
  track_data_study <- data.frame(file_name = (gsub(pattern = paste0(bigwig_files_path, "/"), replacement = "", x = bigwig_files, fixed = T)),
                                  bigWig = bigwig_files)
  
  track_data_study <- track_data_study %>% 
    dplyr::mutate(sample_id = gsub(pattern = ".bigwig", replacement = "", x = file_name)) %>% 
    dplyr::filter(sample_id %in% sample_meta_clean$sample_id)
    
  track_data_study <-track_data_study %>% 
    dplyr::left_join(var_genotype %>% dplyr::select(sample_id, qtl_group, DS), by = c("sample_id")) %>% 
    dplyr::left_join(scaling_factor_data) %>% 
    dplyr::rename(scaling_factor = scaling_factors) %>% 
    dplyr::mutate(track_id = qtl_group) %>% 
    dplyr::mutate(colour_group = as.character(DS)) %>% 
    dplyr::select(sample_id, scaling_factor, bigWig, track_id, colour_group, qtl_group) %>% 
    dplyr::filter(qtl_group==qtl_group_in)
  
  # Generate the output path 
  signal_name <- paste0(gsub(pattern = ":", replacement = "_", x = ss_oi$molecular_trait_id), "___", ss_oi$variant)
  path_plt = file.path(output_dir, signal_name)
  if (!dir.exists(path_plt)){
    dir.create(path_plt, recursive = TRUE)
  }


  # Extract the QTLs of exons according to gene and variant of interest
  nom_exon_cc_sumstats_filt <- nom_exon_cc_sumstats_all %>% 
    dplyr::filter(variant == ss_oi$variant, molecular_trait_object_id == ss_oi$gene_id) %>% 
    dplyr::filter(rsid == rsid[1]) %>% # if variant has more than one rsid keep only the first unique rsid 
    dplyr::mutate(exon_end = as.numeric(gsub(pattern = ".*\\_", replacement = "", x = molecular_trait_id))) %>% 
    dplyr::mutate(exon_start = gsub(pattern = "_[^_]+$", replacement = "", x = molecular_trait_id)) %>% 
    dplyr::mutate(exon_start = as.numeric(gsub(pattern = ".*\\_", replacement = "", x = exon_start))) %>% 
    dplyr::mutate(exon_row_num = dplyr::row_number()) %>% 
    dplyr::mutate(interval = ci.value * se) %>% 
    dplyr::mutate(p_fdr = p.adjust(pvalue, method = "fdr"))
  
  if (nrow(nom_exon_cc_sumstats_filt) == 0) {
    message("Weirdly there are no exon summary statistics for this credible set: GENE_ID:", ss_oi$gene_id, ", Variant:", ss_oi$variant )
    next
  }
  
  nom_exon_granges <- list()
  nom_exon_granges_cc <- list()
  nom_exon_granges[[paste0("GENE:",ss_oi$gene_name)]] = GenomicRanges::GRanges(
    seqnames = nom_exon_cc_sumstats_filt$chromosome,
    ranges = IRanges::IRanges(start = nom_exon_cc_sumstats_filt$exon_start, end = nom_exon_cc_sumstats_filt$exon_end),
    strand = ifelse(test = ss_oi$strand == 1, yes = "+", no = "-"),
    mcols = data.frame(exon_id = nom_exon_cc_sumstats_filt$molecular_trait_id, 
                       gene_id = nom_exon_cc_sumstats_filt$gene_id))
  
  
  nom_exon_cc_sumstats_ss_oi = nom_exon_cc_sumstats_filt %>% 
    dplyr::filter(molecular_trait_id == ss_oi$molecular_trait_id)
  
  nom_exon_granges_cc[[paste0("GENE:",ss_oi$gene_name)]] = GenomicRanges::GRanges(
    seqnames = nom_exon_cc_sumstats_ss_oi$chromosome,
    ranges = IRanges::IRanges(start = nom_exon_cc_sumstats_ss_oi$exon_start, end = nom_exon_cc_sumstats_ss_oi$exon_end),
    strand = ifelse(test = ss_oi$strand == 1, yes = "+", no = "-"),
    mcols = data.frame(exon_id = nom_exon_cc_sumstats_ss_oi$molecular_trait_id, 
                       gene_id = nom_exon_cc_sumstats_ss_oi$gene_id))
  
  exons_to_plot <- nom_exon_granges
  exon_cdss_to_plot <- nom_exon_granges_cc
  
  
  if (!is.null(mane_transcript_gene_map_file)) {
    MANE_transcript_oi <- mane_transcript_gene_map %>% dplyr::filter(gene_id %in% ss_oi$gene_id) %>% dplyr::pull(transcript_id)
    mane_transcript_exons <-  make_transcript_exon_granges(gff = gtf_ref, transcript_ids = MANE_transcript_oi, add_gene = TRUE)
    # mane_transcript_exons_cdss <-  make_transcript_exon_granges_ccds(gff = gtf_ref, transcript_ids = MANE_transcript_oi, add_gene = TRUE)

    # mane_transcript_exons_df <- mane_transcript_exons[[1]] %>% BiocGenerics::as.data.frame()
    exons_to_plot <- append(exons_to_plot, mane_transcript_exons)
    exon_cdss_to_plot <- append(exon_cdss_to_plot, mane_transcript_exons)
  }
  
  plot_rel_height = 3
  
  message(" ## Extracting coverage data")
  coverage_data_list = tryCatch(wiggleplotr::extractCoverageData(exons = exons_to_plot, 
                                                        cdss = exon_cdss_to_plot, 
                                                        plot_fraction = 0.2,
                                                        track_data = track_data_study), 
  error = function(e) {
    message(" ## Problem with generating coverage_data wiggleplotr")
    message(e)
  })

  if (!exists("coverage_data_list")) {
    message(" ERROR: !exists")
    next
  }
  if (all(is.na(coverage_data_list)) | length(coverage_data_list) == 0) {
    message(" ERROR: is.na(coverage_data_list) | length(coverage_data_list)")
    next
  }

  message(" ## Prepare transcript structure for plotting")
  tx_structure_df = prepareTranscriptStructureForPlotting(exon_ranges = coverage_data_list$tx_annotations$exon_ranges, 
                                                          cds_ranges = coverage_data_list$tx_annotations$cds_ranges, 
                                                          transcript_annotations = coverage_data_list$plotting_annotations)
  
  message(" ## Generate coverage data plots (wiggleplotr)")
  wiggle_plots <- wiggleplotr::plotCoverageData(coverage_data_list, alpha = 1, 
                                                fill_palette = wiggleplotr::getGenotypePalette(), 
                                                coverage_type = "line", return_subplots_list = TRUE, 
                                                show_legend = TRUE)
    
  intron_ss_oi_vert_lines = tx_structure_df %>% 
    dplyr::filter(transcript_id == ss_oi$molecular_trait_id, feature_type == "exon") 
  intron_ss_oi_vert_lines <- c(intron_ss_oi_vert_lines[1,] %>% dplyr::pull(end), 
                               intron_ss_oi_vert_lines[2,] %>% dplyr::pull(start))

  coverage_plot <- wiggle_plots$coverage_plot + ggplot2::geom_vline(xintercept = intron_ss_oi_vert_lines, alpha = 0.5, color = "lightgrey")
  exon_plot <- wiggle_plots$tx_structure + ggplot2::geom_vline(xintercept = intron_ss_oi_vert_lines, alpha = 0.5, color = "lightgrey")
  if (nrow(nom_exon_cc_sumstats_filt) > 0) {
    beta_plot <- generate_beta_plot(transcript_struct_df_loc = tx_structure_df, 
                                    nom_exon_cc_sumstats_filt_loc = nom_exon_cc_sumstats_filt,
                                    limits = coverage_data_list$limits,
                                    ss_oi = ss_oi)
    beta_plot <- beta_plot + ggplot2::geom_vline(xintercept = intron_ss_oi_vert_lines, alpha = 0.5, color = "lightgrey")
    merged_plot <- cowplot::plot_grid(coverage_plot, beta_plot, exon_plot , align = "v", axis = "lr", rel_heights = c(3, 3, plot_rel_height), ncol = 1)
  } else {
    merged_plot <- cowplot::plot_grid(coverage_plot, exon_plot , align = "v", axis = "lr", rel_heights = c(3, plot_rel_height), ncol = 1)
  }
  
  filename_plt = paste0("cov_plot_", signal_name,".pdf")
  ggplot2::ggsave(path = path_plt, filename = filename_plt, plot = merged_plot, device = "pdf", width = 10, height = 8)
  message(" ## Saved: ", filename_plt)

  # permute the rows so that it becomes anonymous
  coverage_data_list$coverage_df <- coverage_data_list$coverage_df[sample(nrow(coverage_data_list$coverage_df)),]
  
  message(" ## Prepare box plot data")
  # BOXPLOTS START HERE
  norm_exp_df_oi <- norm_exp_df %>% dplyr::filter(phenotype_id %in% exon_quant_pheno_ids$phenotype_id)
  norm_exp_df_oi <- tibble::column_to_rownames(.data = norm_exp_df_oi,var = "phenotype_id")
  norm_exp_df_oi <- norm_exp_df_oi %>% base::t() %>% 
    GenomicRanges::as.data.frame() %>% 
    tibble::rownames_to_column(var = "sample_id")
  norm_exp_df_oi <- norm_exp_df_oi %>% 
    tidyr::pivot_longer(cols = -sample_id, names_to="tx_id", values_to = "norm_exp")
  
  track_data_study_box <- track_data_study %>% 
    dplyr::mutate(genotype_text = as.factor(colour_group)) %>% 
    dplyr::filter(qtl_group %in% qtl_group_in) %>%
    dplyr::mutate(condition_name = qtl_group) %>% 
    dplyr::mutate(gene_name = ss_oi$gene_name) %>% 
    dplyr::mutate(snp_id = ss_oi$variant) 
  
  track_data_study_box <- norm_exp_df_oi %>%  
    dplyr::left_join(track_data_study_box, by = "sample_id") %>% 
    dplyr::mutate(is_significant = tx_id == ss_oi$molecular_trait_id)
  
  message(" ## Reading nominal summary stats with seqminer")
  nom_cc_sumstats <- nom_exon_cc_sumstats_all %>% 
    dplyr::filter(variant %in% variant_regions_vcf$variant)
  
  # Keep only 1 rsid per variant per molecular_trait_id
  nom_cc_sumstats <- nom_cc_sumstats %>% 
    dplyr::group_by(molecular_trait_id, variant) %>% 
    dplyr::filter(rsid == rsid[1]) %>% 
    dplyr::ungroup()
  
  nom_cc_sumstats_filt <- nom_cc_sumstats %>% 
    dplyr::select(molecular_trait_id, pvalue, beta, se, maf) %>% 
    dplyr::rename(tx_id = molecular_trait_id)
  
  track_data_study_box_wrap <- track_data_study_box %>% 
    dplyr::left_join(nom_cc_sumstats_filt, by = "tx_id") %>% 
    dplyr::mutate(stats_text = paste0("Pval: ", pvalue, "			BETA: ", beta, 
                                      "\nSE: ", se)) %>% 
    dplyr::mutate(tx_id_with_stats = paste0(tx_id, "\n", stats_text))
  
  message(" ## Plotting box plots")

  boxplot_facet <- ggplot2::ggplot(track_data_study_box_wrap, 
                                   ggplot2::aes(x = genotype_text, 
                                                y = norm_exp, 
                                                color = is_significant, 
                                                group = genotype_text)) + 
    ggplot2::facet_wrap(~ tx_id_with_stats, dir = "v") + 
    ggplot2::geom_boxplot(outlier.shape = NA) + 
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
    ggplot2::ylab(paste0("Expression of exons of ", ss_oi$group_id, " gene (", ss_oi$gene_name, ")")) +
    ggplot2::xlab(paste0(track_data_study_box_wrap$snp_id[1], "    (MAF:", track_data_study_box_wrap$maf[1], ")")) + 
    ggplot2::theme_light() + 
    ggplot2::theme(strip.text.x = ggplot2::element_text(colour = "grey10"), strip.background = ggplot2::element_rect(fill = "grey85"))
  
  filename_plt_box_facet = paste0("box_facet_plot_", signal_name,".pdf")
  save_width_height = round(sqrt(length(exon_quant_pheno_ids$phenotype_id))*3)
  ggplot2::ggsave(path = path_plt, filename = filename_plt_box_facet, plot = boxplot_facet, device = "pdf", width = save_width_height, height = save_width_height)
  
  track_data_study_box_wrap_for_RDS <- track_data_study_box_wrap %>%
    dplyr::select(genotype_text, norm_exp, is_significant, tx_id, pvalue, beta, se, snp_id, maf)
  
  track_data_study_box_wrap_for_RDS <- track_data_study_box_wrap_for_RDS[sample(nrow(track_data_study_box_wrap_for_RDS)),]
  
  tx_str_df <- tx_structure_df %>% dplyr::mutate(limit_max = max(coverage_data_list$limits))
  Rds_list <- list(coverage_data_list = coverage_data_list, ss_oi = ss_oi)
  Rds_list[["nom_exon_cc_sumstats_df"]] <- nom_exon_cc_sumstats_filt
  Rds_list[["box_plot_wrap_df"]] <- track_data_study_box_wrap_for_RDS
  Rds_plot_file_name <- paste0(path_plt, "/plot_data_", signal_name,".Rds")
  saveRDS(object = Rds_list, compress = "gzip", file = Rds_plot_file_name)
  
  tar_path = "./plot_data_tsv"
  if (!dir.exists(tar_path)){
    dir.create(tar_path, recursive = TRUE)
  }

  gzfile = gzfile(paste0(tar_path, "/coverage_df_", signal_name, ".tsv.gz"), "w")
  write.table(x = coverage_data_list$coverage_df, file = gzfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
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

  nom_exon_cc_sumstats_filt <- nom_exon_cc_sumstats_filt %>% 
    mutate(rsid = stringr::str_trim(rsid))
  
  gzfile = gzfile(paste0(tar_path, "/nom_exon_cc_", signal_name, ".tsv.gz"), "w")
  write.table(x = nom_exon_cc_sumstats_filt, file = gzfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  close(gzfile)
  
  signal_name <- gsub(pattern = "&", replacement = "\\&", x = signal_name)
  
  filename_all_plt_data_tar = paste0(path_plt, "/plot_data_", signal_name,".tar.gz")
  # prev_wd <- getwd()
  # setwd(path_plt)
  tar(tarfile = filename_all_plt_data_tar, files = "plot_data_tsv",
      compression = "gzip")
  unlink("plot_data_tsv", recursive = TRUE)
  # setwd(prev_wd)
  
  if (!individual_boxplots) {
    next
  }
  
  for (tx_id_sel in track_data_study_box$tx_id %>% BiocGenerics::unique()) {
    nom_cc_sumstats_filt <- nom_cc_sumstats %>% 
      dplyr::filter(molecular_trait_id == tx_id_sel, variant %in% ss_oi$variant) %>% 
      dplyr::slice(1)
    
    labels <- c(paste0("P-value: ", nom_cc_sumstats_filt$pvalue, "			BETA: ", nom_cc_sumstats_filt$beta, 
                       "\nStd. Err: ", nom_cc_sumstats_filt$se, "				MAF: ", nom_cc_sumstats_filt$maf))
    
    track_data_study_box_intron <- track_data_study_box %>% dplyr::filter(tx_id == tx_id_sel)
    
    box_plot <- ggplot2::ggplot(track_data_study_box_intron, 
                                ggplot2::aes(x = genotype_text, 
                                             y = norm_exp, 
                                             color = condition_name, 
                                             group = genotype_text)) + 
      ggplot2::geom_boxplot(outlier.shape = NA) + 
      ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
      ggplot2::ylab(paste0(track_data_study_box_intron$tx_id[1], " usage")) +
      ggplot2::xlab(track_data_study_box_intron$snp_id[1]) + 
      ggplot2::theme_light() + 
      ggplot2::labs(subtitle = labels) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(colour = "grey10"), strip.background = ggplot2::element_rect(fill = "grey85"))
    
    filename_plt_box = paste0("box_plot_", gsub(pattern = ":", replacement = "_", x = tx_id_sel), "&", ss_oi$variant, "&", ss_oi$gene_id,".pdf")
    ggplot2::ggsave(path = path_plt, filename = filename_plt_box, plot = box_plot, device = "pdf", width = 7, height = 5)
  }
}
