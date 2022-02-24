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
  make_option(c("-q", "--qtl_group"), type="character", default="gene_counts",
              help="Quantification method. Possible values: gene_counts, leafcutter, txrevise, transcript_usage, exon_counts and HumanHT-12_V4 [default \"%default\"]", metavar = "type"),
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
              help="Path to the normalised usage matrix", metavar = "type")
)

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, cqn, data.table")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("readr"))
library(seqminer)

make_transcript_exon_granges <- function(gff, transcript_ids) {
  exon_list <- list()
  for (transcript_id in transcript_ids) {
    exon_list[[transcript_id]] <- gff[(base::gsub("\\..*","",SummarizedExperiment::elementMetadata(gff)[,"transcript_id"]) == transcript_id)]
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

make_pseudo_exons <- function(df_introns, ps_exon_len = 50){
  uniq_intron_ids <- df_introns %>% dplyr::pull(phenotype_id)
  introns_oi_exon1 <- df_introns %>% dplyr::mutate(ps_exon_start = intron_start - ps_exon_len, ps_exon_end = intron_start)
  introns_oi_exon2 <- df_introns %>% dplyr::mutate(ps_exon_start = intron_end, ps_exon_end = intron_end + ps_exon_len)
  df_introns <- BiocGenerics::rbind(introns_oi_exon1, introns_oi_exon2) %>% dplyr::arrange(phenotype_id)
  exon_list <- list()
  for (intron_id in uniq_intron_ids) {
    df_introns_sub <- df_introns %>% dplyr::filter(phenotype_id == intron_id)
    exon_list[[intron_id]] <- GenomicRanges::GRanges(
      seqnames = df_introns_sub$chromosome,
      ranges = IRanges::IRanges(start = df_introns_sub$ps_exon_start, end = df_introns_sub$ps_exon_end),
      strand = ifelse(test = df_introns_sub$strand == 1, yes = "+", no = "-"),
      mcols = data.frame(intron_id = df_introns_sub$phenotype_id, 
                         group_id = df_introns_sub$group_id, 
                         gene_id = df_introns_sub$gene_id)
    )
  }
  exon_list <- rlist::list.clean(exon_list, function(x) length(x) == 0L, recursive = TRUE)
  return(exon_list)
}

make_connected_components_from_cs <- function(susie_all_df, z_threshold = 3, cs_size_threshold = 10) {
  # Filter the credible sets by Z-score and size
  susie_filt <- susie_all_df %>%
    dplyr::group_by(cs_uid) %>%
    dplyr::mutate(max_abs_z = max(abs(z))) %>%
    dplyr::filter(max_abs_z > z_threshold, cs_size < cs_size_threshold) %>%
    dplyr::ungroup()
  
  # make the ranges object in order to find overlaps
  cs_ranges = GenomicRanges::GRanges(
    seqnames = susie_filt$chromosome,
    ranges = IRanges::IRanges(start = susie_filt$position, end = susie_filt$position),
    strand = "*",
    mcols = data.frame(cs_uid = susie_filt$cs_uid, variant_id = susie_filt$variant, gene_id = susie_filt$molecular_trait_id)
  )
  
  # find overlaps and remove the duplicated
  olaps <-  GenomicRanges::findOverlaps(cs_ranges, cs_ranges, ignore.strand = TRUE) %>%
    GenomicRanges::as.data.frame() %>%
    dplyr::filter(queryHits <= subjectHits)
  
  # change variant sharing into credible set sharing 
  # not to have multiple connected components of variants but credible sets
  olaps <- olaps %>% dplyr::mutate(cs_uid_1 = cs_ranges$mcols.cs_uid[queryHits], cs_uid_2 = cs_ranges$mcols.cs_uid[subjectHits])
  edge_list <- olaps %>% dplyr::select(cs_uid_1, cs_uid_2) %>% BiocGenerics::unique() %>% base::as.matrix()
  
  # make the graph of connected components
  g <- igraph::graph_from_edgelist(edge_list, directed = F)
  g_cc <- igraph::components(g)
  
  # turn connected components graph into data frame
  cc_df <- data.frame(cc_membership_no = g_cc$membership, 
                      cs_uid = g_cc$membership %>% names()) %>% 
    dplyr::mutate(molecular_trait_id = base::sub(pattern = "_[^_]+$", replacement = "",x = base::gsub(".*\\%","",cs_uid)))
  
  return(cc_df)
}




#Debugging
if (FALSE) {
  opt = list()
  opt$n = "Alasoo_2018"
  opt$f = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/susie/Alasoo_2018_leafcutter_macrophage_naive.purity_filtered.txt.gz"
  opt$s = "/Users/kerimov/Work/GitHub/SampleArcheology/studies/cleaned/Alasoo_2018.tsv"
  opt$p = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/leafcutter_metadata.txt.gz"
  opt$q = "macrophage_naive"
  opt$o = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/plots"
  opt$v = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/vcf/Alasoo_2018.filtered.vcf.gz"
  opt$b = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/bigwig/"
  opt$m = "/Users/kerimov/Work/temp_files/leafcutter_data/MANE_transcript_gene_map.txt"
  opt$g = "/Users/kerimov/Work/temp_files/leafcutter_data/Homo_sapiens.GRCh38.105.gtf"
  opt$u = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/qtl_group_split_norm/Alasoo_2018.macrophage_naive.tsv"
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
message("######### norm_usage_matrix  : ", norm_usage_matrix_path)


# coverage_plot_output_dir = paste0(output_dir, "/coverage_plots")
# if (!dir.exists(coverage_plot_output_dir)){
#   dir.create(coverage_plot_output_dir, recursive = TRUE)
# }
# 
# box_plot_output_dir = paste0(output_dir, "/box_plots")
# if (!dir.exists(box_plot_output_dir)){
#   dir.create(box_plot_output_dir, recursive = TRUE)
# }

message(" ## Reading GTF file")
gtf_ref <- rtracklayer::import(gtf_file_path, 
                               colnames = c("type", "gene_id", "gene_name", "gene_biotype", 
                                            "transcript_id", "transcript_name","transcript_biotype", "exon_number", "exon_id"),
                               feature.type = c("exon"))
message(" ## Reading sample_metadata file")
sample_metadata <- readr::read_tsv(sample_meta_path)

message(" ## Reading mane_transcript_gene_map file")
mane_transcript_gene_map <- readr::read_tsv(mane_transcript_gene_map_file)

message(" ## Reading susie file")
susie_naive <- readr::read_tsv(file = susie_file_path, col_types = "cccicccccdddddddd")

message(" ## Reading normalised usage matrix")
norm_exp_df <- readr::read_tsv(norm_usage_matrix_path)

message(" ## Reading leafcutter metadata file")
leafcutter_metadata <- readr::read_tsv(phenotype_meta_path) %>% 
  dplyr::filter(!is.na(gene_id)) %>% 
  dplyr::filter(gene_count == 1)


if (is.null(study_name)) { 
  assertthat::has_name(sample_metadata, "study" )
  study_name <- sample_metadata$study[1] 
}

message(" ## Building Connected-Components")
susie_naive <- susie_naive %>% dplyr::mutate(cs_uid = paste0(study_name, "%", cs_id))
susie_naive_cc <- make_connected_components_from_cs(susie_all_df = susie_naive)

susie_naive_filt_cc <- susie_naive %>% 
  dplyr::filter(molecular_trait_id %in% susie_naive_cc$molecular_trait_id) %>% 
  dplyr::left_join(susie_naive_cc %>% dplyr::select(cc_membership_no, molecular_trait_id))

susie_naive_high_pip_var <- susie_naive_filt_cc %>% 
  dplyr::group_by(cc_membership_no) %>% 
  dplyr::arrange(-pip) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(-pip)

susie_high_pip_with_gene <- susie_naive_high_pip_var %>% 
  dplyr::left_join(leafcutter_metadata %>% dplyr::select(-chromosome), by = c("molecular_trait_id" = "phenotype_id")) %>% 
  dplyr::filter(!is.na(gene_id))

message(" ## Starting to plot")

for (index in 1:nrow(susie_high_pip_with_gene)) {
  message("index: ", index)
  ss_oi = susie_high_pip_with_gene[index,]
  
  # get all the intons in leafcutter cluster
  cluster_introns <- leafcutter_metadata %>% dplyr::filter(group_id %in% ss_oi$group_id)
  ps_exons <- make_pseudo_exons(df_introns = cluster_introns)
  
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
    left_join(sample_metadata %>% dplyr::select(sample_id, genotype_id, qtl_group))
  
  bigwig_files <- list.files(bigwig_files_path, full.names = T)
  
  track_data_alasoo <- data.frame(file_name = (gsub(pattern = paste0(bigwig_files_path, "/"), replacement = "", x = bigwig_files)),
                                  bigWig = bigwig_files)
  track_data_alasoo <- track_data_alasoo %>% dplyr::mutate(sample_id = gsub(pattern = ".bigwig", replacement = "", x = file_name))
  
  track_data_alasoo <-track_data_alasoo %>% 
    dplyr::left_join(var_genotype %>% dplyr::select(sample_id, qtl_group, DS)) %>% 
    dplyr::mutate(scaling_factor = 1, track_id = qtl_group) %>% 
    dplyr::mutate(colour_group = as.character(DS)) %>% 
    dplyr::select(sample_id, scaling_factor, bigWig, track_id, colour_group, qtl_group)
  
  ps_exons_cdss <- make_pseudo_exons(leafcutter_metadata %>% filter(phenotype_id %in% ss_oi$molecular_trait_id))
  
  MANE_transcript_oi <- mane_transcript_gene_map %>% dplyr::filter(gene_id %in% ss_oi$gene_id) %>% dplyr::pull(transcript_id)
  mane_transcript_exons <-  make_transcript_exon_granges(gff = gtf_ref, transcript_ids = MANE_transcript_oi)
  
  path_plt = paste0(output_dir, paste0("/",gsub(pattern = ":", replacement = "_", x = ss_oi$molecular_trait_id), "&", ss_oi$variant, "&", ss_oi$gene_id))
  if (!dir.exists(path_plt)){
    dir.create(path_plt, recursive = TRUE)
  }
  
  plot_rel_height = ifelse(length(ps_exons)+1 <= 5, 3, length(ps_exons)+1) 
  plt_coverage <- wiggleplotr::plotCoverage(exons = append(ps_exons, mane_transcript_exons), cdss = ps_exons_cdss,
                                            track_data = track_data_alasoo %>% dplyr::filter(qtl_group==qtl_group_in), 
                                            heights = c(2, plot_rel_height),
                                            fill_palette = wiggleplotr::getGenotypePalette(), 
                                            coverage_type = "line")
  
  filename_plt = paste0("cov_plot_", gsub(pattern = ":", replacement = "_", x = ss_oi$molecular_trait_id), "&", ss_oi$variant, "&", ss_oi$gene_id,".pdf")
  ggplot2::ggsave(path = path_plt, filename = filename_plt, plot = plt_coverage, device = "pdf", width = 10, height = 8)

  
  norm_exp_df_oi <- norm_exp_df %>% dplyr::filter(phenotype_id %in% cluster_introns$phenotype_id)
  norm_exp_df_oi <- tibble::column_to_rownames(.data = norm_exp_df_oi,var = "phenotype_id")
  norm_exp_df_oi <- norm_exp_df_oi %>% base::t() %>% 
    GenomicRanges::as.data.frame() %>% 
    tibble::rownames_to_column(var = "sample_id")
  norm_exp_df_oi <- norm_exp_df_oi %>% 
    tidyr::pivot_longer(cols = -sample_id, names_to="intron_id", values_to = "norm_exp")
  
  track_data_alasoo_box <- track_data_alasoo %>% 
    dplyr::mutate(genotype_text = as.factor(colour_group)) %>% 
    dplyr::filter(qtl_group %in% qtl_group_in) %>%
    dplyr::mutate(condition_name = qtl_group) %>% 
    dplyr::mutate(gene_name = ss_oi$gene_name) %>% 
    dplyr::mutate(snp_id = ss_oi$variant) 
  
  track_data_alasoo_box <- norm_exp_df_oi %>%  
    dplyr::left_join(track_data_alasoo_box) %>% 
    dplyr::mutate(is_significant = intron_id == ss_oi$molecular_trait_id)
  
  message(" ## Plotting box plots")
  boxplot_facet <- ggplot2::ggplot(track_data_alasoo_box, 
                                   ggplot2::aes(x = genotype_text, 
                                                y = norm_exp, 
                                                color = is_significant, 
                                                group = genotype_text)) + 
    ggplot2::facet_wrap(~ intron_id, dir = "v") + 
    ggplot2::geom_boxplot(outlier.shape = NA) + 
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
    ggplot2::ylab(paste0("Usage of introns in cluster ", ss_oi$group_id, "\n(", ss_oi$gene_name, " gene region)")) +
    ggplot2::xlab(track_data_alasoo_box$snp_id[1]) + 
    ggplot2::theme_light() + 
    # ggplot2::scale_color_manual(values = seqUtils::conditionPalette(), guide=FALSE) +
    ggplot2::theme(strip.text.x = ggplot2::element_text(colour = "grey10"), strip.background = ggplot2::element_rect(fill = "grey85"))
  
  filename_plt_box_facet = paste0("box_facet_plot_", gsub(pattern = ":", replacement = "_", x = ss_oi$molecular_trait_id), "&", ss_oi$variant, "&", ss_oi$gene_id,".pdf")
  ggplot2::ggsave(path = path_plt, filename = filename_plt_box_facet, plot = boxplot_facet, device = "pdf", width = 12, height = 10)
  
  for (intron_id_sel in track_data_alasoo_box$intron_id %>% BiocGenerics::unique()) {
    box_plot <- ggplot2::ggplot(track_data_alasoo_box %>% dplyr::filter(intron_id == intron_id_sel), 
                                ggplot2::aes(x = genotype_text, 
                                             y = norm_exp, 
                                             color = condition_name, 
                                             group = genotype_text)) + 
      ggplot2::geom_boxplot(outlier.shape = NA) + 
      ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
      ggplot2::ylab(paste0(track_data_alasoo_box$intron_id[1], " usage")) +
      ggplot2::xlab(track_data_alasoo_box$snp_id[1]) + 
      ggplot2::theme_light() + 
      # ggplot2::scale_color_manual(values = seqUtils::conditionPalette(), guide=FALSE) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(colour = "grey10"), strip.background = ggplot2::element_rect(fill = "grey85"))
    
    filename_plt_box = paste0("box_plot_", gsub(pattern = ":", replacement = "_", x = intron_id_sel), "&", ss_oi$variant, "&", ss_oi$gene_id,".pdf")
    ggplot2::ggsave(path = path_plt, filename = filename_plt_box, plot = box_plot, device = "pdf", width = 7, height = 5)
    
  }
}
