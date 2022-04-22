message(" ## Loading libraries: optparse")
suppressPackageStartupMessages(library("optparse"))

#Parse command-line options
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  make_option(c("-r", "--rds_file"), type="character", default=NULL,
              help="Rds object path containing plot data", metavar = "type"),
  make_option(c("-t", "--tar_file"), type="character", default=NULL,
              help="tar object path containing plot data", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default="./normalised_results/",
              help="Path to the output directory. [default \"%default\"]", metavar = "type"),
  make_option(c("-e", "--exon_summ_stats"), type="character", default=NULL,
              help="Full exon summary statistics file to extract exon betas", metavar = "type")
)

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, cqn, data.table")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ggplot2"))

conf.level = 0.95
ci.value <- -qnorm( ( 1 - conf.level ) / 2 )

dataTrackTheme <- function(){
  theme = theme(axis.text.x = element_blank(), 
                axis.title.x = element_blank(), 
                axis.ticks.x = element_blank(),
                plot.margin=unit(c(0.1,1,0.1,1),"line"),
                legend.position="none",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.text.y = element_text(colour = "grey10"),
                strip.background = element_rect(fill = "grey85"))
  return(theme)
}

makeCoveragePlot_loc <- function(coverage_df, limits, 
                                 alpha = 1, 
                                 fill_palette = c("#a1dab4","#41b6c4","#225ea8"), 
                                 coverage_type = "area", 
                                 show_genotype_legend = FALSE, 
                                 vert_lines = NULL,
                                 vert_line_alpha = 0.3){
  coverage_max_y = max(coverage_df$coverage) + max(coverage_df$coverage) * 0.05
  
  #Plot coverage over a region
  coverage_plot = ggplot(coverage_df, aes_(~bins, ~coverage, group = ~sample_id, alpha = ~alpha)) + 
    geom_blank() +
    theme_light()
  #Choose between plotting a line and plotting area
  if(coverage_type == "line"){
    coverage_plot = coverage_plot + 
      geom_line(aes_(colour = ~colour_group), alpha = alpha, position = "identity") 
  } else if (coverage_type == "area"){
    coverage_plot = coverage_plot + 
      geom_area(aes_(fill = ~colour_group), alpha = alpha, position = "identity")
  } else if (coverage_type == "both"){
    coverage_plot = coverage_plot + 
      geom_area(aes_(fill = ~colour_group), alpha = alpha, position = "identity") +
      geom_line(aes_(colour = ~colour_group), alpha = alpha, position = "identity") 
  } else{
    stop("Coverage type not supported.")
  }
  coverage_plot = coverage_plot +
    facet_grid(track_id~.) +
    dataTrackTheme() + 
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, coverage_max_y)) +
    coord_cartesian(xlim = limits) +
    scale_color_manual(values = fill_palette) +
    scale_fill_manual(values = fill_palette) +
    ylab("FPM")
  
  # Add genotype legend if show_genotype_legend
  if(show_genotype_legend) {
    coverage_plot = coverage_plot + 
      labs(colour = NULL) +
      theme(legend.justification = c(1, 1), 
            legend.position = "right", 
            legend.direction = "vertical", 
            legend.title.align = 0, legend.box = "vertical",
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.margin=margin(t=0, r=0, b=0, l=0))  
  }
  
  if (!is.null(vert_lines)) {
    coverage_plot <- coverage_plot + ggplot2::geom_vline(xintercept = vert_lines, alpha = vert_line_alpha, color = "lightgrey")
  }
  return(coverage_plot)
}

plotTranscriptStructure_loc <- function(exons_df, 
                                        limits = NA, 
                                        connect_exons = TRUE,  
                                        xlabel = "Distance from gene start (bp)", 
                                        transcript_label = TRUE, 
                                        vert_lines = NULL,
                                        vert_line_alpha = 0.3){
  
  #Extract the position for plotting transcript name
  transcript_annot = dplyr::group_by_(exons_df, ~transcript_id) %>% 
    dplyr::filter_(~feature_type == "exon") %>%
    dplyr::arrange_('transcript_id', 'start') %>%
    dplyr::filter(row_number() == 1)
  
  #Create a plot of transcript structure
  plot = ggplot(exons_df) + geom_blank()
  if(connect_exons){ #Print line connecting exons
    plot = plot + geom_line(aes_(x = ~start, y = ~transcript_rank, group = ~transcript_rank, color = ~feature_type))
  }
  plot = plot + 
    geom_rect(aes_(xmin = ~start, 
                   xmax = ~end, 
                   ymax = ~transcript_rank + 0.25, 
                   ymin = ~transcript_rank - 0.25, 
                   fill = ~feature_type)) + 
    theme_light() +
    theme(plot.margin=unit(c(0,1,1,1),"line"), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(colour = "grey10"),
          strip.background = element_rect(fill = "grey85")) +
    xlab(xlabel) +
    facet_grid(type~.) +
    scale_y_continuous(expand = c(0.2,0.15)) +
    scale_fill_manual(values = c("#2c7bb6","#abd9e9")) + 
    scale_colour_manual(values = c("#2c7bb6","#abd9e9"))
  if(all(!is.na(limits))){
    plot = plot + scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(xlim = limits)
  }
  if(transcript_label){
    plot = plot + geom_text(aes_(x = ~start, 
                                 y = ~transcript_rank + 0.30, 
                                 label = ~transcript_label), 
                            data = transcript_annot, hjust = 0, vjust = 0, size = 4)
    
  }
  
  if (!is.null(vert_lines)) {
    plot <- plot + ggplot2::geom_vline(xintercept = vert_lines, alpha = vert_line_alpha, color = "lightgrey")
  }
  return(plot)
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
if (TRUE) {
  opt = list()
  opt$r = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/plots_data_generate/18_333143_334742_clu_4553_+&chr18_334742_C_T&ENSG00000158270/plot_data_18_333143_334742_clu_4553_+&chr18_334742_C_T&ENSG00000158270.Rds"
  opt$t = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/plots_data_generate/18_333143_334742_clu_4553_+&chr18_334742_C_T&ENSG00000158270/plot_data_18_333143_334742_clu_4553_+&chr18_334742_C_T&ENSG00000158270.tar.gz"
  opt$e = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/sumstats/Alasoo_2018_exon_macrophage_naive.all.tsv.gz"
  opt$o = "/Users/kerimov/Work/temp_files/leafcutter_data/Alasoo_2018_17Feb/plots_from_data"
}

rds_file_path = opt$r
tar_file_path = opt$t
nominal_exon_sumstats_path = opt$e
output_dir = opt$o

message("######### Options: ######### ")
message("######### Working Directory  : ", getwd())
message("######### rds_file_path      : ", rds_file_path)
message("######### tar_file_path      : ", tar_file_path)
message("######### exon_sumstats      : ", nominal_exon_sumstats_path)
message("######### output_dir         : ", output_dir)


useRds = TRUE
if (is.null(rds_file_path) && is.null(tar_file_path)) {
  message("No RDS or TAR file is provided. Please provide one!")
  quit(save = "no")
}

if (!is.null(rds_file_path) && !is.null(tar_file_path)) {
  message("Both RDS or TAR files are provided. Will use Rds file: ", rds_file_path)
}

if (!is.null(rds_file_path) && is.null(tar_file_path)) {
  message("Plotting using Rds file: ", rds_file_path)
}

if (is.null(rds_file_path) && !is.null(tar_file_path)) {
  message("Plotting using TAR file: ", tar_file_path)
  useRds = FALSE
}

message(" ## Starting to plot")
if (useRds) {
  Rds_data_list <- readr::read_rds(file = rds_file_path)  
  
  transcript_struct_df <- Rds_data_list$tx_str_df
  limits_loc <- c(1, transcript_struct_df$limit_max[1])
  coverage_df <- Rds_data_list$coverage_plot_df
  box_plot_wrap <- Rds_data_list$box_plot_wrap_df
  ss_oi <- Rds_data_list$ss_oi
  
} else { # using TAR file to plot 
  tar_files <- untar(tarfile = tar_file_path, list = T)
  path_to_export = file.path(output_dir, gsub("\\..*","",basename(tar_file_path)))
  untar(tarfile = tar_file_path, exdir = path_to_export)
  tar_files <- file.path(path_to_export, tar_files)
  transcript_struct_df = read_tsv(tar_files[which(grepl(pattern = "tx_str", x = tar_files))], col_types = "ciiiiccdcd", na = character())
  coverage_df = read_tsv(tar_files[which(grepl(pattern = "coverage_df", x = tar_files))], col_types = "ccidf")
  limits_loc = c(1, transcript_struct_df$limit_max[1])
  box_plot_wrap = read_tsv(tar_files[which(grepl(pattern = "box_plot", x = tar_files))])
  ss_oi = read_tsv(tar_files[which(grepl(pattern = "ss_oi", x = tar_files))])
  
  unlink(path_to_export, recursive = TRUE)
}



intron_ss_oi_vert_lines = transcript_struct_df %>% 
  dplyr::filter(transcript_id == ss_oi$molecular_trait_id, feature_type == "exon") 
intron_ss_oi_vert_lines <- c(intron_ss_oi_vert_lines[1,] %>% pull(end), intron_ss_oi_vert_lines[2,] %>% pull(start))

coverage_plot = makeCoveragePlot_loc(coverage_df = coverage_df, 
                                     limits = limits_loc, 
                                     alpha = 1, 
                                     fill_palette = wiggleplotr::getGenotypePalette(), 
                                     coverage_type = "line", 
                                     show_genotype_legend = TRUE,
                                     vert_lines = intron_ss_oi_vert_lines)

exon_plot <- plotTranscriptStructure_loc(exons_df = transcript_struct_df, 
                                         limits = limits_loc, 
                                         vert_lines = intron_ss_oi_vert_lines)

plot_rel_height = ifelse(length(transcript_struct_df$transcript_id %>% unique()) <= 5, 3, length(transcript_struct_df$transcript_id %>% unique())) 

sumstat_colnames <- c("molecular_trait_id", "chromosome", "position", 
                      "ref", "alt", "variant", "ma_samples", "maf", 
                      "pvalue", "beta", "se", "type", "ac", "an", 
                      "r2", "molecular_trait_object_id", "gene_id", 
                      "median_tpm", "rsid")

if(!is.null(nominal_exon_sumstats_path)) {
  message(" ## Reading exon summary statistics")
  variant_regions_ss <- ss_oi %>% 
    dplyr::select(variant, chromosome, position) %>% 
    dplyr::mutate(region = paste0(chromosome,":", position, "-", position))
  
  nom_exon_cc_sumstats <- seqminer::tabix.read.table(nominal_exon_sumstats_path, variant_regions_ss$region) 
  colnames(nom_exon_cc_sumstats) <- sumstat_colnames
  
  # Extract the QTLs of exons according to gene and variant of interest
  nom_exon_cc_sumstats_filt <- nom_exon_cc_sumstats %>% 
    dplyr::filter(variant == ss_oi$variant, molecular_trait_object_id == ss_oi$gene_id) %>% 
    dplyr::filter(rsid == rsid[1]) %>% # if variant has more than one rsid keep only the first unique rsid 
    dplyr::mutate(n_uniq_pheno = length(unique(molecular_trait_id))) %>% 
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
  
}

if (exists("nom_exon_cc_sumstats_filt") && nrow(nom_exon_cc_sumstats_filt) > 0) {
  beta_plot <- generate_beta_plot(transcript_struct_df_loc = transcript_struct_df, 
                                  nom_exon_cc_sumstats_filt_loc = nom_exon_cc_sumstats_filt, 
                                  limits = limits_loc, 
                                  vert_lines = intron_ss_oi_vert_lines)
  
  merged_plot <- cowplot::plot_grid(coverage_plot, beta_plot, exon_plot , align = "v", axis = "lr", rel_heights = c(3, 3, plot_rel_height), ncol = 1)
} else {
  merged_plot <- cowplot::plot_grid(coverage_plot, exon_plot , align = "v", axis = "lr", rel_heights = c(3, plot_rel_height), ncol = 1)
}

signal_name <- paste0(gsub(pattern = ":", replacement = "_", x = ss_oi$molecular_trait_id), "&", ss_oi$variant, "&", ss_oi$gene_id)
path_plt = file.path(output_dir, signal_name)
if (!dir.exists(path_plt)){
  dir.create(path_plt, recursive = TRUE)
}
filename_plt = paste0("cov_plot_", signal_name,".pdf")
ggplot2::ggsave(path = path_plt, filename = filename_plt, plot = merged_plot, device = "pdf", width = 10, height = 8)

box_plot_wrap <- box_plot_wrap %>% 
  dplyr::mutate(stats_text = paste0("Pval: ", pvalue, "			BETA: ", beta, 
                                    "\nSE: ", se)) %>% 
  dplyr::mutate(intron_id_with_stats = paste0(intron_id, "\n", stats_text))

boxplot_facet <- ggplot2::ggplot(box_plot_wrap, 
                                 ggplot2::aes(x = genotype_text, 
                                              y = norm_exp, 
                                              color = is_significant, 
                                              group = genotype_text)) + 
  ggplot2::facet_wrap(~ intron_id_with_stats, dir = "v") + 
  ggplot2::geom_boxplot(outlier.shape = NA) + 
  ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
  ggplot2::ylab(paste0("Usage of introns in cluster ", ss_oi$group_id, "\n(", ss_oi$gene_name, " gene region)")) +
  ggplot2::xlab(paste0(box_plot_wrap$snp_id[1], "    (MAF:", box_plot_wrap$maf[1], ")")) + 
  ggplot2::theme_light() + 
  ggplot2::theme(strip.text.x = ggplot2::element_text(colour = "grey10"), strip.background = ggplot2::element_rect(fill = "grey85"))

filename_plt_box_facet = paste0("box_facet_plot_", signal_name,".pdf")
ggplot2::ggsave(path = path_plt, filename = filename_plt_box_facet, plot = boxplot_facet, device = "pdf", width = 10, height = 11)

