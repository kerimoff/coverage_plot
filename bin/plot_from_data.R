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
  make_option(c("-i", "--individual_boxplots"), type="logical", 
              action = "store_true", default=FALSE,
              help="Flag to generate individual boxplots", metavar = "type")
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
                                        xlabel = "Distance from region start (bp)",
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
  opt$r = "/Users/kerimov/Work/temp_files/debug/generate_from_data_all_quants/tx/ENST00000082468&chr6_26417623_C_T&ENSG00000026950/plot_data_ENST00000082468&chr6_26417623_C_T&ENSG00000026950.Rds"
  opt$t = "/Users/kerimov/Work/temp_files/debug/generate_from_data_all_quants/tx/ENST00000082468&chr6_26417623_C_T&ENSG00000026950/plot_data_ENST00000082468&chr6_26417623_C_T&ENSG00000026950.tar.gz"
  opt$o = "/Users/kerimov/Work/temp_files/debug/generate_from_data_all_quants/plots_generated"
  opt$individual_boxplots = FALSE
}

rds_file_path = opt$r
tar_file_path = opt$t
output_dir = opt$o
individual_boxplots <- opt$individual_boxplots

message("######### Options: ######### ")
message("######### Working Directory  : ", getwd())
message("######### rds_file_path      : ", rds_file_path)
message("######### tar_file_path      : ", tar_file_path)
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
  
  coverage_data_list <- Rds_data_list$coverage_data_list
  nom_exon_cc_sumstats_df <- Rds_data_list$nom_exon_cc_sumstats_df
  box_plot_wrap <- Rds_data_list$box_plot_wrap_df
  ss_oi <- Rds_data_list$ss_oi
  limits_loc <- coverage_data_list$limits
  
  tx_structure_df = prepareTranscriptStructureForPlotting(exon_ranges = coverage_data_list$tx_annotations$exon_ranges, 
                                                          cds_ranges = coverage_data_list$tx_annotations$cds_ranges, 
                                                          transcript_annotations = coverage_data_list$plotting_annotations)
  
  # intron_ss_oi_vert_lines = tx_structure_df %>% 
  #   dplyr::filter(transcript_id == ss_oi$molecular_trait_id, feature_type == "exon") 
  # intron_ss_oi_vert_lines <- c(intron_ss_oi_vert_lines[1,] %>% dplyr::pull(end), 
  #                              intron_ss_oi_vert_lines[nrow(intron_ss_oi_vert_lines),] %>% dplyr::pull(start))
  
  
  wiggle_plots <- wiggleplotr::plotCoverageData(coverage_data_list, alpha = 1, 
                                                fill_palette = wiggleplotr::getGenotypePalette(), 
                                                coverage_type = "line", return_subplots_list = TRUE, 
                                                show_legend = TRUE)
  
  # exon_plot <- wiggle_plots$tx_structure + ggplot2::geom_vline(xintercept = intron_ss_oi_vert_lines, alpha = 0.5, color = "lightgrey")
  exon_plot <- wiggle_plots$tx_structure 
  
  coverage_data_list$coverage_df <- coverage_data_list$coverage_df %>% dplyr::filter(!is.na(coverage))
  # coverage_plot <- wiggle_plots$coverage_plot + ggplot2::geom_vline(xintercept = intron_ss_oi_vert_lines, alpha = 0.5, color = "lightgrey")
  coverage_plot <- wiggle_plots$coverage_plot
  
  # plot_rel_height = ifelse(length(tx_structure_df$transcript_id %>% unique()) <= 5, 3, length(transcript_struct_df$transcript_id %>% unique()))
  
  plot_rel_height = ifelse(length(tx_structure_df$transcript_id %>% unique())-1 <= 5, 3, length(tx_structure_df$transcript_id %>% unique())) 
  plot_rel_height = ifelse(plot_rel_height > 20, 20, plot_rel_height)  
} else { # using TAR file to plot 
  tar_files <- untar(tarfile = tar_file_path, list = T)
  path_to_export = file.path(output_dir, gsub("\\..*","",basename(tar_file_path)))
  untar(tarfile = tar_file_path, exdir = path_to_export)
  tar_files <- file.path(path_to_export, tar_files)
  
  tx_structure_df = read_tsv(tar_files[which(grepl(pattern = "tx_str", x = tar_files))], col_types = "ciiiiccdcd", na = character())
  coverage_df = read_tsv(tar_files[which(grepl(pattern = "coverage_df", x = tar_files))], col_types = "ccidf")
  limits_loc = c(1, tx_structure_df$limit_max[1])
  box_plot_wrap = read_tsv(tar_files[which(grepl(pattern = "box_plot", x = tar_files))])
  ss_oi = read_tsv(tar_files[which(grepl(pattern = "ss_oi", x = tar_files))])
  nom_exon_cc_sumstats_df = read_tsv(tar_files[which(grepl(pattern = "nom_exon_cc", x = tar_files))])
  
  unlink(path_to_export, recursive = TRUE)
  
  # intron_ss_oi_vert_lines = tx_structure_df %>% 
  #   dplyr::filter(transcript_id == ss_oi$molecular_trait_id, feature_type == "exon") 
  # intron_ss_oi_vert_lines <- c(intron_ss_oi_vert_lines[1,] %>% dplyr::pull(end), 
  #                              intron_ss_oi_vert_lines[nrow(intron_ss_oi_vert_lines),] %>% dplyr::pull(start))
  
  coverage_plot = makeCoveragePlot_loc(coverage_df = coverage_df, 
                                       limits = limits_loc, 
                                       alpha = 1, 
                                       fill_palette = wiggleplotr::getGenotypePalette(), 
                                       coverage_type = "line", 
                                       show_genotype_legend = TRUE,
                                       vert_lines = NULL)
  
  exon_plot <- plotTranscriptStructure_loc(exons_df = tx_structure_df, 
                                           limits = limits_loc, 
                                           vert_lines = NULL)

  plot_rel_height = ifelse(length(tx_structure_df$transcript_id %>% unique())-1 <= 5, 3, length(tx_structure_df$transcript_id %>% unique())) 
  plot_rel_height = ifelse(plot_rel_height > 20, 20, plot_rel_height)  
}

if (exists("nom_exon_cc_sumstats_df") && nrow(nom_exon_cc_sumstats_df) > 0) {
  beta_plot <- generate_beta_plot(transcript_struct_df_loc = tx_structure_df, 
                                  nom_exon_cc_sumstats_filt_loc = nom_exon_cc_sumstats_df, 
                                  limits = limits_loc, 
                                  ss_oi = ss_oi,
                                  vert_lines = NULL)
  
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
                                    "\nSE: ", se)) 

if (box_plot_wrap %>% hasName("tx_id")) {
    box_plot_wrap <- box_plot_wrap %>% dplyr::mutate(intron_id = tx_id)
}
box_plot_wrap <- box_plot_wrap %>% dplyr::mutate(intron_id_with_stats = paste0(intron_id, "\n", stats_text))

boxplot_facet <- ggplot2::ggplot(box_plot_wrap, 
                                 ggplot2::aes(x = genotype_text %>% as.factor(), 
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

if (individual_boxplots) {
  message(" ## Plotting individual boxplots")
  
  for (intron_id_sel in box_plot_wrap$intron_id %>% BiocGenerics::unique()) {
    nom_cc_sumstats_filt <- box_plot_wrap %>% dplyr::filter(intron_id == intron_id_sel) 
    
    labels <- c(paste0(unique(nom_cc_sumstats_filt$stats_text), "				MAF: ", nom_cc_sumstats_filt$maf[1]))
    signifincance_color = ifelse(unique(nom_cc_sumstats_filt$is_significant), "#00BFC4", "#F8766D")
    box_plot <- ggplot2::ggplot(nom_cc_sumstats_filt, 
                                ggplot2::aes(x = genotype_text, 
                                             y = norm_exp, 
                                             group = genotype_text)) + 
      ggplot2::geom_boxplot(outlier.shape = NA, color = signifincance_color) + 
      ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5, color = signifincance_color) + 
      ggplot2::ylab(paste0(nom_cc_sumstats_filt$intron_id[1], "\nCounts Per Million")) +
      ggplot2::xlab(nom_cc_sumstats_filt$snp_id[1]) + 
      ggplot2::theme_light() + 
      ggplot2::labs(subtitle = labels) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(colour = "grey10"), strip.background = ggplot2::element_rect(fill = "grey85"))
    
    filename_plt_box = paste0("box_plot_", gsub(pattern = ":", replacement = "_", x = intron_id_sel), "&", ss_oi$variant, "&", ss_oi$gene_id,".pdf")
    ggplot2::ggsave(path = path_plt, filename = filename_plt_box, plot = box_plot, device = "pdf", width = 5, height = 5)
  }

}
