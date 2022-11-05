message(" ## Loading libraries: optparse")
suppressPackageStartupMessages(library("optparse"))

#Parse command-line options
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  make_option(c("-f", "--finemap_susie"), type="character", default=NULL,
              help="Purity filtered susie output. Tab separated file", metavar = "type"),
  make_option(c("-p", "--phenotype_meta"), type="character", default=NULL,
              help="Phenotype metadata file. Tab separated file", metavar = "type"),
  make_option(c("-q", "--qtl_group"), type="character", default=NULL,
              help="The selected qtl_group in the study", metavar = "type"),
  make_option(c("--quant_method"), type="character", default=NULL,
              help="The quantification method used", metavar = "type"),
  make_option(c("-c", "--chunk_size"), type="integer", default=20,
              help="Perform analysis in chunks. Eg value 5 10 would indicate that phenotypes are split into 10 chunks and the 5th one of those will be processed. [default \"%default\"]", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default="./batched_susie/",
              help="Path to the output directory. [default \"%default\"]", metavar = "type"),
  make_option(c("-d", "--debug_mode"), type="logical", 
              action = "store_true", default=FALSE,
              help="If to run the script in debug mode", metavar = "type")
)

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, cqn, data.table")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))


#Debugging
if (FALSE) {
  opt = list()
  opt$f = "/Users/kerimov/Work/temp_files/debug/new/Alasoo_2018_leafcutter_macrophage_IFNg.purity_filtered.txt.gz"
  opt$p = "/Users/kerimov/Work/temp_files/debug/new/phenotype_metadata.txt.gz"
  opt$q = "macrophage_IFNg"
  opt$o = "batched_susie/"
  opt$chunk_size = 20
  opt$debug_mode = TRUE
}

susie_file_path = opt$f
phenotype_meta_path = opt$p
output_dir = opt$o
qtl_group_in = opt$q
debug_mode = opt$debug_mode
chunk_size = opt$chunk_size
quant_method = opt$quant_method

message("######### Options: ######### ")
message("######### Working Directory  : ", getwd())
message("######### qtl_group          : ", qtl_group_in)
message("######### quant_method       : ", quant_method)
message("######### susie_file_path    : ", susie_file_path)
message("######### phenotype_meta_path: ", phenotype_meta_path)
message("######### output_dir         : ", output_dir)
message("######### debug_mode         : ", debug_mode)
message("######### chunk_size         : ", chunk_size)


message(" ## Reading susie_purity_filtered file")
susie_purity_filtered <- readr::read_tsv(file = susie_file_path, col_types = "cccdcccccdddddddd")

phenotype_metadata_col_types = ifelse(test = quant_method == "leafcutter", yes = "cccccddiccidddddddd", no = "cccccddicciddi")
if (quant_method == "tx" || quant_method == "txrev") {
  phenotype_metadata_col_types = "cccccddiccid"
}

message(" ## Reading leafcutter metadata file")
phenotype_metadata <- readr::read_tsv(phenotype_meta_path, col_types = phenotype_metadata_col_types) 

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
  dplyr::ungroup() %>% 
  dplyr::filter(!is.na(position))

if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

if (debug_mode) {
  message(" ## Slicing 10 highest pip credible set variants for debug_mode")
  highest_pip_vars_per_cs <- highest_pip_vars_per_cs %>% 
    dplyr::arrange(-pip) %>% 
    dplyr::slice_head(n = 10)
  filename = file.path(output_dir, "single_batch_debug_mode.tsv")
  readr::write_tsv(highest_pip_vars_per_cs, file = filename)
} else {
  n_chunks = (nrow(highest_pip_vars_per_cs) / chunk_size) %>% ceiling()
  for (i in 1:n_chunks) {
    message(" ## Writing batch number: ", i, "_", n_chunks)
    chunk = highest_pip_vars_per_cs[((i-1)*chunk_size+1): (min((i)*chunk_size, nrow(highest_pip_vars_per_cs))),]
    filename = file.path(output_dir, paste0("batch_",i,"_",n_chunks,".tsv"))
    readr::write_tsv(chunk, file = filename)
  }
}



