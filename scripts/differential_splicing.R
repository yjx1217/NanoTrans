#!/usr/bin/env Rscript
##############################################################
#  script: differential_splicing.R
#  author: yld
#  last edited: 2023.12.22
#  description: differential splicing
#  example: 
##############################################################

# ============================================================================ #
# load

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DRIMSeq))

# ============================================================================ #
# parse args

option_list <- list(
  make_option(
    c("--inputs_dir"),
    type = "character",
    default = NULL,
    help = "location of alternative splicing quant tables",
    metavar = "character"
  ),
  make_option(
    c("--threads"),
    type = "integer",
    default = 1,
    help = "number of threads",
    metavar = "integer"
  ),
  make_option(
    c("--contrast"),
    type = "character",
    default = "treated,control",
    help = "contrast specification for differential expression detection",
    metavar = "character"
  ),
  make_option(
    c("--batch_id"),
    type = "character",
    default = NULL,
    help = "dataset name",
    metavar = "character"
  )   
)

opt <- parse_args(OptionParser(option_list = option_list))

# ============================================================================ #
# logic

# inputs
## as_count_table, contrast, batch_id

# acts
##  extract sample ids (eg. vir1.rep1_vir1_rep1) from cols of as_count_table
##  

#TODO add a pseudocount of 1 to each event in each sample
# ============================================================================ #
# test

if(FALSE) {
    opt <- list(
        inputs_dir = "/public/yangludong/projects/test_nanotrans_updated_231122/proj_arabdopsis_linkraw/03.Isoform_Expression_and_Splicing_Comparison/Batch_Dataset2/all_samples_combined/Batch_Dataset2.VIRc-vir1_comparison.differential_isoform_splicing/",
        contrast = "vir1,VIRc",
        batch_id = "Batch_Dataset2",
        threads = 4
    )
}
# ============================================================================ #
# functions

run_drimseq <- function(
  count_table, # path of {es/ir/alt5/alt3} quant csv file
  as_type  # gene or isoform to specific id type of input table
) {

  # ids
  dt_ids <- count_table[, c("feature_id", "isoform_ids")] %>% unique()

  # build countData | cols: feature_id (as_coords), gene_id, counts cols..., isoform_ids.
  count_table <- count_table %>% dplyr::rename(gene_id = coordinate) 
  
  # build colData | rownames corresponding to colnames of count matrix
  col_data <- data.frame(condition = colnames(count_table)[-c(1:2,ncol(count_table))]) %>% 
    tidyr::separate(condition, c("tmp", "condition", "batch"), sep = "_") %>% 
    dplyr::select(-tmp) %>% 
    mutate(condition = factor(condition, levels = rev(treatment_conditions)))  # the first level as ref group
  
  if (sum(is.na(col_data$condition)) > 0) print("Please verify that names of treatment conditions are equal to master table's 'comparison_group'.")
  rownames(col_data) <- colnames(count_table)[-c(1:2, ncol(count_table))] 
  samples <- col_data %>% tibble::rownames_to_column("sample_id")
  
  # build dmDSdata
  data <- dmDSdata(counts = as.data.frame(count_table), samples = samples)
  
  # filter
  filtered <- dmFilter(data, min_samps_gene_expr = 4, min_samps_feature_expr = 2, min_gene_expr = 15, min_feature_expr = 5)
  
  # design matrix
  conditions <- samples$condition
  design_full <- model.matrix(~ conditions, data = samples(filtered))
  
  # run dm
  set.seed(123)
  
  d <- dmPrecision(filtered, design = design_full, BPPARAM=BiocParallel::MulticoreParam(numThread))
  d <- dmFit(d, design = design_full, verbose = 1, BPPARAM=BiocParallel::MulticoreParam(numThread))
  d <- dmTest(d, coef = colnames(design_full)[2], verbose = 1, BPPARAM=BiocParallel::MulticoreParam(numThread))
  
  # tidy res
  res <- merge(proportions(d),results(d,level="feature"), by=c("feature_id","gene_id")) %>% 
    left_join(dt_ids, by = "feature_id") %>% 
    mutate(comparison_group_pair = paste0(treatment_conditions[1], "/", treatment_conditions[2])) %>% 
    mutate(alt_splicing_class = as_type) %>% 
    dplyr::relocate(isoform_ids, alt_splicing_class, comparison_group_pair,lr, df, pvalue, adj_pvalue, .after = gene_id) %>% 
    arrange(adj_pvalue) %>% 
    as.data.table()
  setnames(res, old = "gene_id", new = "coordinate")

}

# ============================================================================ #
# prep

# set threads
numThread <- opt$threads

# set output dir
print(opt$batch_id)
base_dir <- getwd()
output_dir <- file.path(
  base_dir, 
  opt$batch_id,
  "all_samples_combined", 
  paste0(opt$batch_id, "_differential_splicing_output")
)

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

print("The plots & differential splicing table will be generated at the path: ")
print(output_dir)

# sep conditions
treatment_conditions <- unlist(strsplit(opt$contrast, split = ","))

# read quant table & set groups
as_quant_files <- list.files(path = opt$inputs_dir, pattern = "diffsplice.(ir|es|alt5|alt3).events.quant.tsv", full.names = TRUE)

## test if the files exsits
stopifnot(length(as_quant_files) > 0)
#tryCatch({
#  if(length(as_quant_files) == 0){
#    stop("The splicing quant tables do not exist.")
#  }
#  }, error = function(e) {
#    q(status = 1, save = "no")
#  })

as_quant_list <- lapply(as_quant_files, fread)
as_list_names <- gsub("(diffsplice.|.events.quant.tsv)", "", basename(as_quant_files))
names(as_quant_list) <- as_list_names
#required_cols <- c("isoform_id", "gene_id", "gene_name") %in% colnames(tidy_count_table)
#if(!all(required_cols)) {
#  print("Required columns do not exist. Please check the whether input quant table is tidy table.")
#}

# ============================================================================ #
# run DRIMSeq

dmres <- lapply(
  1:length(as_quant_list), 
  function(x) run_drimseq(as_quant_list[[x]], as_list_names[x])
)
dt_dmres <- do.call(rbind, dmres)

# ============================================================================ #
# save
fwrite(dt_dmres, file = paste0(output_dir, "/", opt$batch_id, ".", treatment_conditions[1], "_vs_", treatment_conditions[2], "_comparison.differential_isoform_splicing.all.events.quant.tidy.txt"), sep = "\t")

