#!/usr/bin/env Rscript
##############################################################
#  script: differential_isoform_usage.R
#  author: yld
#  last edited: 2023.12.18
#  description: differential isoform usage
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
    c("--tidy_count_table"),
    type = "character",
    default = NULL,
    help = "tidy count table based on flair quantify results",
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
    c("--adj_p_value_cutoff"),
    type = "double",
    default = 0.05,
    help = "adjusted p-value cutoff",
    metavar = "double"
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
  ),
  make_option(
    c("--interested_genes"),
    type = "character",
    default = NULL,
    help = "Genes that want to be highlighted, gene id",
    metavar = "character"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# ============================================================================ #
# test

if(FALSE) {
  opt <- list(
    tidy_count_table = "/public/yangludong/projects/test_nanotrans_updated_231122/proj_arabdopsis_tworeps/02.Isoform_Clustering_and_Quantification/Batch_Dataset2/all_samples_combined/Batch_Dataset2.all_samples_combined.counts_matrix.tidy.txt",
    contrast = "vir1,VIRc",
    read_counts_cutoff = 5,
    batch_id = "Batch_Dataset2",
    interested_genes = "FLC"
  )
}
  
# ============================================================================ #
# prep

# set output dir
print(opt$batch_id)
base_dir <- getwd()
output_dir <- file.path(
  base_dir, 
  opt$batch_id,
  "all_samples_combined", 
  paste0(opt$batch_id, "_differential_expression_output")
)

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

print("The plots & DTU table will be generated at the path: ")
print(output_dir)

# read quant table & set groups
tidy_count_table <- tryCatch({
  fread(opt$tidy_count_table)
}, error = function(err) {
  print("The quant table does not exist.")
})

required_cols <- c("isoform_id", "gene_id", "gene_name") %in% colnames(tidy_count_table)
if(!all(required_cols)) {
  print("Required columns do not exist. Please check the whether input quant table is tidy table.")
}

# sep conditions
treatment_conditions <- unlist(strsplit(opt$contrast, split = ","))

count_table <- tidy_count_table %>% 
    dplyr::select(-gene_name) %>% 
    dplyr::relocate(gene_id, .before = isoform_id) 
setnames(count_table, "isoform_id", "feature_id")
# set threads
numThread <- opt$threads

# ids
dt_ids <- tidy_count_table[, c("isoform_id", "gene_id", "gene_name")] %>% unique()

# ============================================================================ #
# calculate differential isoforms usage by DRIMseq

# build colData | rownames corresponding to colnames of count matrix
col_data <- data.frame(condition = colnames(count_table)[-c(1:2)]) %>% 
  tidyr::separate(condition, c("condition", "batch"), sep = "\\.") %>% 
  mutate(condition = factor(condition, levels = rev(treatment_conditions)))  # the first level as ref group
if (sum(is.na(col_data$condition)) > 0) print("Please verify that names of treatment conditions are equal to master table's 'comparison_group'.")
rownames(col_data) <- colnames(count_table)[-c(1:2)] 
samples <- col_data %>% tibble::rownames_to_column("sample_id")

# run 
## build dm data
data <- dmDSdata(counts = as.data.frame(count_table), samples = samples)
## filtering 
filtered <- dmFilter(data, min_samps_gene_expr = 4, min_samps_feature_expr = 2, min_gene_expr = 15, min_feature_expr = 5)
## design matrix
conditions <- samples$condition
design_full <- model.matrix(~ conditions, data = samples(filtered))
## run dm
set.seed(123)
d <- dmPrecision(filtered, design = design_full, BPPARAM=BiocParallel::MulticoreParam(numThread))
d <- dmFit(d, design = design_full, verbose = 1, BPPARAM=BiocParallel::MulticoreParam(numThread))
d <- dmTest(d, coef = colnames(design_full)[2], verbose = 1, BPPARAM=BiocParallel::MulticoreParam(numThread))
## tidyr results
res <- merge(proportions(d),results(d,level="feature"), by=c("feature_id","gene_id")) %>%
  left_join(unique(dt_ids[, c("isoform_id", "gene_name")]), by = c("feature_id" = "isoform_id")) %>% 
  mutate(comparison_group_pair = paste0(treatment_conditions[1], "/", treatment_conditions[2])) %>% 
  dplyr::relocate(comparison_group_pair,lr, df, pvalue, adj_pvalue, .after = gene_name) %>% 
  arrange(adj_pvalue) %>% 
  as.data.table()
setnames(res, old = "feature_id", new = "isoform_id")

# ============================================================================ #
# save

# DTU table
fwrite(res, file = paste0(output_dir, "/", opt$batch_id, ".", treatment_conditions[1], "_vs_", treatment_conditions[2], "_comparison.differential_isoform_usage.tidy.txt"))

# drimseq data
save(d, file = paste0(output_dir, "/", opt$batch_id, ".", treatment_conditions[1], "_vs_", treatment_conditions[2], "_comparison.differential_isoform_usage.drimseq.rdt"))
load("/public/yangludong/projects/test_nanotrans_updated_231122/proj_arabdopsis_tworeps/03.Isoform_Expression_and_Splicing_Comparison/Batch_Dataset2/all_samples_combined/Batch_Dataset2_differential_expression_output/Batch_Dataset2.vir1_vs_VIRc_comparison.differential_isoform_usage.drimseq.rdt")
# ============================================================================ #
# option | plot DTU for interested genes

if(!is.null(opt$interested_genes)) {
  genes_to_label <- unlist(strsplit(opt$interested_genes, split = ","))
  if(genes_to_label %in% names(d)) {
    pdf(paste0(output_dir, "/", opt$batch_id, ".", treatment_conditions[1], "_vs_", treatment_conditions[2], "_comparison.differential_isoform_usage_interested_genes.pdf"))
      plotProportions(d, gene_id = genes_to_label, group_variable = "condition", plot_type = "ribbonplot")
    dev.off()
  }
}

  
