#!/usr/bin/env Rscript
##############################################################
#  script: differential_expression.R
#  author: yld
#  last edited: 2023.12.15
#  description: calculate differentially expressed genes/isoforms
#  example: 
##############################################################

# ============================================================================ #
# load

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))

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
    c("--master_table"),
    type = "character",
    default = NULL,
    help = "name of the sample table file",
    metavar = "character"
  ),
  make_option(
    c("--read_counts_cutoff"),
    type = "integer",
    default = NULL,
    help = "minimum reads coverage for the genes across all samples at least in one condition",
    metavar = "integer"
  ),
  make_option(
    c("--log2foldchange_cutoff"),
    type = "double",
    default = "1",
    help = "log2foldchange cutoff",
    metavar = "double"
  ),
  make_option(
    c("--adj_p_value_cutoff"),
    type = "double",
    default = "0.05",
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
    help = "Genes that want to be highlighted",
    metavar = "character"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# ============================================================================ #
# test

if (FALSE) {
  opt <- list(
    tidy_count_table = "/public/yangludong/projects/test_nanotrans_updated_240108/proj_gsa_PRJCA017047_downs_10reads/02.Isoform_Clustering_and_Quantification/downs10reads_PRJCA017047/all_samples_combined/downs10reads_PRJCA017047.all_samples_combined.counts_matrix.tidy.txt",
    contrast = "mu,wt",
    read_counts_cutoff = 5,
    batch_id = "Batch_debug",
    master_table = ""
  )
}

# ============================================================================ #
# functions

# TEST args for run_deseq2
if (FALSE) {
  count_table = gene_count_table
  read_counts_cutoff = 5
  groups = treatment_conditions
  ids_type = "gene"
  outdir = output_dir
}


run_deseq2 <- function(
  count_table, # tidy count table based on flair quantity results, retain 1 id col only
  read_counts_cutoff,  # loweast read coverage on the gene / isoform in all samples in at least one condition
  groups, # treatment_conditions
  ids_type,  # gene or isoform to specific id type of input table
  outdir
) {

  # build colData | rownames corresponding to colnames of count matrix
  col_data <- data.frame(sample_ids = colnames(count_table)[-1]) %>% 
    left_join(master_table[, c("#sample_id", "comparison_group", "replicate_id")], by = c("sample_ids" = "#sample_id")) %>% 
    mutate(condition = comparison_group) %>% 
    mutate(condition = factor(condition, levels = rev(groups)))  # set the 1th level as the ref group
  if (sum(is.na(col_data$condition)) > 0) print("Please verify that names of treatment conditions are equal to master table's 'comparison_group'.")
  rownames(col_data) <- colnames(count_table)[-1]
  
  # filter count data
  ## Column groups
  col_group1 <- rownames(col_data[which(col_data[["condition"]] == groups[1]), ])
  col_group2 <- rownames(col_data[which(col_data[["condition"]] == groups[2]), ])
  
  ## Only retain genes/isoforms whose reads coverage greater than 'cutoff' in all samples under at least one condition
  filtered_count_table <- count_table[
    rowSums(count_table[, ..col_group1] > read_counts_cutoff) == length(col_group1) |
    rowSums(count_table[, ..col_group2] > read_counts_cutoff) == length(col_group2)
  ]

  # build count matrix for running dds
  filtered_count_matrix <- filtered_count_table %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames("ids") %>% 
    as.matrix()

  # test are they consistent order between rownames of coldata & colnames of count matrix
  #if(any(rownames(col_data) != colnames(filtered_count_matrix))) stop("An error occurred.")

  # build DESeqDataSet
  dds <- DESeqDataSetFromMatrix(filtered_count_matrix, colData = col_data, design = ~ condition)  

  # run
  dds <- DESeq(dds)
  #dds$condition <- relevel(dds$condition, ref = treatment_conditions[2])
  res <- results(dds, contrast = c("condition", treatment_conditions))
  print(nrow(as.data.frame(res)))
  #resLFC <- lfcShrink(dds, coef = paste0("condition_", treatment_conditions[1], "_vs_", treatment_conditions[2])) # type = apeglm
  
  # tidy res
  restidy <- as.data.frame(res) %>% 
    tibble::rownames_to_column("ids") %>% 
    as.data.table() %>% 
    mutate(comparison_group_pair = paste0(treatment_conditions[1], "/", treatment_conditions[2])) %>% 
    left_join(count_table, by = "ids")

  # filter res | ordered by adj.p & keep sig. only
  filtered_res <- restidy %>% 
    na.omit() %>% 
    arrange(padj) %>% 
    .[padj < 0.05, ]

  # quanlity control
  se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                             colData=colData(dds))
  pcaData <- DESeq2::plotPCA(DESeqTransform(se), intgroup = c("condition"), returnData = TRUE)

  # plots
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1],"% variance")) +
    ylab(paste0("PC2: ", percentVar[2],"% variance")) + 
    ggrepel::geom_text_repel(label = pcaData$name, show.legend = FALSE) +
    coord_fixed() +
    theme_bw()
  p_hist_pval <- as.data.frame(restidy) %>% 
    ggplot() + 
    aes_string(x="pvalue") +
    geom_histogram() +
    theme_classic() +
    ggtitle("pvalue distribution")

  # write & save
  if(dir.exists(outdir)) {
    # save figs
    pdf(paste0(outdir, "/", ids_type, "_deseq2_QCplots_", treatment_conditions[1], "_vs_", treatment_conditions[2], ".pdf"))
      print(p_pca)
      print(p_hist_pval)
    dev.off()
    # save res table
    if(ids_type == "gene") {
      restidy_addcounts <- restidy %>% 
        left_join(unique(dt_ids[, c("gene_id", "gene_name")]), by = c("ids" = "gene_id")) %>% 
        dplyr::relocate(gene_name, comparison_group_pair, .after = ids)
      fwrite(restidy_addcounts, file = paste0(outdir, "/", opt$batch_id, ".", treatment_conditions[1], "_vs_", treatment_conditions[2], "_comparison.differential_gene_expression.tidy.txt"), sep = "\t")
    } else {
      restidy_addcounts <- restidy %>% 
        left_join(dt_ids, by = c("ids" = "isoform_id")) %>% 
        dplyr::relocate(gene_id, gene_name, comparison_group_pair, .after = ids) 
      fwrite(restidy_addcounts, file = paste0(outdir, "/", opt$batch_id, ".", treatment_conditions[1], "_vs_", treatment_conditions[2], "_comparison.differential_isoform_expression.tidy.txt"), sep = "\t")
    } 

  } else {
    stop("The outdir does not exists. Please check it.")
  }

  # return
  restidy_addcounts
}

plot_volcano <- function(
  tidy_res,  # tidy res generated from run_deseq2
  ids_type,  # gene or isoform
  log2foldchange_cutoff,  # minimum log2foldchange, twosides
  adj_p_value_cutoff,  # minimum adjusted p value
  outdir
) {

  comparison_group_pair_tag <- sub("/" , " vs. ", unique(tidy_res$comparison_group_pair))

  tidy_res <- tidy_res %>% 
    mutate(class = ifelse(log2FoldChange > 0 , "up", "down")) %>% 
    mutate(class = ifelse(
      abs(log2FoldChange) < log2foldchange_cutoff | 
      padj > adj_p_value_cutoff,
      "no",
      class
    ))
  
  if("gene_id" %in% colnames(tidy_res)) {
    tidy_res$label <- paste0(tidy_res$gene_name, ":", tidy_res$ids)
  } else {
    tidy_res$label <- tidy_res$gene_name
  }
  
  #tidy_res$label[tidy_res$class == "no"] <- NA
  
  volcanoplot_palette <- c("blue", "red", "grey")
  names(volcanoplot_palette) <- c("down", "up", "no")
  
  p <- ggplot(data = tidy_res,
         aes(
           x = log2FoldChange,
           y = -log10(padj),
           color = class,
         )) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_point(
      data = tidy_res %>% as.data.table() %>% .[gene_name %in% genes_to_label, ],
      shape = 21, stroke = 0.5, size = 0.8
    ) +
    geom_vline(
      xintercept = c(-log2foldchange_cutoff, log2foldchange_cutoff),
      col = "snow3",
      linetype = "dashed",
      linewidth = 0.3
    ) +
    geom_hline(
      yintercept = -log10(adj_p_value_cutoff),
      col = "snow3",
      linetype = "dashed",
      linewidth = 0.3
    ) +
    geom_text_repel(
      data = tidy_res %>% as.data.table() %>% arrange(pvalue) %>% .[1:20,],
      aes(label = label), 
      min.segment.length = 0.5, segment.size = 0.2,
      max.overlaps = 15, bg.color = "white",
      size = 1.2, max.time = 6
    ) +
    geom_label_repel(  # highlight the interested genes
      data = tidy_res %>% as.data.table() %>% .[gene_name %in% genes_to_label, ], 
      aes(label = label), size = 1.2, label.padding = 0.1, box.padding = 0.05
    ) + 
    annotate(
      geom = "text",
      x = -Inf, y = -log10(adj_p_value_cutoff), label = paste0("adj. p-value = ", adj_p_value_cutoff),
      vjust = 1.1, hjust = -0.1, size = 2, color = "black"
    ) +
    scale_x_continuous(limits = c(-6, 6), breaks = seq(-5, 5, 2)) +
    scale_color_manual(values = volcanoplot_palette) +
    labs(
      subtitle = comparison_group_pair_tag,
      x = "log2(fold change)",
      y = "-log10(adj. p-value)"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 7),
      axis.text = element_text(size = 7),
      panel.grid = element_blank(),
      axis.ticks = element_line(linewidth = 0.3),
      aspect.ratio = 1,
      legend.position = "none"
    )
  ggsave(
    paste0(outdir, "/", opt$batch_id, ".", treatment_conditions[1], "_vs_", treatment_conditions[2], "_comparison.differential_", ids_type, "_volcanoplot.pdf"),
    p,
    width = 2.5, height = 2.5
  )
}

# ============================================================================ #
# prep

# set output dir

base_dir <- getwd()
output_dir <- file.path(
  base_dir, 
  opt$batch_id, 
  "all_samples_combined", 
  paste0(opt$batch_id, "_differential_expression_output")
)
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

print("The plots & DEG table will be generated at the path: ")
print(output_dir)

# read the master table
master_table <- tryCatch({
  fread(opt$master_table)
}, error = function(err) {
  print("The master table does not exist.")
})

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

# ============================================================================ #
# prep | separate quant table to isoform table and gene table

dt_ids <- tidy_count_table[, c("isoform_id", "gene_id", "gene_name")] %>% unique()

gene_count_table <- tidy_count_table[, -c("isoform_id", "gene_name")] %>% 
  group_by(gene_id) %>% 
  summarise_all(sum) %>% 
  as.data.table()
setnames(gene_count_table, old = "gene_id", new = "ids")
#colnames(gene_count_table)[colnames(gene_count_table) == "gene_id"] <- "ids"

isoform_count_table <- tidy_count_table[, -c("gene_id", "gene_name")] %>% 
  group_by(isoform_id) %>% 
  summarise_all(sum) %>% 
  as.data.table()
setnames(isoform_count_table, old = "isoform_id", new = "ids")
#colnames(isoform_count_table)[colnames(isoform_count_table) == "isoform_id"] <- "ids"

# ============================================================================ #
# run deseq2 & save res table

res_gene    <- run_deseq2(gene_count_table,    opt$read_counts_cutoff, treatment_conditions, "gene",    output_dir)
res_isoform <- run_deseq2(isoform_count_table, opt$read_counts_cutoff, treatment_conditions, "isoform", output_dir)

# ============================================================================ #
# plot

# interested genes to highlight
if(length(opt$interested_genes) > 0) {
  genes_to_label <- unlist(strsplit(opt$interested_genes, split = ","))
} else {
  genes_to_label <- ""
}

plot_volcano(res_gene,    "gene",    opt$log2foldchange_cutoff, opt$adj_p_value_cutoff, output_dir)
plot_volcano(res_isoform, "isoform", opt$log2foldchange_cutoff, opt$adj_p_value_cutoff, output_dir)
