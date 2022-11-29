#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
library(ggrepel)


option_list <- list(
  make_option(
    c("--input"),
    type = "character",
    default = NULL,
    help = "input summary table",
    metavar = "character"
  ),
  make_option(
    c("--group_a"),
    type = "character",
    default = NULL,
    help = "name of group a in comparison",
    metavar = "character"
  ),
  make_option(
    c("--group_b"),
    type = "character",
    default = NULL,
    help = "name of group b in comparison",
    metavar = "character"
  ),
  make_option(
    c("--prefix"),
    type = "character",
    default = NULL,
    help = "output prefix",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# # for testing only
# opt$input <- "Batch_TENT5A.all_samples_combined.polya_profiling.summary.txt"
# opt$group_a <- "TENT5A_WT"
# opt$group_b <- "TENT5A_MUT"
# opt$prefix <- "Batch_TENT5A.polyA_length"


data <- read.table(opt$input, header = TRUE, row.names = NULL, stringsAsFactors = FALSE, quote="", sep = "\t")
data$gene_isoform_tag <- paste0(data$gene_name, ":", data$isoform_id)
data$comparison_group <- factor(data$comparison_group, levels = c(opt$group_a, opt$group_b))

data_group_a <- subset(data, comparison_group == opt$group_a)
data_group_b <- subset(data, comparison_group == opt$group_b)


# based on median
median_length_wilcox_stat <- wilcox.test(data_group_a$median_length, data_group_b$median_length)
test_name <- median_length_wilcox_stat$method
p_value <- median_length_wilcox_stat$p.value
ggplot(data, aes(x=comparison_group, y=median_length, fill=comparison_group)) +
  geom_violin(alpha=0.3) +
  geom_boxplot(outlier.shape=NA, width=0.2) +
  ylab("median poly(A) tail length (bp)") +
  ggtitle(paste0("compared groups: ", opt$group_b, " vs. ", opt$group_a, "\ncompared quantity: median length\nstat test: ", test_name,"\nP=", p_value)) +
  theme_bw() 
ggsave(paste0(opt$prefix, ".median_length.violin_plot.pdf"))

# density plot
ggplot(data, aes(x=median_length, color=comparison_group)) +
  geom_density() +
  xlab("median poly(A) tail length (bp)") +
  ylab("density") +
  ggtitle(paste0("compared groups: ", opt$group_b, " vs. ", opt$group_a, "\ncompared quantity: median length")) +
  theme_bw() 
ggsave(paste0(opt$prefix, ".median_length.density_plot.pdf"))


# based on mean
mean_length_wilcox_stat <- wilcox.test(data_group_a$mean_length, data_group_b$mean_length)
test_name <- mean_length_wilcox_stat$method
p_value <- mean_length_wilcox_stat$p.value
ggplot(data, aes(x=comparison_group, y=mean_length, fill=comparison_group)) +
    geom_violin(alpha=0.3) +
    geom_boxplot(outlier.shape=NA, width=0.2) +
    ylab("mean poly(A) tail length (bp)") +
    ggtitle(paste0("compared groups: ", opt$group_b, " vs. ", opt$group_a, "\ncompared quantity: mean length\nstat test: ", test_name,"\nP=", p_value)) +
    theme_bw() 
ggsave(paste0(opt$prefix, ".mean_length.violin_plot.pdf"))
  
# density plot
ggplot(data, aes(x=mean_length, color=comparison_group)) +
  geom_density() +
  xlab("mean poly(A) tail length (bp)") +
  ylab("density") +	  
  ggtitle(paste0("compared groups: ", opt$group_b, " vs. ", opt$group_a, "\ncompared quantity: mean length")) +
  theme_bw() 
ggsave(paste0(opt$prefix, ".mean_length.density_plot.pdf"))
  
