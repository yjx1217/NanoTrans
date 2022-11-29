#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)


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
    c("--top_n"),
    type = "integer",
    default = NULL,
    help = "name of top entries to show",
    metavar = "integer"
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
# setwd("/Users/yjx/Projects/2022_NanoTrans/Project_Dataset2")
# opt$input <- "Batch_arabidopsis.rna_modification.majority_direction_kmer_diffmod.table.tidy.txt"
# opt$group_a <- "VIRc"
# opt$group_b <- "vir1"
# opt$top_n <- 20
# opt$prefix <- "Batch_arabidopsis.filtered"


data <- read.table(opt$input, header = TRUE, row.names = NULL, stringsAsFactors = FALSE, quote="", sep = "\t")
data_lower = subset(data, data$diff_mod_rate < 0)
data_higher = subset(data, data$diff_mod_rate > 0)

data_lower_top = head(data_lower, n=opt$top_n)
data_higher_top = head(data_higher, n=opt$top_n)


ggplot(data_lower_top, aes(x=gene_name, y=data_lower_top$diff_mod_rate, label = kmer)) +
  geom_point(aes(color = data_lower_top$pval)) +
  coord_cartesian(ylim = c(min(data_lower_top$diff_mod_rate)*1.1, max(data_lower_top$diff_mod_rate)*0.9), clip = "off") +
  # geom_text_repel() +
  geom_text(angle=45) +
  xlab("gene name") +
  ylab(paste0("differential modification rate: ", opt$group_a, "-", opt$group_b)) +
  scale_color_gradientn(limits=range(0, 1),
                       colours=c("red", "yellow"),
                       name="p_value", na.value="#101010") + 
  ggtitle(paste0("Top ", opt$top_n, " kmers with lower RNA modification rates in ", opt$group_a, " than ", opt$group_b)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
ggsave(paste0(opt$prefix, "diff_RNA_mod.", opt$group_a, "-", opt$group_b, ".negative.top", opt$top_n, ".pdf"))


ggplot(data_higher_top, aes(x=gene_name, y=data_higher_top$diff_mod_rate, label = kmer)) +
  geom_point(aes(color = data_higher_top$pval)) +
  coord_cartesian(ylim = c(min(data_higher_top$diff_mod_rate)*0.9, max(data_higher_top$diff_mod_rate)*1.1), clip = "off") +
  # geom_text_repel() +
  geom_text(angle=45) +
  xlab("gene name") +
  ylab(paste0("differential modification rate: ", opt$group_a, "-", opt$group_b)) +
  scale_color_gradientn(limits=range(0, 1),
                        colours=c("red", "yellow"),
                        name="p_value", na.value="#101010") + 
  ggtitle(paste0("Top ", opt$top_n, " kmers with higher RNA modification rates in ", opt$group_a, " than ", opt$group_b)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 
ggsave(paste0(opt$prefix, "diff_RNA_mod.", opt$group_a, "-", opt$group_b, ".positive.top", opt$top_n, ".pdf"))



