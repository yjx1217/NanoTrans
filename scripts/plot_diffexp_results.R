#!/usr/bin/env Rscript
library(optparse)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)


option_list <- list(
  make_option(
    c("--input_data"),
    type = "character",
    default = NULL,
    help = "input DEG tidy table",
    metavar = "character"
  ),
  make_option(
    c("--sample_a"),
    type = "character",
    default = NULL,
    help = "name of sample a in comparison",
    metavar = "character"
  ),
  make_option(
    c("--sample_b"),
    type = "character",
    default = NULL,
    help = "name of sample b in comparison",
    metavar = "character"
  ),
  make_option(
    c("--comparison_pair_by_file"),
    type = "character",
    default = NULL,
    help = "whether to extract the comparison pair information by the input file",
    metavar = "character"
  ),

  make_option(
    c("--log2foldchange_cutoff"),
    type = "double",
    default = "1",
    help = "log2foldchange cutoff",
    metavar = "character"
  ),
  make_option(
    c("--adj_p_value_cutoff"),
    type = "double",
    default = "0.05",
    help = "adjusted p-value cutoff",
    metavar = "character"
  ),
  make_option(
    c("--output_prefix"),
    type = "character",
    default = NULL,
    help = "output_prefix",
    metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
log2foldchange_bound <- 5


data <-
  read.table(
    opt$input_data,
    header = TRUE,
    row.names = NULL,
    stringsAsFactors = FALSE,
    quote="",
    sep = "\t"
  )


comparison_group_pair_tag <- "NA"
if (opt$comparison_pair_by_file == "yes") {
  comparison_group_pair_tag <- data$comparison_group_pair[[1]]
} else {
  comparison_group_pair_tag <- paste0(opt$sample_b, "/", opt$sample_a);
}

data_lite <- data[,1:7]
#data_lite$adj_p_value <- p.adjust(data$p_value, method="BH")



write.table(
  data_lite,
  file = paste0(opt$output_prefix, ".lite.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

data_lite_for_volcanoplot <-  data_lite[data_lite[["adj_p_value"]] != 9999, ]
data_lite_for_volcanoplot$class <- "no"
data_lite_for_volcanoplot$class[data_lite_for_volcanoplot$log2foldchange > opt$log2foldchange_cutoff &
                                   data_lite_for_volcanoplot$adj_p_value < opt$adj_p_value_cutoff] <-
  "up"
data_lite_for_volcanoplot$class[data_lite_for_volcanoplot$log2foldchange < -opt$log2foldchange_cutoff &
                                   data_lite_for_volcanoplot$adj_p_value < opt$adj_p_value_cutoff] <-
  "down"

if(all(is.na(data_lite_for_volcanoplot[["isoform_id"]]))) {
  data_lite_for_volcanoplot$label <- data_lite_for_volcanoplot$gene_name
} else {
  data_lite_for_volcanoplot$label <-
    paste0(data_lite_for_volcanoplot$gene_name, ":", data_lite_for_volcanoplot$isoform_id)
}

data_lite_for_volcanoplot$label[data_lite_for_volcanoplot$class == "no"] <-
  NA

volcanoplot_palette <- c("blue", "red", "grey")
names(volcanoplot_palette) <- c("down", "up", "no")

p <- ggplot(data = data_lite_for_volcanoplot,
       aes(
         x = log2foldchange,
         y = -log10(adj_p_value),
         col = class,
         label = label
       )) +
  geom_point(alpha = 0.3) +
  geom_text_repel(size = 2) +
  geom_vline(
    xintercept = c(-opt$log2foldchange_cutoff, opt$log2foldchange_cutoff),
    col = "snow3",
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = -log10(opt$adj_p_value_cutoff),
    col = "snow3",
    linetype = "dashed",
    linewidth = 0.3
  ) +
  annotate(
    geom = "text",
    x = -Inf, y = -log10(opt$adj_p_value_cutoff), label = paste0("adj. p-value = ", opt$adj_p_value_cutoff),
    vjust = 1.1, hjust = -0.1, size = 2.5, color = "black"
  ) +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 2)) +
  scale_color_manual(values = volcanoplot_palette) +
  labs(
    subtitle = comparison_group_pair_tag,
    x = "log2(fold change)",
    y = "-log10(adj. p-value)"
  ) +
  #ggtitle(paste0("foldchange calculation: ", comparison_group_pair_tag, "\nadj_p_value_cutoff: ", opt$adj_p_value_cutoff, "\nlog2foldchange_cutoff: ", opt$log2foldchange_cutoff, ",", -opt$log2foldchange_cutoff)) +
  theme_bw() +
  theme(
    text = element_text(size = 9),
    panel.grid = element_blank(),
    axis.ticks = element_line(linewidth = 0.3),
    aspect.ratio = 1,
    legend.position = "none"
  )
ggsave(paste0(opt$output_prefix, ".volcanoplot.pdf"), p, width = 4, height = 4)



