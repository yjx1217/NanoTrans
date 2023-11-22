# ============================================================================ #
# info

# date: 231011
# desc: barplot for frequency of various m6a motifs

# ============================================================================ #
# load

library(optparse)
library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)

# ============================================================================ #
# parse arguments

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "diffmod tidy table",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "output name",
    metavar = "character"
  )
)

opt_parser <- OptionParser(
    option_list = option_list,
    usage = "Barplot of m6a motifs' frequency"
)
opt <- parse_args(opt_parser)

# ============================================================================ #
# test 

if(FALSE) {
  opt <- list(
    input = "results/difmod-filtered-bydiffrate.tsv",
    output = "results/04c_barplot-freq-m6a-motifs.pdf"
  )
}

# ============================================================================ #
# read

difmod <- fread(opt$input, header = TRUE, fill = TRUE, quote = "", sep = "\t")

D <- c("G", "T", "A")
R <- c("G", "A")
H <- c("C", "T", "A")
df.motifs <- expand.grid(D, R, "A", "C", H)
m6a.motifs <- do.call(paste0, df.motifs[1:ncol(df.motifs)]) 

# ============================================================================ #
# prep plot data table

freq.kmer <- difmod[, .N, kmer] %>% 
  mutate(ifClassic = ifelse(kmer %in% m6a.motifs, "DRACH", "Non-DRACH")) %>% 
  arrange(N) %>% 
  tail(10) %>% 
  mutate(kmer = factor(kmer, levels = kmer))


# ============================================================================ #
# barplot for frequency of diff modified kmers

text.size = 8
bar.freq.kmers <- ggplot(freq.kmer, aes(x = N, y = kmer, group = ifClassic)) +
    geom_col(aes(fill = ifClassic), width = 0.6) +
    labs(title = "Frequency of differential modified kmers") +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.15), add = c(0, 0)),
      position = "top"
    ) +
    scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
    geom_text(
      #data = subset(data, count < 8),
      aes(x = N, y = kmer, group = ifClassic, label = N),
      hjust = -0.2,
      size = 3
    ) +
    scale_fill_manual(values = c("steelblue", "grey")) +
    theme(
      # Set the default size for all text
      text = element_text(size = text.size),
      axis.text = element_text(size = text.size),
      # Set background color to white
      panel.background = element_rect(fill = "white"),
      # Set the color and the width of the grid lines for the horizontal axis
      panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
      # Remove tick marks by setting their length to 0
      axis.ticks.length = unit(0, "mm"),
      # Remove the title for both axes
      axis.title = element_blank(),
      # Only left line of the vertical axis is painted in black
      axis.line.y.left = element_line(color = "black"),
      legend.title = element_blank(),
      legend.key.size = unit(0.4, "cm")
    )

# ============================================================================ #
# save 

ggsave(opt$output, bar.freq.kmers, width = 4.5, height = 2.8)
