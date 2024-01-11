# ============================================================================ #
# info

# date: 231116
# desc: characteristics of kmers
#         sequence logo for differentially modified k-mers 
# auth: yld
# proj: NanoTrans

# ============================================================================ #
# load

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggseqlogo))
suppressPackageStartupMessages(library(dplyr))

# ============================================================================ #
# parse args

option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character",
        default = NULL,
        help = "differentially modified k-mers table",
        metavar = "character"
    ),
    make_option(
        c("-s", "--output_seqlogo"),
        type = "character",
        default = NULL,
        help = "output name of the kmer seqlogo, type = pdf",
        metavar = "character"
    ),
    make_option(
        c("-f", "--output_frequency"),
        type = "character",
        default = NULL,
        help = "output name of the kmer seqlogo, type = pdf",
        metavar = "character"
    )   
)

opt <- parse_args(OptionParser(option_list = option_list))

# ============================================================================ #
# run

# read 
difmods <- fread(opt$input) 

# canonical m6A motifs
D <- c("G", "T", "A")
R <- c("G", "A")
H <- c("C", "T", "A")
df.motifs <- expand.grid(D, R, "A", "C", H)
m6a.motifs <- do.call(paste0, df.motifs[, 1:5])  # DRACH

# ============================================================================ #
# plot | seq logo 

kmer_seqlogo <- ggseqlogo(difmods[["kmer"]], method = 'prob') +
  ggplot2::scale_x_continuous(breaks = 1:5, labels = c(-2:2))

# ============================================================================ #
# plot | barplot for frequency of diff modified kmers

# select top10 freq kmers
freq_kmer <- difmods[, .N, kmer] %>% 
  mutate(ifClassic = ifelse(kmer %in% m6a.motifs, "DRACH", "Non-DRACH")) %>% 
  arrange(N) %>% 
  tail(10) %>% 
  mutate(kmer = factor(kmer, levels = kmer))

# plot
text.size = 8
bar_freqkmers <- ggplot(freq_kmer, aes(x = N, y = kmer, group = ifClassic)) +
    geom_col(aes(fill = ifClassic), width = 0.6) +
    labs(title = "Frequency of differential modified kmers") +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.15), add = c(0, 0)),
      position = "top"
    ) +
    scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
    geom_text(
      aes(x = N, y = kmer, group = ifClassic, label = N),
      hjust = -0.2,
      size = 3
    ) +
    scale_fill_manual(values = c("steelblue", "grey")) +
    theme(
      text = element_text(size = text.size),
      axis.text = element_text(size = text.size),
      panel.background = element_rect(fill = "white"),
      panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
      axis.ticks.length = unit(0, "mm"),
      axis.title = element_blank(),
      axis.line.y.left = element_line(color = "black"),
      legend.title = element_blank(),
      #legend.position = c(0.8, 0.2),
      #legend.background = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      aspect.ratio = 0.8
    )

# ============================================================================ #
# save plots

#pdf(file = opt$output, width = 4.5, height = 2.8)
#  kmer_seqlogo
#  bar_freqkmers
#dev.off()
pdf(file = opt$output_seqlogo,  width = 4, height = 3)
  kmer_seqlogo
dev.off()

pdf(file = opt$output_frequency, width = 4, height = 3)
  bar_freqkmers
dev.off()

