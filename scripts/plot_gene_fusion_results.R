#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
library(grid)
library(dplyr)   

option_list <- list(
  make_option(
    c("--input"),
    type = "character",
    default = NULL,
    help = "input gene fusion table",
    metavar = "character"
  ),
  make_option(
    c("--gtf"),
    type = "character",
    default = NULL,
    help = "input gtf file",
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
# opt$input <- "SGNex_A549_rep1.fusion_gene.transcripts.tsv"
# opt$gtf <- "ref.genome.gtf"
# opt$prefix <- "Batch_SGNex.gene_fusion"

gtf_data <- read.table(opt$gtf, header = FALSE, row.names = NULL, stringsAsFactors = FALSE, quote="", sep = "\t")
colnames(gtf_data) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
gtf_data$gene_id <- sapply(as.character(gtf_data$attribute), function(x) sub(".*gene_id\\s\"([^;]+)\";.*", "\\1", x))

fusion_gene_data <- read.table(opt$input, header = TRUE, row.names = NULL, stringsAsFactors = FALSE, quote="", sep = "\t")


########## multiplot #############
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=2, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

for (r in 1:nrow(fusion_gene_data)) {
  fusion_part_by_name <- unlist(strsplit(fusion_gene_data$fusion_gene_by_name[r], ":"))
  fusion_part_by_id <- unlist(strsplit(fusion_gene_data$fusion_gene_by_id[r], ":"))
  
  part1_name <- fusion_part_by_name[1]
  part1_id <- fusion_part_by_id[1]
  part1_chr <- fusion_gene_data$part1_chr[r]
  part1_breakpoint <- fusion_gene_data$part1_breakpoint[r]
  part1_strand <- fusion_gene_data$part1_strand[r]
  
  part2_name <- fusion_part_by_name[2]
  part2_id <- fusion_part_by_id[2]
  part2_chr <- fusion_gene_data$part2_chr[r]
  part2_breakpoint <- fusion_gene_data$part2_breakpoint[r]
  part2_strand <- fusion_gene_data$part2_strand[r]
  
  part1_gtf <- subset(gtf_data, gene_id == part1_id)
  part1_gtf$transcript_id <- sapply(as.character(part1_gtf$attribute), function(x) sub(".*transcript_id\\s\"([^;]+)\";.*", "\\1", x))
  part1_gtf$fusion_tag <- "part1"

  part2_gtf <- subset(gtf_data, gene_id == part2_id)
  part2_gtf$transcript_id <- sapply(as.character(part2_gtf$attribute), function(x) sub(".*transcript_id\\s\"([^;]+)\";.*", "\\1", x))
  part2_gtf$fusion_tag <- "part2"
  
  
  part1_structure <- part1_gtf[, c(1, 3:5, 7, 11, 12)]
  names(part1_structure) <- c("chr", "feature", "start", "end", "strand", "transcript_id", "fusion_tag")
  part1_idx <- which(part1_structure$feature == "transcript")
  part1_s <- part1_idx + 1
  part1_e <- c(part1_idx[-1]-1, nrow(part1_structure))
  part1_g <- lapply(seq_along(part1_s), function(i) {
    x <- part1_structure[part1_s[i]:part1_e[i],]
    return(x)
    }) %>% do.call(rbind, .)
  
  part1_exons <- part1_g[part1_g$feature == "exon", ]
  part1_gene <- part1_structure[part1_structure$feature == "gene", ]
  part1_transcripts <- part1_structure[part1_structure$feature == "transcript",]
  part1_fused_region <- part1_gene
  part1_fused_region$feature <- "fused_region"
  part1_fused_region$end <- part1_breakpoint
  
  part2_structure <- part2_gtf[, c(1, 3:5, 7, 11, 12)]
  names(part2_structure) <- c("chr", "feature", "start", "end", "strand", "transcript_id", "fusion_tag")
  part2_idx <- which(part2_structure$feature == "transcript")
  part2_s <- part2_idx + 1
  part2_e <- c(part2_idx[-1]-1, nrow(part2_structure))
  part2_g <- lapply(seq_along(part2_s), function(i) {
    x <- part2_structure[part2_s[i]:part2_e[i],]
    return(x)
  }) %>% do.call(rbind, .)
  
  part2_exons <- part2_g[part2_g$feature == "exon", ]
  part2_gene <- part2_structure[part2_structure$feature == "gene", ]
  part2_transcripts <- part2_structure[part2_structure$feature == "transcript",]
  part2_fused_region <- part2_gene
  part2_fused_region$feature <- "fused_region"
  part2_fused_region$end <- part2_breakpoint
  
  part1_transcript_arrow_direction <- "last"
  if (part1_strand == "-") {
    part1_transcript_arrow_direction <- "first"
  }
  part2_transcript_arrow_direction <- "last"
  if (part2_strand == "-") {
    part2_transcript_arrow_direction <- "first"
  }
part1_plot <- ggplot(data=part1_transcripts) +
    geom_rect(data=part1_fused_region, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="lightgrey", alpha=0.5) +
    geom_segment(aes(x=start, xend=end, y=transcript_id, yend=transcript_id, color=transcript_id), arrow = arrow(length = unit(0.15, "cm"), ends=part1_transcript_arrow_direction, type = "closed"), size=1) + 
    geom_segment(data=part1_exons, aes(x=start, xend=end, y=transcript_id, yend=transcript_id, color=transcript_id, alpha = 0.6), size=6) + 
    xlab(paste0("chromosome:", part1_chr, " genomic coordinate (bp)")) + 
    ylab("transcript ID") +
    labs(title=paste0(part1_name, "|", part1_id)) + 
    theme_bw() +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1)) 
part2_plot <- ggplot(data=part2_transcripts) +
  geom_rect(data=part2_fused_region, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="lightgrey", alpha=0.5) +
  geom_segment(aes(x=start, xend=end, y=transcript_id, yend=transcript_id, color=transcript_id), arrow = arrow(length = unit(0.15, "cm"), ends=part2_transcript_arrow_direction, type = "closed"), size=1) + 
  geom_segment(data=part2_exons, aes(x=start, xend=end, y=transcript_id, yend=transcript_id, color=transcript_id, alpha = 0.6), size=6) + 
  xlab(paste0("chromosome:", part2_chr, " genomic coordinate (bp)")) + 
  ylab("transcript ID") +
  labs(title=paste0(part2_name, "|", part2_id)) + 
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1)) 

pdf(paste0(opt$prefix, ".", part1_name, "_", part2_name, ".gene_fusion_plot.pdf"), width=8, height=6)
print(multiplot(part1_plot, part2_plot))
dev.off()
  
  
}

