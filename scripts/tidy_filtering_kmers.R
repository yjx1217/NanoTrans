# ============================================================================ #
# info

# date: 231019
# desc: filter diffmods
# auth: yld
# proj: NanoTrans

# ============================================================================ #
# load

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

# ============================================================================ #
# parse arguments

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = NULL,
    help = "diffmod table",
    metavar = "character"
  ),
  make_option(
    c("-f", "--fdirection"),
    type = "character",
    default = NULL,
    help = "direction of filtering diff mod rate"
  ),
  make_option(
    c("-p", "--pval"),
    type = "numeric",
    default = NULL,
    help = "cutoff of pvalue",
    metavar = "numeric"
  ),
  make_option(
    c("-d", "--difrate"),
    type = "numeric",
    default = NULL,
    help = "cutoff of diff mod rate",
    metavar = "numeric"
  ),
  make_option(
    c("-l", "--lowfreq"),
    type = "integer",
    default = NULL,
    help = "cutoff of kmers frequency",
    metavar = "integer"
  ),
  make_option(
    c("-o", "--outtsv"),
    type = "character",
    default = NULL,
    help = "output name",
    metavar = "character"
  )
)

opt_parser <- OptionParser(
    option_list = option_list,
    usage = "Filtering differentially modified kmers"
)
opt <- parse_args(opt_parser)

if(FALSE) {
    opt <- list(
       #input = "results/arabidopsis/Batch_Dataset2.rna_modification.majority_direction_kmer_diffmod.table.tidy.txt"
       input =  "results/hek293/Batch_Dataset3.rna_modification.majority_direction_kmer_diffmod.table.tidy.txt",
       fdirection = "twoside",  # higher, lower, twoside
       outtsv = "results/hek293/difmod-filtered-bydiffrate.tsv"
    )
}

# ============================================================================ #
# vars

difmods <- fread(opt$input) 
colnames.difmod <- colnames(difmods)
colnames.new <- colnames.difmod %>% 
  sub("^diff_mod_rate.*", "diff_mod_rate", ., perl = TRUE) %>% 
  sub("^pval.*", "pval", ., perl = TRUE) %>% 
  sub("^z_score.*", "z_score", ., perl = TRUE)
colnames(difmods) <- colnames.new

# motifs
D <- c("G", "T", "A")
R <- c("G", "A")
H <- c("C", "T", "A")
df.motifs <- expand.grid(D, R, "A", "C", H)
m6a.motifs <- do.call(paste0, df.motifs[, 1:5])  # DRACH

# ============================================================================ #
# filtered only by diff mod rate

p <- opt$pval
lowfreq <- opt$lowfreq
difrate <- opt$difrate

# diff modified direction
if(opt$fdirection == "twoside") {
  difmods.f <- difmods[pval <= p] %>% 
    .[abs(diff_mod_rate) >= difrate] 
} else if(opt$fdirection == "higher") {
  difmods.f <- difmods[pval <= p] %>% 
    .[diff_mod_rate >= difrate] 
} else if(opt$fdirection == "lower") {
  difmods.f <- difmods[pval <= p] %>% 
    .[diff_mod_rate <= -difrate]
} else {
  stop("Arg -f is invalid!")
}

lowfreq.kmers <- difmods.f[, .N ,kmer][N <= lowfreq, kmer]
difmods.f <- difmods.f[grep("^..A", kmer)] %>% 
  .[! kmer %in% lowfreq.kmers] 

fwrite(difmods.f, opt$outtsv, sep = "\t")

