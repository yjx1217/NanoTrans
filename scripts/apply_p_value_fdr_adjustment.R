#!/usr/bin/env Rscript

# library(optparse)

# option_list <- list(
#   make_option(
#     c("--input"),
#     type = "character",
#     default = NULL,
#     help = "input tsv table",
#     metavar = "character"
#   ),
#   make_option(
#     c("--output"),
#     type = "character",
#     default = NULL,
#     help = "input tsv table",
#     metavar = "character"
#   )
# )


# opt_parser <- OptionParser(option_list = option_list)
# opt <- parse_args(opt_parser)

# data <-
#   read.table(
#     opt$input_data,
#     header = TRUE,
#     row.names = NULL,
#     stringsAsFactors = FALSE,
#     quote="",
#     sep = "\t"
#   )


# data$adj_p_value <- p.adjust(data$p_value, method = "BH")

# write.table(
#   data,
#   file = opt$output,
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE
# )

args <- commandArgs(trailingOnly = TRUE)
data <-
  read.table(
    args[1],
    header = TRUE,
    row.names = NULL,
    stringsAsFactors = FALSE,
    quote="",
    sep = "\t"
  )


data$adj_p_value <- p.adjust(data$p_value, method = "BH")

write.table(
  data,
  file = args[2],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
