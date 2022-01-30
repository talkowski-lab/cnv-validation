#!/usr/bin/env R
# Author: asanchis@broadinstitute.org

library(data.table)
library(optparse)

option_list = list(
  make_option(c("-f", "--files"), type="character", default=NULL,
              help="files to merge", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

files <- fread(opt$files, header=F)
out <- opt$out

merged <- do.call(rbind, lapply(files$V1, fread))

merged_sort <- merged[order(merged$CHROM, merged$START, merged$END),]

write.table(merged_sort, out, sep="\t", quote=F, row.names=F)
