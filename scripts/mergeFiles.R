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

print("Defining variable")
files <- fread(opt$files, header=F)
out <- opt$out

print("Merging")
merged <- do.call(merge, lapply(file, fread))

print("Writing")
write.table(merged, out, sep="\t", quote=F, row.names=F)
