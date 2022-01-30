#!/usr/bin/env R
# Author: asanchis@broadinstitute.org

library(data.table)
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
)
 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

df <- fread(opt$input)
samples <- fread(opt$samples, header = F)
out <- opt$out
chrom <- unique(df$chr)

df2 <- subset(df, sample %in% samples$V1)

df2 %>% 
	group_by(id) %>% 
	summarise(chr = chr, start = start, end = end, id = id, svtype = svtype, samples = paste0(unique(sample), collapse=",")) %>%
	unique() -> df3

write.table(df3, out, sep="\t", quote=F, row.names=F, col.names=F)
