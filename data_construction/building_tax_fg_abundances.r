#pairing dada2 output with map.
rm(list=ls())
source('paths.r')

#load data.----
esv <- readRDS(dada2_esv_table.path)
map <- read.table(map.path, header = T)
srr <- read.csv(srr_key.path)

#merge map and srr info.
srr$sample_name <- as.character(srr$sample_name)
srr$Name <- substr(srr$sample_name,6,nchar(srr$sample_name))
map <- merge(map, srr)
map$run <- as.character(map$run)


