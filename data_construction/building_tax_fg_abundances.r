#pairing dada2 output with map.
rm(list=ls())
source('paths.r')

#set output path.----
output.path <- variance_analysis_product.path

#load data.----
esv <- readRDS(dada2_esv_table.path)
map <- read.table(map.path, header = T)
srr <- read.csv(srr_key.path)
tax <- readRDS(dada2_tax_table.path)

#remove leading x__ from taxonomy assignments.----
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}
colnames(tax) <- tolower(colnames(tax))

#merge map and srr info.----
srr$sample_name <- as.character(srr$sample_name)
srr$Name <- substr(srr$sample_name,6,nchar(srr$sample_name))
map <- merge(map, srr)
map$run <- as.character(map$run)

#drop observations that do not map to bacteria or archea from tax/esv tables.----
to_drop <- rownames(tax[is.na(tax[,1]),])
tax <- tax[!(rownames(tax) %in% to_drop),]
esv <- esv[,!(colnames(esv) %in% to_drop)]
tax <- data.frame(tax)

#get read depth.----
seq_depth <- rowSums(esv)

#start counting read by tax scale starting at phylum.----
tax.reads <- list()
for(i in 2:ncol(tax)){
  tax_level <- colnames(tax)[i]
  taxa <- unique(tax[,i])
  taxa.sum <- list()
  for(j in 1:length(taxa)){
    tax.sub <- tax[,i]
    if( is.na(taxa[j])){esv_cols <- tax.sub[is.na(tax.sub)]    } 
    if(!is.na(taxa[j])){esv_cols <- tax.sub[tax.sub %in% taxa[j]]}
    taxa.sum[[j]] <- rowSums(as.matrix(esv[,colnames(esv) %in% names(esv_cols)]))
  }
  taxa.sum <- do.call(cbind,taxa.sum)
  colnames(taxa.sum) <- taxa
  tax.reads[[i-1]] <- taxa.sum
  names(tax.reads)[i-1] <- tax_level
}


#Save for variance analysis.----
output <- list(map,tax.reads, seq_depth)
names(output) <- c('map','tax.reads','seq_depth')
saveRDS(output,output.path)
