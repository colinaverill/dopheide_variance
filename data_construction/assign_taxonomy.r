#Assign taxonomy using dada2 in parallel.
#This script assumesyou have a taxonomy table where:
#1. the row names are sample names.
#2. the column names are the actual unique sequences.
#clear environment, source paths.
rm(list=ls())
source('paths.r')
#functions.----
tic = function() {assign("timer", Sys.time(), envir=.GlobalEnv)}
toc = function() print(Sys.time()-timer)

#specify output path.----
tax_output_path <- dada2_tax_table.path

#Load data.----
otu <- readRDS(dada2_esv_table.path)
to_assign <- colnames(otu) #grab sequences to assign.

#Everything from here below *should* just run and save where you told it to.
#download greengenes.----
greengenes.path <- paste0(small_dir,'gg_13_8_train_set_97.fa')
if(file.exists(greengenes.path)) {
  cat('using previously downloaded green genes database.')
} else {
  #download greengenes reference database.
  cat('downloading green genes...\n')
  gg.download.link <- 'https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz?download=1'
  cmd <- paste0('curl ',gg.download.link,' > ',greengenes.path)
  system(cmd)
  cat('greengenes download complete.\n')
}

#assign taxonomy.----
tic()
cat('Assigning taxonomy using the RDP Classifier...\n')
out <- dada2::assignTaxonomy(to_assign,greengenes.path,multithread = T, tryRC=T)
cat('Taxonomy assignment complete. ')
toc()

#remove greengenes download.----
cmd <- paste0('rm -f ',greengenes.path)
system(cmd)

#save output as your taxonomy file.----
saveRDS(out, tax_output_path)
