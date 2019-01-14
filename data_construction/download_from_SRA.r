#Downloading Tedersoo et al. 2014 ITS sequences from the SRA database.
#clear environment, source paths.
rm(list=ls())
source('paths.r')
library(SRAdb)

#output path.
output.dir <- seq_dir

#SRP project number.
SRP_id <- 'SRP148718'

#functions.----
tic = function() {assign("timer", Sys.time(), envir=.GlobalEnv)}
toc = function() print(Sys.time()-timer)

#Get SRA SQLite table. This takes a while.----
if(file.exists('SRAmetadb.sqlite') == T){
  srafile = 'SRAmetadb.sqlite'
  cat('SRA sqlite table already downloaded and unzipped. fantastic.\n')
}
if(file.exists("SRAmetadb.sqlite") == F){
  tic()
  srafile = getSRAdbFile()
  cat('SRA sqlite table downloaded and unzipped! ');toc()
}
con = dbConnect(RSQLite::SQLite(),srafile)

#get SRR codes associated with given SRP.----
conversion <- sraConvert(c(SRP_id), sra_con = con)

#Download the sequences.----
cat('Begin downloading sequences...\n'); tic()
getSRAfile(conversion$run, con, fileType = 'fastq',destDir = output.dir)
cat('Download complete.\n'); toc()

#unzip all fastq.gz files.----
cat('unzipping files...\n'); tic()
cmd <- paste0('gunzip ',output.dir,'/*.gz')
system(cmd)
cat('files unzipped! '); toc()

#clean up. Remove SRAdb sqlite file.----
system('rm SRAmetadb.sqlite')
cat('Script complete.\n')

#end script.----