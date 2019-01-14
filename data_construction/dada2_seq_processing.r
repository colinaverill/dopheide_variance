#process dopheide sequences using dada2.
rm(list=ls())
library(dada2)
source('paths.r')
#functions.----
tic = function() {assign("timer", Sys.time(), envir=.GlobalEnv)}
toc = function() print(Sys.time()-timer)


#set output path.----
output_filepath <- dada2_esv_table.path
 track_filepath <- dada2_sequence_tracking.path

#path to sequences
seq.dir <- seq_dir

# set fwd reverse truncaltion length based on quality profile plots.
truncation_setting <- c(250, 150)

#grab fwd and reverse reads, sample names.----
fwd <- sort(list.files(seq.dir, pattern = '_1.fastq', full.names = T))
rev <- sort(list.files(seq.dir, pattern = '_2.fastq', full.names = T))
sample.names <- substr(basename(fwd),1, nchar(basename(fwd)[1]) - 8)

#Filter samples.----
cat('Begin quality filtering...\n')
tic()
filtFs <- file.path(seq.dir,'filtered',basename(fwd))
filtRs <- file.path(seq.dir,'filtered',basename(rev))
out <- filterAndTrim(fwd, filtFs, rev, filtRs, maxN = 0, maxEE = c(2,2), truncQ = 2, minLen = 50, rm.phix = TRUE, 
                     truncLen = truncation_setting, compress = TRUE, multithread = T)
cat('Quality filtering complete.\n')
toc()

#Learn error rates.----
cat('Learning error rates...\n');tic()
errF <- learnErrors(filtFs, multithread = T)
errR <- learnErrors(filtRs, multithread = T)
cat('Error rates learned.\n');toc()

#Deprelicate reads.----
tic()
cat('Dereplicating reads...\n')
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
cat('Reads dereplicated.\n')
toc()

#Sample inference with dada2 algorithm.----
tic()
cat('Performing sample inference with dada2 algorithm...\n')
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
cat('Sample inference complete.\n')
toc()

#Merge paired reads, make sequence table.----
tic()
cat('Merging forward and reverese reads...\n')
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = T)
seqtab <- makeSequenceTable(mergers)
cat('Reads merged.\n')
toc()
dim(seqtab)

#Remove chimeras.----
tic()
cat('Removing chimeras...\n')
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
cat('Chimeras removed.\n')
toc()

#Track reads through pipeline.----
getN <- function(x) sum(getUniques(x))
out_sub <- out[rownames(out) %in% paste0(sample.names,'_1.fastq'),]
track <- cbind(out_sub, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names
head(track)

#clean up intermediate files.----
cmd <- paste0('rm -rf ',seq_dir,'filtered/')
system(cmd)

#Save output.----
saveRDS(seqtab.nochim, output_filepath)
saveRDS(track        ,  track_filepath)
cat('Analysis complete.\n')
#end script.

