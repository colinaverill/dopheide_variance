#Figuring out primer issues.
rm(list=ls())
source('paths.r')

#first, check for primers.
seq.files <- list.files(seq_dir)
seq.file <- seq.files[3]
fwd.primer <- 'GTGCCAGCMGCNGCGG'
rev.primer <- 'TCGTTG'   #just a fragment
fwd.primer.rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(fwd.primer)))
rev.primer.rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(rev.primer)))

#read lines.
#foo <- readLines(paste0(seq_dir,seq.file))
#grab <- seq(2, length(foo), by = 4)
#foo.dna <- foo[grab]
#test <- foo.dna[grep(rev.primer.rc, foo.dna)]

#lets run bbduk on one sample.
bbduk.path <- '/home/caverill/dopheide_variance/bbmap/bbduk.sh'
read1.fq.path <- paste0(seq_dir,seq.files[1])
read2.fq.path <- paste0(seq_dir,seq.files[2])

#quality filter fastq files using qiime.
sample.name <- seq.files[2]
sample.name <- substr(sample.name,1,nchar(sample.name)-8)
output.dir1 <- 'q.trim.L/'
output.path1 <- paste0(seq_dir,output.dir1,sample.name,'_1.fastq')
output.path2 <- paste0(seq_dir,output.dir1,sample.name,'_2.fastq')
cmd <- paste0(bbduk.path,
              ' literal=',fwd.primer,
              ' ktrim=l k=10',
              ' in=',read1.fq.path,
              #' in2=',read2.fq.path,
              ' out=',output.path1,
              #' out2=',output.path2,
              ' ordered=t')
system(cmd)
#now trim 3' end since all "forward" primers are 28bp long.
output.dir2 <- 'q.trim.R/'
sample.path <- paste0(seq.path,output.dir1,sample.name,'.fastq')
output.path <- paste0(seq.path,output.dir2,sample.name,'.fastq')
cmd <- paste0(bbduk.path,
              ' literal=',fwd.primer,
              ' ktrim=l k=10',
              ' in1=',read1.fq.path,
              ' in2=',read2.fq.path,
              ' out1=',output.path1,
              ' out2=',output.path2,
              ' ordered=t')
system(cmd)
