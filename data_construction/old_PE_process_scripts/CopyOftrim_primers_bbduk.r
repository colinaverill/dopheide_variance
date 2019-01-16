#Subset to 16S sequences, trim primers.
#Note bbduk requires java and only works on the scc.
rm(list=ls())
source('paths.r')

#declare primer sequences.----
fwd.primer <- 'GTGCCAGCNGCNGCGG' #first N is actually a M, but bbsuk only handles N.
rev.primer <- 'GGGTTNCGNTCGTTG'

#Setup intermediate directories.----
#seq_dir is defined in paths.r
sub_16S.fwd <- paste0(seq_dir,'sub_16S.fwd/')
sub_16S.rev <- paste0(seq_dir,'sub_16S.rev/')
trim.L <- paste0(seq_dir,'trim.L/')
trim.R <- paste0(seq_dir,'trim.R/')
dir.list <- list(sub_16S.fwd,sub_16S.rev,trim.L,trim.R)
for(i in 1:length(dir.list)){system(paste0('mkdir -p ',dir.list[[i]]))}

#Set paths to bbmap functions, sequence files and sample names.----
   bbduk.path <- '/projectnb/talbot-lab-data/caverill/dopheide_variance/bbmap/bbduk.sh'
reformat.path <- '/projectnb/talbot-lab-data/caverill/dopheide_variance/bbmap/reformat.sh'
seq.files <- list.files(seq_dir)
seq.files <- seq.files[grep('.fastq',seq.files)]
sample.names <- substr(seq.files,1,nchar(seq.files) - 8)
sample.names <- unique(sample.names)


#get all permutations fwd-rev primers. This was for troubleshooting.
all_bases = c('A', 'T', 'C', 'G')
all.fwd <- data.table::tstrsplit(fwd.primer, '', fixed = TRUE)
all.fwd <- lapply(all.fwd, function(x) if (x == 'N') all_bases else x)
all.fwd <- Reduce(paste0, do.call(data.table::CJ, all.fwd))
all.rev <- data.table::tstrsplit(rev.primer, '', fixed = TRUE)
all.rev <- lapply(all.rev, function(x) if (x == 'N') all_bases else x)
all.rev <- Reduce(paste0, do.call(data.table::CJ, all.rev))
all.fwd.rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(all.fwd)))
all.rev.rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(all.rev)))

#read lines, 'manually' look for primers. This was also for troubleshooting.
#Read 1 is 5'->3', starts with reverse primer.
#Read 2 is 5'->3', starts with forward primer.
#minimum length DNA sequence observed is 35 characters, max is 301.
#This means we should left trim each primer out of fwd/rev pairs before downstream analysis.
#foo <- readLines(paste0(seq_dir,seq.files[1]))
#grab <- seq(2, length(foo), by = 4)
#foo.dna <- foo[grab]
#test <- foo.dna[grep(paste(all.rev,collapse="|"), foo.dna)];length(test)

#Match reads to 16S primers based on forward read.----
output.dir <- paste0(basename(sub_16S.fwd),'/')
for(i in 1:length(sample.names)){
  read1.fq.path <- paste0(seq_dir,sample.names[i],'_1.fastq')
  read2.fq.path <- paste0(seq_dir,sample.names[i],'_2.fastq')
  sample.name <- sample.names[i]
  output.path1 <- paste0(seq_dir,output.dir,sample.name,'_1.fastq')
  output.path2 <- paste0(seq_dir,output.dir,sample.name,'_2.fastq')
  unmatch.path <- paste0(seq_dir,output.dir,sample.name,'_unmatch.fastq')
  cmd <- paste0(bbduk.path,
                ' literal=',rev.primer,
                ' k=10 copyundefined mm=f',
                ' in1=',read1.fq.path,
                ' in2=',read2.fq.path,
                ' out1=',output.path1,
                ' out2=',output.path2,
                ' outm=',unmatch.path,
                ' ordered=t skipr2')
  system(cmd)
  #everything we actually want is interleaved within the unmatch file. fix this.
  system(paste0('rm -f ',output.path1))
  system(paste0('rm -f ',output.path2))
  cmd <- paste0(reformat.path,
                ' in=',unmatch.path,
                ' out1=',output.path1,
                ' out2=',output.path2)
  system(cmd)
  system(paste0('rm -f ',unmatch.path))
}

#Match reads to 16S primers based on reverse read.----
ref_dir <- paste0(seq_dir,output.dir)
output.dir <- paste0(basename(sub_16S.rev),'/')
for(i in 1:length(sample.names)){
  read1.fq.path <- paste0(ref_dir,sample.names[i],'_1.fastq')
  read2.fq.path <- paste0(ref_dir,sample.names[i],'_2.fastq')
  sample.name <- sample.names[i]
  output.path1 <- paste0(seq_dir,output.dir,sample.name,'_1.fastq')
  output.path2 <- paste0(seq_dir,output.dir,sample.name,'_2.fastq')
  unmatch.path <- paste0(seq_dir,output.dir,sample.name,'_unmatch.fastq')
  cmd <- paste0(bbduk.path,
                ' literal=',fwd.primer,
                ' k=10 copyundefined mm=f',
                ' in1=',read1.fq.path,
                ' in2=',read2.fq.path,
                ' out1=',output.path1,
                ' out2=',output.path2,
                ' outm=',unmatch.path,
                ' ordered=t skipr1')
  system(cmd)
  #everything we actually want is interleaved within the unmatch file. fix this.
  system(paste0('rm -f ',output.path1))
  system(paste0('rm -f ',output.path2))
  cmd <- paste0(reformat.path,
                ' in=',unmatch.path,
                ' out1=',output.path1,
                ' out2=',output.path2)
  system(cmd)
  system(paste0('rm -f ',unmatch.path))
}


#Left trim reverse primer from forward read.----
ref_dir <- paste0(seq_dir,output.dir)
output.dir <- paste0(basename(trim.L),'/')
for(i in 1:length(sample.names)){
  read1.fq.path <- paste0(ref_dir,sample.names[i],'_1.fastq')
  read2.fq.path <- paste0(ref_dir,sample.names[i],'_2.fastq')
  output.path1 <- paste0(seq_dir,output.dir,sample.names[i],'_1.fastq')
  output.path2 <- paste0(seq_dir,output.dir,sample.names[i],'_2.fastq')
  cmd <- paste0(bbduk.path,
                ' literal=',rev.primer,
                ' ktrim=l k=10 copyundefined mm=f',
                ' in1=',read1.fq.path,
                ' in2=',read2.fq.path,
                ' out1=',output.path1,
                ' out2=',output.path2,
                ' ordered=t skipr2')
  system(cmd)
}

#Left trim forward primer from reverse read.----
ref_dir <- paste0(seq_dir,output.dir)
output.dir <- paste0(basename(trim.R),'/')
for(i in 1:length(sample.names)){
  read1.fq.path <- paste0(ref_dir,sample.names[i],'_1.fastq')
  read2.fq.path <- paste0(ref_dir,sample.names[i],'_2.fastq')
  output.path1 <- paste0(seq_dir,output.dir,sample.names[i],'_1.fastq')
  output.path2 <- paste0(seq_dir,output.dir,sample.names[i],'_2.fastq')
  cmd <- paste0(bbduk.path,
                ' literal=',fwd.primer,
                ' ktrim=l k=10 copyundefined mm=f',
                ' in1=',read1.fq.path,
                ' in2=',read2.fq.path,
                ' out1=',output.path1,
                ' out2=',output.path2,
                ' ordered=t skipr1')
  system(cmd)
}

#You all trimmed up! move files and cleanup.----
#move final file up a directory and rename to "trim_fastq_files"
cmd <- paste0('mv ',trim.R,' ',trim_dir)
system(cmd)
#remove intermediate files.
system(paste0('rm -rf ',sub_16S.fwd))
system(paste0('rm -rf ',sub_16S.rev))
system(paste0('rm -rf ',trim.L))

#end script.----
cat('16S subsetting and primer trimminmg complete.\n')

