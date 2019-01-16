#paths.r for dopeheide variance analysis.
master_dir <- '/projectnb/talbot-lab-data/caverill/dopheide_data/'
host <- system('hostname', intern=T)
if(host == 'pecan2'){
  master_dir <- '/fs/data3/caverill/dopheide_data/'
}

#raw sequence data.----
 seq_dir <- paste0(master_dir,'raw_seq_data/')
trim_dir <- paste0(master_dir,'trim_fastq_files/') #post bbduk trimming.
cmd <- paste0('mkdir -p ',seq_dir)
system(cmd)

#small data directory.----
small_dir <- paste0(master_dir,'small_data/')
cmd <- paste0('mkdir -p ',small_dir)
system(cmd)

#dada2 output----
        dada2_esv_table.path <- paste0(small_dir,'dada2_esv_table.rds')
dada2_sequence_tracking.path <- paste0(small_dir,'dada2_sequence_tracking.rds')
        dada2_tax_table.path <- paste0(small_dir,'dada2_tax_table.rds')
