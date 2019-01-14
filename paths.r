#paths.r for dopeheide variance analysis.
master_dir <- '/projectnb/talbot-lab-data/caverill/dopheide_data'
host <- system('hostname', intern=T)
if(host == 'pecan2'){
  master_dir <- '/fs/data3/caverill/dopheide_data/'
}


#raw sequence data.
seq_dir <- paste0(master_dir,'raw_seq_data/')
cmd <- paste0('mkdir -p ',seq_dir)
system(cmd)
