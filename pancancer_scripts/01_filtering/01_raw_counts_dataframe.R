
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<=1){
  stop('Usage: Rscript 01_raw_counts_dataframe.R [cancer1] [cancer2] ...\nAt least 2 cancers must be supplied (i.e.)')
} else{
  print(args)
}

project <- 'TARGET'
dir.create(paste0('../../',project,'_PAN'), showWarnings = F)
setwd(paste0('../../',project,'_PAN'))
dir.create('01_filtering', showWarnings = F)


# create combined raw_counts df and patient list
patient_df <- data.frame(type=c())
pan_df <- data.frame()
for (cancer in args){
  print(paste0('appending ', cancer))
  cfile <- paste0('../',project,'_',cancer,'/RNA-seq_datasets/',cancer,'_raw_counts.csv')
  df <- read.csv(cfile,row.names=1)
  p_names <- data.frame(colnames(df), type=cancer, row.names=1)
  patient_df <- rbind(patient_df, p_names)
  if ((nrow(pan_df)==0)|(ncol(pan_df)==0)){
    pan_df <- df
  } else{
    pan_df <- merge(pan_df,df,by=0)
    rownames(pan_df) <- pan_df[,1]
    pan_df <- pan_df[,-1]
  }
}

# make rnaseq counts directory
dir.create("RNA-seq_datasets")
# write raw counts table
write.csv(pan_df,file=paste0("RNA-seq_datasets/",paste(args,collapse='_'),"_raw_counts.csv"),quote=F)
write.csv(patient_df, file=paste0('01_filtering/', paste(args,collapse='_'),'_patient_data.csv'),quote=F)
