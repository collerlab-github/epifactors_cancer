library(dplyr)


project <- 'TARGET' # TCGA or TARGET
if (project=='TCGA'){
  cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                   'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                   'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
  end <- 'PFI'
} else if (project=='TARGET'){
  cancer_list <- c('AML','NBL','OS','WT')
  end <- 'OS'
} else{
  quit(save='no')
}################

measures_df <- data.frame('cancer'=c(),'rank'=c(), 'silhouette'=c(),'cophenetic'=c()) 

for (c in cancer_list){
  measures <- read.csv(paste0('../../',project,'_',c,'/03_nmf/multi_rank_measures.csv'))
  measures$cancer <- c
  measures_df <- rbind(measures_df,measures[,c('cancer','rank','silhouette.consensus','cophenetic')])
}
colnames(measures_df)[3] <- 'silhouette'

write.csv(measures_df,paste0('../../Summary_tables/',project,'_nmf_cluster_metrics.csv'))

