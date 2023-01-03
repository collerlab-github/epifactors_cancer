setwd('../../Summary_tables/')

library(dplyr)


project <- 'TCGA' # TCGA or TARGET
if (project=='TCGA'){
  cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                   'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                   'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
} else if (project=='TARGET'){
  cancer_list <- c('AML','NBL','OS','WT')
} else{
  quit(save='no')
}

rank <- 2

################

# store median survival values
median_df <- as.data.frame(matrix(0,nrow=length(cancer_list),ncol=3))
colnames(median_df) <- c('OS','DSS','PFI')
rownames(median_df) <- cancer_list
# read in each cancer's survival curve p-values
for (c in cancer_list) {
  # file
  surv_file <- paste0('../',project,'_',c,'/04_clinical_analysis/Rank_',rank,'/',project,'_',c,'_survival_info.csv')
  # read and fill in cancer survival p-value
  surv_df <- read.csv(surv_file,row.names=1)
  median_df[c,'OS'] <- median(surv_df$OS.time, na.rm = T)
  median_df[c,'DSS'] <- median(surv_df$DSS.time,na.rm = T)
  median_df[c,'PFI'] <- median(surv_df$PFI.time,na.rm = T)
}

write.csv(median_df,paste0(project,'_median_survival_times.csv'))
