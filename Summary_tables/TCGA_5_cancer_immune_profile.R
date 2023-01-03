library(dplyr)
library(readxl)
library(writexl)

immune <- read_excel('../../metadata/TCGA_immune_profile.xlsx')

cancers <- c('ACC','KIRC','LGG','LIHC','LUAD')
top5_df <- data.frame('id'=c(),'cancer'=c(),'cluster'=c())
for (cancer in cancers){
  info <- read_excel('../../nmf_clinical_analysis_figures/top_five_cancer_survival_info.xlsx', sheet = paste0('TCGA_',cancer,'_survival_info'))
  info$cancer <- cancer
  top5_df <- rbind(top5_df, info[,c(1,11,8)])
}
colnames(top5_df)[1] <- 'id'
top5_immune <- merge(top5_df, immune, by.x=1, by.y=1)

write_xlsx(top5_immune, '../../nmf_clinical_analysis_figures/top_five_immune_profile.xlsx')
