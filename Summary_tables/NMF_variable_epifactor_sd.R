library(dplyr)
library(ggplot2)


# cancers <- c('ACC','LGG','LUAD')
TCGA_cancers <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                  'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                  'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
TARGET_cancers <- c('AML','NBL','OS','WT')


# TCGA
cancer_df <- data.frame()

for (cancer in TCGA_cancers){
  epigene_file <- list.files(path=paste0('../../TCGA_',cancer,'/RNA-seq_datasets/'),pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'), full.names = T)
  epigene_df <- read.csv(epigene_file, row.names=1)
  epigene_sds <- data.frame(apply(epigene_df,1,sd))
  epigene_sds$gene <- rownames(epigene_sds)
  epigene_sds$cancer <- cancer
  epigene_sds <- epigene_sds[,c(2,3,1)]
  colnames(epigene_sds)[3] <- 'sd'
  rownames(epigene_sds) <- NULL
  cancer_df <- rbind(cancer_df, epigene_sds)
}

write.csv(cancer_df, '../../Summary_tables/TCGA_NMF_input_variable_epifactors.csv', row.names=F)

# TARGET
cancer_df <- data.frame()

for (cancer in TARGET_cancers){
  epigene_file <- list.files(path=paste0('../../TARGET_',cancer,'/RNA-seq_datasets/'),pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'), full.names = T)
  epigene_df <- read.csv(epigene_file, row.names=1)
  epigene_sds <- data.frame(apply(epigene_df,1,sd))
  epigene_sds$gene <- rownames(epigene_sds)
  epigene_sds$cancer <- cancer
  epigene_sds <- epigene_sds[,c(2,3,1)]
  colnames(epigene_sds)[3] <- 'sd'
  rownames(epigene_sds) <- NULL
  cancer_df <- rbind(cancer_df, epigene_sds)
}

write.csv(cancer_df, '../../Summary_tables/TARGET_NMF_input_variable_epifactors.csv', row.names=F)


