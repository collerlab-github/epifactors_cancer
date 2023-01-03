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

end <- 'PFI'


overlap_df <- data.frame(cancer=c(),overlap=c())
for (c in cancer_list){
  nmf_genes <- c(read.csv(paste0('../',project,'_',c,'/03_nmf/Rank_2/nmf_lee_rank2_feature1_genes.csv'))[,1], read.csv(paste0('../',project,'_',c,'/03_nmf/Rank_2/nmf_lee_rank2_feature2_genes.csv'))[,1])
  prognostic_file <- paste0('../',project,'_',c,'/05_gene_clinical_prediction_analysis/',c,'_significant_cox_cluster_differential_genes_',end,'.csv')
  prognostic_genes <- read.csv(prognostic_file, row.names=1)
  prognostic_genes <- rownames(prognostic_genes)[which(prognostic_genes$adjpval_cluster<0.05)]
  overlap <- sum(nmf_genes %in% prognostic_genes)/length(nmf_genes)
  overlap_df <- rbind(overlap_df, c(cancer=c,overlap=overlap))
}
colnames(overlap_df) <- c('cancer','overlap')

write.csv(overlap_df,paste0(project,'_nmf_single_gene_overlap_',end,'.csv'))
