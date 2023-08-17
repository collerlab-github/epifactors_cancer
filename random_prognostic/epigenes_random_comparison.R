


library(dplyr)
library(ggplot2)
library(reshape2)

setwd('../../Summary_tables/random_prognostic/')
cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                 'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                 'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')

cancer_df <- data.frame(matrix(nrow = length(cancer_list), ncol=7), row.names=cancer_list)
colnames(cancer_df) <- c('cancer','sig_epigenes','analysis_epigenes','fraction_epigenes', 'sig_random','analysis_random','fraction_random')
for (cancer in cancer_list){
  sig_epi <- read.csv(paste0('../../TCGA_',cancer,'/05_gene_clinical_prediction_analysis/',cancer,'_significant_cox_cluster_differential_genes_PFI.csv'), row.names=1)
  epigene_file <- list.files(path=paste0('../../TCGA_',cancer,'/RNA-seq_datasets/'),pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'), full.names = T)
  analysis_epigenes <- nrow(read.csv(epigene_file))
  sig_rand <- read.csv(paste0('../../TCGA_',cancer,'/random_prognostic/',cancer,'_significant_cox_cluster_differential_genes_PFI.csv'), row.names=1)
  n_epi <- nrow(sig_epi[which(sig_epi[,2]<=0.05),2,drop=F])
  n_rand <- nrow(sig_rand[which(sig_rand[,2]<=0.05),2, drop=F])
  cancer_df[cancer, ] <- c(cancer,n_epi, analysis_epigenes, n_epi/analysis_epigenes, n_rand, nrow(sig_rand), n_rand/nrow(sig_rand))
}

write.csv(cancer_df, 'epigenes_v_random_prognostic_comparison.csv', row.names=F)
cancer_df$fraction_epigenes <- as.numeric(cancer_df$fraction_epigenes)
cancer_df$fraction_random <- as.numeric(cancer_df$fraction_random)
ggplot(cancer_df, aes(x=fraction_epigenes, y=fraction_random))+
  geom_point()+
  xlim(0,max(cancer_df$fraction_epigenes))+
  ylim(0,max(cancer_df$fraction_random))+
  ggtitle('Epigenes vs Random PFI Prognostic Genes')
ggsave('epigenes_v_random_prognostic_scatterplot.png')

cancer_df <- melt(data = cancer_df, measure.vars = c('fraction_epigenes','fraction_random'), variable.name = 'gene_set', value.name = 'prognostic_genes')
cancer_df <- cancer_df %>% arrange(gene_set,desc(prognostic_genes))
cancer_df$cancer <- factor(cancer_df$cancer, levels=unique(cancer_df$cancer))
ggplot(cancer_df, aes(x=cancer, y=prognostic_genes, fill=gene_set))+
  geom_bar(position='dodge', stat='identity')+
  ggtitle('Epigenes vs Random PFI Prognostic Genes')+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust = 1))
ggsave('epigenes_v_random_prognostic_barplot.png', width=10)
