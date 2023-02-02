library(dplyr)
library(stringr)
library(reshape2)
library(readxl)

project <- 'TCGA'
cancer <- 'LUAD'

setwd(paste0('../../',project,'_',cancer,'/06_diff_gene_info/'))
# load data
cnv_df <- read.csv(paste0(cancer,'_low_high_cna_all.csv'))
table(cnv_df$total_num_of_samples_high)
table(cnv_df$total_num_of_samples_low)
cnv_df$total_num_of_samples_high[is.na(cnv_df$total_num_of_samples_high)] <- names(table(cnv_df$total_num_of_samples_high))
cnv_df$profiled_samples_high[is.na(cnv_df$profiled_samples_high)] <- 0
cnv_df$total_num_of_samples_low[is.na(cnv_df$total_num_of_samples_low)] <- names(table(cnv_df$total_num_of_samples_low))
cnv_df$profiled_samples_low[is.na(cnv_df$profiled_samples_low)] <- 0

cnv_df <- cnv_df[which(as.numeric(sub("%", "", cnv_df$freq_low)) >= 10.0 | as.numeric(sub("%", "", cnv_df$freq_high)) >=10.0),]


epigenes <- read.csv("../../metadata/epigenes.csv", sep='\t')[,1]
table(cnv_df$gene_symbol %in% epigenes)
sig_top_nmf <- data.frame(read_excel('../../../manuscript/Supplementary_Data_4.xlsx',sheet=3, skip = 1) %>% filter(cancer_type==cancer))
rownames(sig_top_nmf) <- sig_top_nmf$top_nmf_gene_symbol
table(cnv_df$gene_symbol %in% rownames(sig_top_nmf))

#  fisher's exact test for each gene
pvals <- c()
is_epigene <- c()
is_sig_top_nmf <- c()
for (i in 1:nrow(cnv_df)){
  gene <- cnv_df[i,"gene_symbol"]
  # contingency table of 2 categorical variables for fisher's exact test
  low_mut <- cnv_df[i,'profiled_samples_low']
  low_none <- as.numeric(cnv_df[i,'total_num_of_samples_low']) - low_mut
  high_mut <- cnv_df[i,'profiled_samples_high']
  high_none <- as.numeric(cnv_df[i,'total_num_of_samples_high'])-high_mut
  
  cont_table <- data.frame('high'=c(high_mut,high_none), 'low'=c(low_mut,low_none), row.names=c('cnv','no_cnv'))
  fish <- fisher.test(cont_table)
  # pval <- format(signif(fish$p.value,5),scientific = T)
  pval <- fish$p.value
  pvals <- c(pvals,pval)
  is_epigene <- c(is_epigene, ifelse(gene %in% epigenes, 'yes', NA))
  is_sig_top_nmf <- c(is_sig_top_nmf, ifelse(gene %in% rownames(sig_top_nmf), sig_top_nmf[gene,'gene_set'], NA))
}

BH_pval <- p.adjust(pvals, method='BH')

cnv_df[,c('Fisher_pval', 'Fisher_BH_pval', 'is_epigene','is_signature_top_nmf')] <- list(pvals, BH_pval, is_epigene, is_sig_top_nmf)

# output: all OG columns, fisher pval, BH pval, is_epigene (yes or NA), is_signature_top_nmf (low_up or high_up or NA) look at supp table 4
cnv_df <- arrange(cnv_df, Fisher_pval)
table(cnv_df$Fisher_BH_pval <=0.05)
write.csv(cnv_df,paste0(cancer,'_cnv_fisher.csv'))
