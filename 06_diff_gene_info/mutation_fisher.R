library(dplyr)
library(stringr)
library(reshape2)
library(readxl)

project <- 'TCGA'
cancer <- 'LUAD'

setwd(paste0('../../',project,'_',cancer,'/06_diff_gene_info/'))
# load data
mut_df <- read.csv(paste0(cancer,'_low_high_mut_all.csv'), row.names=1)
table(mut_df$total_num_high)
table(mut_df$total_num_low)
mut_df$total_num_high[is.na(mut_df$total_num_high)] <- names(table(mut_df$total_num_high))
mut_df$num_with_mut_high[is.na(mut_df$num_with_mut_high)] <- 0
mut_df$total_num_low[is.na(mut_df$total_num_low)] <- names(table(mut_df$total_num_low))
mut_df$num_with_mut_low[is.na(mut_df$num_with_mut_low)] <- 0
epigenes <- read.csv("../../metadata/epigenes.csv", sep='\t')[,1]
table(rownames(mut_df) %in% epigenes)
sig_top_nmf <- data.frame(read_excel('../../../manuscript/Supplementary_Data_4.xlsx',sheet=3, skip = 1) %>% filter(cancer_type==cancer))
rownames(sig_top_nmf) <- sig_top_nmf$top_nmf_gene_symbol
table(rownames(mut_df) %in% rownames(sig_top_nmf))

#  fisher's exact test for each gene
pvals <- c()
is_epigene <- c()
is_sig_top_nmf <- c()
for (gene in rownames(mut_df)){
  # contingency table of 2 categorical variables for fisher's exact test
  low_mut <- mut_df[gene,'num_with_mut_low']
  low_none <- as.numeric(mut_df[gene,'total_num_low']) - low_mut
  high_mut <- mut_df[gene,'num_with_mut_high']
  high_none <- as.numeric(mut_df[gene,'total_num_high'])-high_mut
  
  cont_table <- data.frame('high'=c(high_mut,high_none), 'low'=c(low_mut,low_none), row.names=c('mut','no_mut'))
  fish <- fisher.test(cont_table)
  # pval <- format(signif(fish$p.value,5),scientific = T)
  pval <- fish$p.value
  pvals <- c(pvals,pval)
  is_epigene <- c(is_epigene, ifelse(gene %in% epigenes, 'yes', NA))
  is_sig_top_nmf <- c(is_sig_top_nmf, ifelse(gene %in% rownames(sig_top_nmf), sig_top_nmf[gene,'gene_set'], NA))
}

mut_df[,c('Fisher_pval','is_epigene','is_signature_top_nmf')] <- list(pvals, is_epigene, is_sig_top_nmf)

mut_df_final <- mut_df[which(as.numeric(sub("%", "", mut_df$freq_num_with_mut_low)) >= 10.0 | as.numeric(sub("%", "", mut_df$freq_num_with_mut_high)) >=10.0),]

mut_df_final$Fisher_BH_pval <- p.adjust(mut_df_final$Fisher_pval, method='BH')

mut_df_final <- mut_df_final[,c(1:9,12,10,11)]
# output: all OG columns, fisher pval, BH pval, is_epigene (yes or NA), is_signature_top_nmf (low_up or high_up or NA) look at supp table 4
mut_df_final <- arrange(mut_df_final, Fisher_pval)
table(mut_df_final$Fisher_BH_pval <=0.05)
write.csv(mut_df_final,paste0(cancer,'_mutation_fisher.csv'))
