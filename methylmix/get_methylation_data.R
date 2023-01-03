args = commandArgs(trailingOnly = T)
if (length(args) !=1 ){
  stop("Must have one argument. Usage: Rscript get_methylation_data.R [cancer]")
}
library(dplyr)
library(ggplot2)
library(plyr)
library(reshape2)
library(stringr)
library(data.table)


# set inputs
cancer <- args[1]
print(cancer)
project <- 'TCGA'
setwd(paste0('../../',project,'_',cancer))


# # read patients and cluster membership
cluster_membership <- read.csv('nmf_lee_rank2_cluster_membership.csv', row.names=1)
rownames(cluster_membership) <- gsub('\\.','-', rownames(cluster_membership))
sample_id <- substr(rownames(cluster_membership),1,str_length(rownames(cluster_membership))-1)

methyl_file <- '/u/scratch/m/mikechen/epigene_analysis/methylation_data/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv'
# methyl_file <- '/u/scratch/m/mikechen/epigene_analysis/methylation_data/sample.tsv'

methyl_df <- fread(methyl_file, select = c('sample',sample_id))
write.csv(methyl_df, paste0(project,'_',cancer,'_methylation_beta_values_bc.csv'),row.names=F)

if(cancer=='LIHC'){
  print('getting normal samples')
  samples <- read.delim('methylation/gdc_sample_sheet.2022-12-28.tsv')
  normal_samples <- samples %>% filter(Sample.Type=='Solid Tissue Normal')
  sample_id <- substr(normal_samples$Sample.ID,1,str_length(normal_samples$Sample.ID)-1)
  write.csv(sample_id,'normal_samples.txt')
  methyl_df <- fread(methyl_file, select = c('sample',sample_id))
  write.csv(methyl_df, paste0(project,'_',cancer,'_methylation_beta_values_bc_normal.csv'),row.names=F)
}
