# source("http://bioconductor.org/biocLite.R")
# BiocManager::install('TCGAbiolinks',lib="/Users/michaelcheng/Library/R/4.0/library")

library(MethylMix)
library(doParallel)
library(readxl)
library(dplyr)

cancer <- 'LIHC'
project <- 'TCGA'

dir.create(paste0('../../',project,'_',cancer,'/methylmix'), showWarnings = F)
setwd(paste0('../../',project,'_',cancer,'/methylmix/'))

# nmf cluster groups
cluster_membership <- read.csv('../03_nmf/Rank_2/nmf_lee_rank2_cluster_membership.csv', row.names=1)
event_group <- read.csv(paste0('../../metadata/',project,'_nmf_cluster_pfi_prognostic_groups.csv'),row.names=1,header=T)
if (event_group[cancer,1]=='high_pfi'){
  high_samples <- rownames(cluster_membership[which(cluster_membership$cluster==1),,F])
  low_samples <- rownames(cluster_membership[which(cluster_membership$cluster==2),,F])
} else {
  high_samples <- rownames(cluster_membership[which(cluster_membership$cluster==2),,F])
  low_samples <- rownames(cluster_membership[which(cluster_membership$cluster==1),,F])
}


# methylation data
methyl_df <- read.csv(paste0('../methylation/',project,'_',cancer,'_methylation_beta_values_bc.csv'), header=T, row.names=1)
colnames(methyl_df) <- paste0(colnames(methyl_df),'A')
methyl_high <- methyl_df[,high_samples]
methyl_high <- as.matrix(methyl_high)
methyl_low <- methyl_df[,low_samples]
methyl_low <- as.matrix(methyl_low)
rm(methyl_df)

# gene expression data (log2 normalized)
gene_exp <- read.csv('../RNA-seq_datasets/LIHC_log2normalized_counts.csv', row.names=1)
# read in gencode genename index
gencode <- read_excel("../../metadata/gencode.v22.annotation.bed.xlsx",sheet="gencode.v22.annotation.bed", range = cell_cols(4:5),col_names = F)
gencode <- as.data.frame(gencode)

# fix gene names
gencode[,1]<- gsub(';','',gencode[,1])
gencode[,2] <- gsub(';','',gencode[,2])

colnames(gencode) <- c('ENSG','HUGO')
rownames(gencode) <- gencode$ENSG


# convert ensg to hugo
filter_func <- function(df){
  # add corresponding hugo name column to counts data
  df$HUGO <- gencode[match(rownames(df), rownames(gencode)),2]
  # remove duplicate HUGO ids
  dup_ids <- df$HUGO[duplicated(df$HUGO)]
  df[df$HUGO %in% dup_ids,] %>% select(HUGO)
  df2 <- df %>% filter(!duplicated(df$HUGO))
  return(df2)
}

hugo_df <- filter_func(gene_exp)
hugo_to_ensg <- rownames(hugo_df)
names(hugo_to_ensg) <- hugo_df$HUGO
rownames(hugo_df) <- hugo_df$HUGO
hugo_df <- hugo_df[,-ncol(hugo_df)]
colnames(hugo_df) <- gsub('-','\\.',colnames(hugo_df))
hugo_df <- hugo_df[,c(high_samples, low_samples)]
hugo_df <- as.matrix(hugo_df)
rm(gene_exp)

# ClusterProbes: maps methylation probes to genes, and then for each gene, performs hierarchical clustering on the mapped probes (pearson correlation as distance, complete linkage)
# arguments: treatment group (high PFI cluster), control group (low PFI cluster), 
cl <- makeCluster(10)
registerDoParallel(cl)
ClusterProbes(methyl_high, methyl_low, CorThreshold = 0.4)

# run methylmix
MethylMixResults <- MethylMix(methyl_high, hugo_df[,c(high_samples)], methyl_low)
stopCluster(cl)







