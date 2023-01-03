# source("http://bioconductor.org/biocLite.R")
# BiocManager::install('TCGAbiolinks',lib="/Users/michaelcheng/Library/R/4.0/library")
# HOFFMAN2
args = commandArgs(trailingOnly = T)
if (length(args) !=1 ){
  stop("Must have one argument. Usage: Rscript get_methylation_data.R [cancer]")
}

library(MethylMix)
library(doParallel)
library(readxl)
library(dplyr)


cancer <- args[1]
print(cancer)
project <- 'TCGA'

dir.create(paste0('../../',project,'_',cancer,'/methylmix'), showWarnings = F)
setwd(paste0('../../',project,'_',cancer,'/methylmix/'))

# nmf cluster groups
cluster_membership <- read.csv('../nmf_lee_rank2_cluster_membership.csv', row.names=1)
# event_group <- read.csv(paste0('../../metadata/',project,'_nmf_cluster_pfi_prognostic_groups.csv'),row.names=1,header=T)
# if (event_group[cancer,1]=='high_pfi'){
#   high_samples <- rownames(cluster_membership[which(cluster_membership$cluster==1),,F])
#   low_samples <- rownames(cluster_membership[which(cluster_membership$cluster==2),,F])
# } else {
#   high_samples <- rownames(cluster_membership[which(cluster_membership$cluster==2),,F])
#   low_samples <- rownames(cluster_membership[which(cluster_membership$cluster==1),,F])
# }

# methylation data for tumor patients
cl <- makeCluster(4)
registerDoParallel(cl)
methyl_df_tumor <- read.csv(paste0('../methylation/',project,'_',cancer,'_methylation_beta_values_bc.csv'), header=T, row.names=1)
colnames(methyl_df_tumor) <- paste0(colnames(methyl_df_tumor),'A')
methyl_df_tumor <- na.omit(methyl_df_tumor)
methyl_df_tumor <- methyl_df_tumor[,which(rownames(cluster_membership)%in%colnames(methyl_df_tumor))]
methyl_df_tumor <- as.matrix(methyl_df_tumor)


methyl_df_normal <- read.csv(paste0('../methylation/',project,'_',cancer,'_methylation_beta_values_bc_normal.csv'), header=T, row.names=1)
colnames(methyl_df_normal) <- paste0(colnames(methyl_df_normal),'A')
methyl_df_normal <- na.omit(methyl_df_normal)
methyl_df_normal <- as.matrix(methyl_df_normal)

# gene expression data (log2 normalized)
gene_exp <- read.csv(paste0('../RNA-seq_datasets/',cancer,'_log2normalized_counts.csv'), row.names=1)
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
hugo_df <- hugo_df[,colnames(methyl_df_tumor)]
hugo_df <- as.matrix(hugo_df)
rm(gene_exp)

# ClusterProbes: maps methylation probes to genes (from illumina annotation), and then for each gene, performs hierarchical clustering on the mapped probes (pearson correlation as distance, complete linkage)
# arguments: treatment group methylation data (high PFI cluster), control group methylation data (low PFI cluster), 


# methyl_tumor_test <- methyl_df_tumor[1:1000,]
# methyl_normal_test <- methyl_df_normal[1:1000,]
res <-  ClusterProbes(methyl_df_tumor, methyl_df_normal, CorThreshold = 0.4)

# Putting everything together
toSave <- list(METcancer = res[[1]], METnormal = res[[2]], hugo_df = hugo_df, ProbeMapping = res$ProbeMapping)
saveRDS(toSave, file = paste0("mapped_input_data.rds"))

# replace name of methylation probes to their corresponding mapped gene names
colnames(res[[1]]) <- colnames(methyl_df_tumor)
colnames(res[[2]]) <- colnames(methyl_df_normal)
# run methylmix
# arguments: treatment group methylation data, gene expression of treatment group, control group methylation data
MethylMixResults <- MethylMix(res[[1]], hugo_df, res[[2]], OutputRoot = "./")
stopCluster(cl)







