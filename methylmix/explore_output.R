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
library(stringr)
library(ggplot2)


cancer <- args[1]
print(cancer)
project <- 'TCGA'

setwd(paste0('../../',project,'_',cancer,'/methylmix/'))


mapped_input_data <- readRDS('mapped_input_data.rds')
MethylMixResults <- readRDS('MethylMix_Results.rds')
methyl_drivers <- data.frame(gene = str_split(MethylMixResults$MethylationDrivers,'---', simplify=T)[,1], num_components=MethylMixResults$NrComponents,
                             row.names=MethylMixResults$MethylationDrivers)
nmf_genes <- c(read.csv('../nmf_lee_rank2_feature1_genes.csv')[,1], read.csv('../nmf_lee_rank2_feature2_genes.csv')[,1])

methyl_nmf_drivers <- methyl_drivers %>% filter(gene %in% nmf_genes)
write.csv(methyl_drivers, 'methylation_drivers.csv')
write.csv(methyl_nmf_drivers, 'methylation_drivers_nmf.csv')
# nmf cluster groups
cluster_membership <- read.csv('../nmf_lee_rank2_cluster_membership.csv', row.names=1)
event_group <- read.csv(paste0('../../metadata/',project,'_nmf_cluster_pfi_prognostic_groups.csv'),row.names=1,header=T)
if (event_group[cancer,1]=='high_pfi'){
  high_samples <- rownames(cluster_membership[which(cluster_membership$cluster==1),,F])
  low_samples <- rownames(cluster_membership[which(cluster_membership$cluster==2),,F])
} else {
  high_samples <- rownames(cluster_membership[which(cluster_membership$cluster==2),,F])
  low_samples <- rownames(cluster_membership[which(cluster_membership$cluster==1),,F])
}

colnames(mapped_input_data$METcancer) <- high_samples
colnames(mapped_input_data$METnormal) <- low_samples

# plot distribution and correlation of mixture componenets
dir.create('plots',showWarnings = F)
for (i in rownames(methyl_nmf_drivers)){
  gene <- methyl_nmf_drivers[i,'gene']
  plots <- MethylMix_PlotModel(i, MethylMixResults, 
                               mapped_input_data$METcancer, 
                               mapped_input_data$hugo_df[,c(high_samples)], 
                               mapped_input_data$METnormal)
  
  
  # plots$MixtureModelPlot
  ggsave(paste0('plots/',i,'_mixture_model_plot.png'), plot=plots$MixtureModelPlot)
  
  # plots$CorrelationPlot
  ggsave(paste0('plots/',i,'_correlation_plot.png'), plot=plots$CorrelationPlot)
}





