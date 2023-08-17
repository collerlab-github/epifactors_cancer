# install.packages("BiocManager")
# BiocManager::install("GSVA")

library(GSVA)
library(readxl)

# read gene expression data
cancer <- 'LGG'
project <- 'validation'
dir.create(paste0('../../',project,'_',cancer,'/03_nmf_survival'), showWarnings = F)
setwd(paste0('../../',project,'_',cancer,'/03_nmf_survival/'))
log2_counts <- read.csv(paste0('../microarray_data/',cancer,'_epigenes_log2_quantile_norm_counts.csv'),
                  row.names=1)


# read top NMF genes from TCGA cancer folder
event_group <- read.csv('../../metadata/TCGA_nmf_cluster_pfi_prognostic_groups.csv', row.names=1, header=T)
if (event_group[cancer,1]=='high_pfi'){
  nmf_high_genes <- read.csv(paste0('../../TCGA_',cancer,'/03_nmf/Rank_2/nmf_lee_rank2_feature1_genes.csv'))[,1]
  nmf_low_genes <- read.csv(paste0('../../TCGA_',cancer,'/03_nmf/Rank_2/nmf_lee_rank2_feature2_genes.csv'))[,1]
} else {
  nmf_high_genes <- read.csv(paste0('../../TCGA_',cancer,'/03_nmf/Rank_2/nmf_lee_rank2_feature2_genes.csv'))[,1]
  nmf_low_genes <- read.csv(paste0('../../TCGA_',cancer,'/03_nmf/Rank_2/nmf_lee_rank2_feature1_genes.csv'))[,1]
}

gs <- list('NMF_high'=nmf_high_genes, 'NMF_low'=nmf_low_genes)

gsva.es <- gsva(as.matrix(log2_counts), gs, verbose=FALSE)
gsva.es <- as.data.frame(t(gsva.es))
gsva.es$NMF_group <- colnames(gsva.es)[apply(gsva.es, MARGIN=1, FUN=which.max)]

write.csv(gsva.es, paste0(cancer,'_NMFgenes_GSVA_scores.csv'))
