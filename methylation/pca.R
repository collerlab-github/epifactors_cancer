
library(PCAtools)
library(cowplot)
library(reshape2)
cancer <- 'LUAD'
setwd(paste0('../../TCGA_',cancer))
betas <- read.csv(paste0('methylation/TCGA_',cancer,'_methylation_beta_values_bc.csv'), row.names=1)
sample_cluster <- read.csv(paste0("03_nmf/Rank_2/nmf_lee_rank2_cluster_membership.csv"), row.names=1, header=T)
rownames(sample_cluster) <- gsub('01A','01',rownames(sample_cluster))
sample_cluster <- sample_cluster[colnames(betas),,F]
colnames(sample_cluster)[1] <- 'Cluster'
# reorder samples based on high/low NMF group
event_group <- read.csv('../metadata/TCGA_nmf_cluster_pfi_prognostic_groups.csv', row.names=1, header=T)
if (event_group[cancer,1]=='high_pfi'){
  sample_cluster[,1] <- factor(sample_cluster[,1], levels=c(2,1))
} else {
  sample_cluster[,1] <- factor(sample_cluster[,1], levels=c(1,2))
  
}

if (file.exists(paste0('methylation/',tolower(cancer),'_clinical_subtypes_xena_060822.csv'))){
  methyl_cluster <- read.csv(paste0('methylation/',tolower(cancer),'_clinical_subtypes_xena_060822.csv'))
  rownames(methyl_cluster) <- gsub('-','\\.',paste0(methyl_cluster$TCGA_id,'-01'))
  sample_cluster$Methyl <- methyl_cluster[rownames(sample_cluster), "Subtype_DNAmeth"]
  sample_cluster$Methyl[is.na(sample_cluster$Methyl)] <- 'na'
} else{
  sample_cluster$Methyl <- 'na'
}

colnames(sample_cluster)[1] <- 'NMF'

p <- pca(betas,metadata=sample_cluster, removeVar=0.1)

# bi-plot
# with just first two pc's
print('Generating PCA visualizations')
p1 = biplot(p, x='PC1', y='PC2', colby='NMF', shape = 'Methyl', lab=NULL, 
       legendPosition = 'none', legendLabSize=18, legendTitleSize = 15, encircle=T, axisLabSize=20)
# ggsave(paste0('methylation/TCGA_',cancer,'_methylation_pca.png'))

p2 = biplot(p, x='PC1', y='PC3', colby='NMF', shape = 'Methyl', lab=NULL, 
       legendPosition = 'none', legendLabSize=18, legendTitleSize = 15, encircle=T, axisLabSize=20)
# ggsave(paste0('methylation/TCGA_',cancer,'_methylation_pca.png'))

p3 = biplot(p, x='PC2', y='PC3', colby='NMF', shape = 'Methyl', lab=NULL, 
                 legendPosition = 'top', legendLabSize=18, legendTitleSize = 15, encircle=T, axisLabSize=20)
if (cancer == 'KIRC'){
  p3 = p3 + guides(shape="none")
}
legend = get_legend(p3)
title <- ggdraw() + draw_label(paste0('TCGA ', cancer,' Methylation PCA'), fontface='bold', x=0, hjust=0, size=24)+ 
  theme(plot.margin = margin(0,0,0,7), plot.background=element_rect(fill = "white",color='white'))
p3 <- p3 + theme(legend.position='none')
# ggsave(paste0('methylation/TCGA_',cancer,'_methylation_pca.png'))

p4 <- plot_grid(plot_grid(title),plot_grid(p1,p2,p3, align = 'h', ncol=3), plot_grid(legend), rel_heights=c(.1,1,.1), ncol=1)
p4
ggsave(paste0('methylation/TCGA_',cancer,'_methylation_pca1.3.png'), p4,width = 13, height = 10)

p1 = biplot(p, x='PC1', y='PC2', colby='NMF', shape = 'Methyl', lab=NULL, 
            legendPosition = 'right', legendLabSize=18, legendTitleSize = 15, encircle=T, axisLabSize=20, title = paste0(cancer, ' Methylation PCA'), titleLabSize = 20)
if (cancer == 'KIRC'){
  p1 = p1 + guides(shape="none")
}
p1
ggsave(paste0('methylation/TCGA_',cancer,'_methylation_pca.png'), p1)
