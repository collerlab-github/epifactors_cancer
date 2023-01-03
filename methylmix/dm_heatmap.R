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
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

cancer <- args[1]
print(cancer)
project <- 'TCGA'

setwd(paste0('../../',project,'_',cancer,'/methylmix/'))


# mapped_input_data <- readRDS('mapped_input_data.rds')
# input data
res <- readRDS('MethylMix_Results.rds')

sample_cluster <- read.csv('../nmf_lee_rank2_cluster_membership.csv', row.names=1)
event_group <- read.csv(paste0('../../metadata/',project,'_nmf_cluster_pfi_prognostic_groups.csv'),row.names=1,header=T)
if (event_group[cancer,1]=='high_pfi'){
  sample_cluster[which(sample_cluster$cluster==1),'cluster'] <- 'high'
  sample_cluster[which(sample_cluster$cluster==2),'cluster'] <- 'low'
  feature_order <- c(low=2,high=1)
} else {
  sample_cluster[which(sample_cluster$cluster==1),'cluster'] <- 'low'
  sample_cluster[which(sample_cluster$cluster==2),'cluster'] <- 'high'
  feature_order <- c(low=1,high=2)
}
sample_cluster <- sample_cluster %>% arrange(desc(cluster))
dm_df <- res$MethylationStates[,rownames(sample_cluster)]

# read in top genes for each nmf feature in order of (low, high) (metagene/geneset)
nmf_genes <- c() 
feat <- c() # list of feature identity of each top gene
# looping through each feature index
for (event_group in names(feature_order)){
  feature <- feature_order[event_group]
  # append nmf genes list
  nmf_genes <- c(nmf_genes,rownames(read.csv(paste0('../nmf_lee_rank2_feature',feature,'_genes.csv'),header=T,row.names=1)))
  # append feat list
  feat <- c(feat,rep(as.character(event_group),length(nmf_genes)-length(feat)))
  print(length(nmf_genes))
}
feature_df <- data.frame(feature = paste0(feat),row.names=nmf_genes)

methyl_drivers <- data.frame(gene = str_split(res$MethylationDrivers,'---', simplify=T)[,1], num_components=res$NrComponents,
                             row.names=res$MethylationDrivers)
methyl_nmf_drivers <- methyl_drivers %>% filter(gene %in% nmf_genes)
methyl_nmf_drivers$cluster <- feature_df[methyl_nmf_drivers$gene,feature]
methyl_nmf_drivers <- arrange(methyl_nmf_drivers, desc(cluster))
write.csv(methyl_drivers, 'methylation_drivers.csv')
write.csv(methyl_nmf_drivers, 'methylation_drivers_nmf.csv')


plot_heatmap <- function(top_df, nmf=F, rclust=T, filename){
  # make png
  png(file=filename, units="in", width=15, height=10, res=500)
  # specify color, scaled expression, column names, row_names
  color_fun <- colorRamp2(c(-1,0,1), c('blue','white','red'))
  # top_df_scaled <- t(scale(t(as.matrix(top_df))))
  column_cluster_order <- c('low','high')
  column_split_names <- factor(paste0('Cluster ', sample_cluster$cluster), levels=paste0('Cluster ',column_cluster_order))
  
  if (nmf){
    row_cluster_order <- c('low','high')
    row_split_names <- factor(methyl_nmf_drivers$cluster, levels=row_cluster_order)
    # heatmap
    h <- Heatmap(as.matrix(top_df), name='DM value', col=color_fun,
                 heatmap_legend_param = list(labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                             title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                             legend_height = unit(2.5, "in")),
                 show_row_names=TRUE,
                 show_column_names=FALSE,
                 row_split=row_split_names, row_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
                 column_split=column_split_names, column_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
                 cluster_row_slices=F, cluster_rows=rclust,
                 cluster_column_slices=F, cluster_columns=T,
                 row_gap=unit(2,'mm'), column_gap=unit(2,'mm'), border = TRUE
    )
  } else{
    # heatmap
    h <- Heatmap(as.matrix(top_df), name='DM value', col=color_fun,
                 heatmap_legend_param = list(labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                             title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                             legend_height = unit(2.5, "in")),
                 show_row_names=FALSE,
                 show_column_names=FALSE,
                 # row_split=row_split_names, row_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
                 column_split=column_split_names, column_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
                 cluster_row_slices=F, cluster_rows=rclust,
                 cluster_column_slices=F, cluster_columns=T,
                 row_gap=unit(2,'mm'), column_gap=unit(2,'mm'), border = TRUE
    )
  }
  
  # add more titles to heatmap
  ht <- draw(h,
             row_title='Methylation Driver Genes', row_title_gp=gpar(fontsize=25, fontfamily='Helvetica'),
             column_title=paste0(cancer,' Patients'), column_title_gp=gpar(fontsize=25, fontfamily='Helvetica'),
  )
  saveRDS(row_order(ht), gsub('.png','_order.RDS',filename))
  dev.off()
}

plot_heatmap(top_df=dm_df, filename=paste0('plots/',cancer,"_dm_values_heatmap.png"))
plot_heatmap(top_df=dm_df[rownames(methyl_nmf_drivers),], nmf=T, filename=paste0('plots/',cancer,"_dm_values_heatmap_nmf.png"))


