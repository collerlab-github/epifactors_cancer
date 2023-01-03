library(ComplexHeatmap)
library(RColorBrewer)
library(ggsci)


prefix <- 'five_cancer'
if(prefix=='five_cancer'){
  pref_title <- 'Five Cancer'
} else if(prefix=='four_cancer'){
  pref_title <- 'Four Cancer'
} else{
  pref_title <- prefix
}

outdir <- paste0('../../nmf_clinical_analysis_figures/',prefix,'/')
# read batch corrected full gene expression patient matrix
gene_exp <- as.data.frame(t(read.csv('../../TCGA_PAN/RNA-seq_datasets/filtered_pancan.csv', row.names=1)))
colnames(gene_exp) <- gsub('-','\\.',colnames(gene_exp))
# top 20 featuers
features <- read.csv(paste0(outdir,'cox-nnet_varimp.csv'))
top20 <- features[1:20,-1]
# test patients
patients <- read.csv(paste0(outdir,'cox-nnet_test_patients.csv'), header=F)[,1]

# PI group
pi_group <- read.csv(paste0(outdir,prefix,'_cox-nnet_survival_curves_final.csv'),row.names=1)
rownames(pi_group) <- patients
pi_group <- pi_group %>% arrange(desc(Group))

# subset gene expression to test patients
top_df <- gene_exp[top20$feature,rownames(pi_group)]


library(circlize)
# color palette from ggsci 
blue <- pal_npg("nrc", alpha = 0.7)(9)[4]
red <- pal_npg("nrc", alpha = 0.7)(9)[8]
plot_heatmap <- function(top_df, rclust=T, filename){
  # make png
  png(file=filename, units="in", width=15, height=10, res=500)
  # specify color, scaled expression, column names, row_names
  color_fun <- colorRamp2(c(-2,0,2), c(blue,'white',red))
  top_df_scaled <- t(scale(t(as.matrix(top_df))))
  column_cluster_order <- c('Low PI','High PI')
  column_split_names <- factor(pi_group$Group, levels=column_cluster_order)
  # row_split_names <- factor(feature_df$feature, levels=row_cluster_order)
  # heatmap
  h <- Heatmap(as.matrix(top_df_scaled), name='Expression', col=color_fun,
               heatmap_legend_param = list(labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                           title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                           legend_height = unit(2.5, "in")),
               show_row_names=TRUE, row_names_gp = gpar(fontsize=9, fontfamily='Helvetica'),
               show_column_names=FALSE,
               # row_split=row_split_names, row_title_gp=gpar(fontsize=20, fontfamily='Helvetica'), 
               column_split=column_split_names, column_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
               cluster_row_slices=F, cluster_rows=rclust,
               cluster_column_slices=F, cluster_columns=T,
               row_gap=unit(2,'mm'), column_gap=unit(2,'mm'), border = TRUE
  )
  # add more titles to heatmap
  ht <- draw(h,
             row_title='Top ML Genes', row_title_gp=gpar(fontsize=25, fontfamily='Helvetica'),
             column_title=paste0(pref_title, ' Patients'), column_title_gp=gpar(fontsize=25, fontfamily='Helvetica'),
  )
  saveRDS(row_order(ht), gsub('.png','_order.RDS',filename))
  dev.off()
}
# plot top genes heatmap with row_clustering
plot_heatmap(top_df=top_df, filename=paste0(outdir,prefix,"_topgenes_heatmap.png"))

