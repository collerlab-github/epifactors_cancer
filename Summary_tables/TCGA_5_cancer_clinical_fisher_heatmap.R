library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(reshape2)
library(RColorBrewer)
library(circlize)

cancers <- c('ACC','KIRC','LGG','LIHC','LUAD')
fisher_df <- matrix(nrow=5,ncol=length(cancers))
rownames(fisher_df) <- c('stage','grade','pathologic_t','pathologic_n','pathologic_m')
colnames(fisher_df) <- cancers
for (cancer in cancers) {
  fishers <- read.csv(paste0('../../TCGA_',cancer,'/04_clinical_analysis/Rank_2/',cancer,'_rank_2_fisher_pvalues.csv'),row.names=1)
  fisher_df[,cancer] <- fishers$x
}

fisher_df[is.na(fisher_df)] <- 1
# heatmap
color_border <- brewer.pal(n=9,"Reds")[c(8,1)]
col_fun <- colorRamp2(c(0,0.049999999,0.05),c(color_border,'grey'))
png("../../Summary_tables/TCGA_NMF_5_cancer_clinical_fisher_heatmap.png", units="in", width=15, height=10, res=500)
h <- Heatmap(fisher_df,
             name = 'pvalue', na_col = 'black',#  cluster_rows = F, cluster_columns = F,
             heatmap_legend_param = list(title='p-value', at=c(0:5)/100, labels=c(c(0:4)/100,'>=0.05'),
                                         labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                         title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                         legend_height = unit(2.5, "in")),
             column_title = 'Cancers',
             column_title_gp = gpar(fontsize=20,fontfamily='Helvetica'),
             column_title_side = "bottom",
             row_title = 'Clinical Information',
             row_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
             row_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
             column_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
             column_names_rot = 60,
             col=col_fun)
# add more titles to heatmap
draw(h,
     column_title="NMF Cluster Clinical Stratification Fisher's P-value", column_title_gp=gpar(fontsize=25, fontfamily='Helvetica')
)
dev.off()

