
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                 'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                 'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
shap_features <- read.csv('../../prognostic_ML_classifier/top_shap_feature_SV.csv', row.names=1)
shap_features <- shap_features[1:20,,drop=F]

for (c in cancer_list){
  
  cox_pval <- read.csv(paste0('../../TCGA_',c,'/05_gene_clinical_prediction_analysis/',c,'_significant_cox_cluster_differential_genes_PFI.csv'), row.names=1)
  shap_features[,c] <- cox_pval[rownames(shap_features),"adjpval_cluster"]
}

shap_features <- shap_features[,-1]
shap_features[is.na(shap_features)] <- 1

# plot heatmap
# color-coded pval_df based on significance between survival curves for each endpoint
# color_function <- colorRamp2(c(0.049999999,0.0500001),c('palegreen','orange'))
color_border <- brewer.pal(n=9,"Greens")[c(8,1)]
color_function <- colorRamp2(c(0,0.049999999,0.05),c(color_border,'orange'))
# function to label heatmap cells
cell_function <- function(j,i,x,y,width,height,fill){
  grid.text(label = round(shap_features[i,j],digits=7),x = x,y = y)
}
png('../../prognostic_ML_classifier/batch_corrected/shap_overlap_single_gene_pfi_pvalues_cluster.png',units='in',width=7,height=7,res=500)
h <- draw(Heatmap(as.matrix(shap_features),
                  # name = 'p-value',
                  heatmap_legend_param = list(title='p-value',at=c(0:5)/100,labels=c(c(0:4)/100,'>=0.05')),
                  column_title = 'Survival Curve P-values',
                  col=color_function,
                  rect_gp=gpar(col='black',lwd=2),
                  # cell_fun = cell_function
)
)
dev.off()
write.csv(shap_features[row_order(h),column_order(h)], '../../prognostic_ML_classifier/batch_corrected/shap_features_cancer_type_pfi_cluster_pval.csv', row.names=T)
