library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

project <- 'TCGA'
# endpoint <- 'PFI'

setwd('../../Summary_tables/')
cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                 'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                 'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')


# create cancer x epigene endpoint pvalue data frame
pval_df <- data.frame()
direction_df <- data.frame()
# rownames(pval_df) <- cancer_list
# colnames(pval_df) <- epigene_list


for (c in cancer_list) {
  # fill in pvalue data frame
  # files for pvalue and cox pvalues
  cancer_file <- paste0('../',project,'_',c,'/EMTome/',c,'_emt_cox_df.csv')
  dir_file <- paste0('../',project,'_',c,'/EMTome/',c,'_emt_direction.csv')
  # read and fill in cancer significant endpoint genes; use adjusted p-value
  if (file.exists(cancer_file)){
    pvals <- read.csv(cancer_file, header=T, row.names=1)['cluster',]
    rownames(pvals) <- c
    pval_df <- rbind(pval_df, pvals)
    direction <- read.csv(dir_file, row.names=1)
    rownames(direction) <- c
    direction_df <- rbind(direction_df, direction)
  }
  else{
    pval_df[c,] <- NA
    direction_df[c,] <- NA
  }
}

write.csv(pval_df, paste0(project,'_EMTome_pvals.csv'), quote=F)
colnames(direction_df) <- 'EMT_direction'
write.csv(direction_df, paste0(project,'_EMTome_direction.csv'), quote=F)

# read grade survival and nmf cluster survival
pval_df[is.na(pval_df)]<- 1
OS_df <- read.csv('TCGA_nmf_sg_cox_cluster_pvalues_comprisons_OS.csv', row.names=1)
OS_df$EMT <- pval_df[rownames(OS_df),'OS']
DSS_df <- read.csv('TCGA_nmf_sg_cox_cluster_pvalues_comprisons_DSS.csv', row.names=1)
DSS_df$EMT <- pval_df[rownames(DSS_df),'DSS']
PFI_df <- read.csv('TCGA_nmf_sg_cox_cluster_pvalues_comprisons_PFI.csv', row.names=1)
PFI_df$EMT <- pval_df[rownames(PFI_df),'PFI']

# function to label heatmap cells

color_function <- colorRamp2(c(0.049999999,0.0500001),c('palegreen','firebrick1'))
# order endpoint df
row_order <- c('PAAD','ESCA','OV','LUSC','GBM','CESC','STAD','BRCA','PRAD','SARC','PCPG',
               'BLCA','TGCT','CRC','SKCM','HNSC','THCA','KIRP','LIHC','LUAD','UCEC','LGG','KIRC','ACC')
# exclude stage
OS_df <- OS_df[row_order,-2]
DSS_df <- DSS_df[row_order,-2]
PFI_df <- PFI_df[row_order,-2]

plot_heatmaps <- function(plot_df, endpoint){
  cell_function <- function(j,i,x,y,width,height,fill){
    grid.text(label = round(plot_df[i,j],digits=7),x = x,y = y)
  }
  pval_color_df <- ifelse(plot_df<0.05,'pval<0.05','pval>=0.05')
  png(paste0(project,'_EMTome_cox_cluster_pvalues_comparisons_',endpoint,'.png'),units='in',width=7,height=7,res=500)
  h <- draw(Heatmap(pval_color_df,
                    name = 'p-value',
                    column_title = paste0('Cox Regression Cluster P-values ', endpoint),
                    cluster_rows = T,
                    col=structure(c('palegreen','firebrick1'),names=c('pval<0.05','pval>=0.05')), 
                    rect_gp=gpar(col='black',lwd=2),
                    cell_fun = cell_function))
  dev.off()
  # cluster heatmap no text
  color_border <- brewer.pal(n=9,"Reds")[c(8,1)]
  color_function <- colorRamp2(c(0,0.049999999,0.05),c(color_border,'grey'))
  png(paste0(project,'_EMTome_cox_cluster_pvalues_comparisons_',endpoint,'_no_text.png'),units='in',width=7,height=7,res=500)
  h <- draw(Heatmap(as.matrix(plot_df), cluster_rows=F, cluster_columns = F,
                    # name = 'p-value',
                    heatmap_legend_param = list(title='p-value',at=c(0:5)/100,labels=c(c(0:4)/100,'>=0.05')),
                    column_title =  paste0('Cox Regression Cluster P-values ', endpoint),
                    col=color_function,
                    rect_gp=gpar(col='black',lwd=2),
                    )
            )
  dev.off()
  write.csv(plot_df[row_order(h),column_order(h)], file = paste0(project,'_EMTome_cox_cluster_pvalues_comprisons_',endpoint,'.csv'),row.names=T, quote=F)
  
}

plot_heatmaps(OS_df, 'OS')
plot_heatmaps(DSS_df, 'DSS')
plot_heatmaps(PFI_df, 'PFI')



