setwd('../../Summary_tables/')

library(xlsx)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

################
#####INPUTS#####
################
project <- 'TCGA' # TCGA or TARGET
if (project=='TCGA'){
  cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                   'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                   'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
  # initialize pvalue dataframe for cox corrected cluster/stage/grade pvalue comparison
  OS_df <- data.frame(NMF=rep(1,length(cancer_list)),stage=rep(1,length(cancer_list)),grade=rep(1,length(cancer_list)),row.names=cancer_list)
}
rank <- 2

################

# consolidate nmf survival curves p-values into one dataframe
DSS_df <- OS_df
PFI_df <- OS_df
# read in each cancer's survival curve p-values
for (c in cancer_list) {
  # file
  dir <- paste0('../',project,'_',c,'/04_clinical_analysis/Rank_',rank,'/')
  print(c)
  dir_files <- list.files(dir, pattern='cox_pval')
  NMF_file <- paste0(dir,dir_files[grepl('cox_pval.csv',dir_files)])
  stage_file <- paste0(dir,dir_files[grepl('stage_only.csv',dir_files)])
  grade_file <- paste0(dir,dir_files[grepl('grade_only.csv',dir_files)])
  # read and fill in cancer survival p-value
  if (file_test('-f',NMF_file)){ # add pvalues to dataframes if possible
    NMF_coxpvals <- read.csv(NMF_file, header=T, row.names=1)
    OS_df[c,'NMF'] <- NMF_coxpvals['cluster','OS']
    DSS_df[c,'NMF'] <- NMF_coxpvals['cluster','DSS']
    PFI_df[c,'NMF'] <- NMF_coxpvals['cluster','PFI']
  }
  if (file_test('-f',stage_file)){ # add pvalues to dataframes if possible
    stage_coxpvals <- read.csv(stage_file, header=T, row.names=1)
    OS_df[c,'stage'] <- stage_coxpvals['stage','OS']
    DSS_df[c,'stage'] <- stage_coxpvals['stage','DSS']
    PFI_df[c,'stage'] <- stage_coxpvals['stage','PFI']
  }
  if (file_test('-f',grade_file)){ # add pvalues to dataframes if possible
    grade_coxpvals <- read.csv(grade_file, header=T, row.names=1)
    OS_df[c,'grade'] <- grade_coxpvals['grade','OS']
    DSS_df[c,'grade'] <- grade_coxpvals['grade','DSS']
    PFI_df[c,'grade'] <- grade_coxpvals['grade','PFI']
  }
}

# write p-value table
write.csv(OS_df, file=paste0(project,'_nmf_sg_cox_cluster_pvalues_comprisons_OS.csv'), row.names=T, quote=F)
write.csv(DSS_df, file=paste0(project,'_nmf_sg_cox_cluster_pvalues_comprisons_DSS.csv'), row.names=T, quote=F)
write.csv(PFI_df, file=paste0(project,'_nmf_sg_cox_cluster_pvalues_comprisons_PFI.csv'), row.names=T, quote=F)


#######################
# cluster_df heatmaps #
#######################
# cluster heatmap (only for TCGA) with and without text
for(endpoint in c('OS','DSS','PFI')){
  endpoint_df = get(paste0(endpoint,'_df'))
  # function to label heatmap cells
  cell_function <- function(j,i,x,y,width,height,fill){
    grid.text(label = round(endpoint_df[i,j],digits=7),x = x,y = y)
  }
  color_function <- colorRamp2(c(0.049999999,0.0500001),c('palegreen','firebrick1'))
  # order endpoint df
  endpoint_df <- endpoint_df[c('PAAD','ESCA','OV','LUSC','GBM','CESC','STAD','BRCA','PRAD','SARC','PCPG',
                               'BLCA','TGCT','CRC','SKCM','HNSC','THCA','KIRP','LIHC','LUAD','UCEC','LGG','KIRC','ACC'),]
  endpoint_color_df <- ifelse(endpoint_df<0.05,'pval<0.05','pval>=0.05')
  png(paste0(project,'_nmf_sg_cox_cluster_pvalues_comparisons_',endpoint,'.png'),units='in',width=7,height=7,res=500)
  h <- draw(Heatmap(endpoint_color_df,
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
  png(paste0(project,'_nmf_sg_cox_cluster_pvalues_comparisons_',endpoint,'_no_text.png'),units='in',width=7,height=7,res=500)
  h <- draw(Heatmap(as.matrix(endpoint_df), cluster_rows=F, cluster_columns = F,
                    # name = 'p-value',
                    heatmap_legend_param = list(title='p-value',at=c(0:5)/100,labels=c(c(0:4)/100,'>=0.05')),
                    column_title =  paste0('Cox Regression Cluster P-values ', endpoint),
                    col=color_function,
                    rect_gp=gpar(col='black',lwd=2),
                    )
  )
  dev.off()
  write.csv(endpoint_df[row_order(h),column_order(h)], file = paste0(project,'_nmf_sg_cox_cluster_pvalues_comprisons_',endpoint,'.csv'),row.names=T, quote=F)
}


