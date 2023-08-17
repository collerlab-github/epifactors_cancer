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
  cancer_list <- c('LGG','PRAD','UCEC','KIRC','CRC',
                   'LIHC','KIRP','LUAD','STAD','ACC')
  pval_df <- data.frame(OS=rep(1,length(cancer_list)),DSS=rep(1,length(cancer_list)),PFI=rep(1,length(cancer_list)),row.names=cancer_list)
} else if (project=='TARGET'){
  cancer_list <- c('AML','NBL','OS','WT')
  pval_df <- data.frame(OS=rep(0,length(cancer_list)),row.names=cancer_list)
} else{
  quit(save='no')
}

rank <- 2

################

# consolidate nmf survival curves p-values into one dataframe
cluster_df <- pval_df
age_df <- pval_df
gender_df <- pval_df
stage_df <- pval_df
grade_df <- pval_df
# read in each cancer's survival curve p-values
for (c in cancer_list) {
  # file
  dir <- paste0('../',project,'_',c,'/04_clinical_analysis/Rank_',rank,'/')
  print(c)
  dir_files <- list.files(dir,pattern='endpoint_pval_')
  cancer_file <- paste0(dir,dir_files[grepl('stage|grade',dir_files)])
  dir_files <- list.files(dir, pattern='cox_pval_')
  cox_file <- paste0(dir,dir_files[grepl('stage|grade',dir_files)])
  print(cancer_file)
  print(cox_file)
  # read and fill in cancer survival p-value
  if (!file_test("-f",cancer_file) | !file_test('-f',cox_file)){ # skip if file doesn't exist
    next
  }
  surv_pvals <- read.csv(cancer_file, header=T, row.names=1)
  cox_pvals <- read.csv(cox_file, header=T, row.names=1)
  pval_df[c,] <- surv_pvals[1,]
  cluster_df[c,] <- cox_pvals['cluster',]
  age_df[c,] <- cox_pvals['age',]
  gender_df[c,] <- cox_pvals['gender',]
  if ('stage' %in% rownames(cox_pvals)){
    stage_df[c,] <- cox_pvals['stage',]
  }
  if ('grade' %in% rownames(cox_pvals)){
    grade_df[c,] <- cox_pvals['grade',]
  }
  
}

# write p-value table
write.csv(pval_df, file = paste0(project,'_nmf_survival_pvalues_sg.csv'),row.names=T, quote=F)
write.csv(cluster_df, file = paste0(project,'_nmf_cox_cluster_pvalues_sg.csv'),row.names=T, quote=F)
write.csv(age_df, file = paste0(project,'_nmf_cox_age_pvalues_sg.csv'),row.names=T, quote=F)
write.csv(gender_df, file = paste0(project,'_nmf_cox_gender_pvalues_sg.csv'),row.names=T, quote=F)
write.csv(stage_df, file = paste0(project,'_nmf_cox_stage_pvalues_sg.csv'),row.names=T, quote=F)
write.csv(grade_df, file = paste0(project,'_nmf_cox_grade_pvalues_sg.csv'),row.names=T, quote=F)

# ####################
# # pval_df heatmaps #
# ####################
# # color-coded pval_df based on significance between surivival curves for each endpoint
# pval_color_df <- ifelse(pval_df<0.05,'pval<0.05','pval>=0.05')
# # function to label heatmap cells
# cell_function <- function(j,i,x,y,width,height,fill){
#   grid.text(label = round(pval_df[i,j],digits=7),x = x,y = y)
# }
# # no text heatmap
# pval_sig <- pval_df
# pval_sig[pval_sig>=0.05] <- NA
# Heatmap(pval_sig,
#         name='p-value',
#         column_title = 'Survival Curve P-values',
#         na_col = "black",
#         )
# # no cluster heatmap
# png(paste0(project,'_nmf_survival_pvalues_no_cluster_sg.png'),units='in',width=7,height=7,res=500)
# Heatmap(pval_color_df,
#         name = 'p-value',
#         column_title = 'Survival Curve P-values',
#         col=structure(c('palegreen','firebrick1'),names=c('pval<0.05','pval>=0.05')), 
#         rect_gp=gpar(col='black',lwd=2),
#         cell_fun = cell_function)
# dev.off()
# 
# # cluster heatmap
# color_function <- colorRamp2(c(0.049999999,0.0500001),c('palegreen','firebrick1'))
# png(paste0(project,'_nmf_survival_pvalues_cluster_sg.png'),units='in',width=7,height=7,res=500)
# h <- draw(Heatmap(as.matrix(pval_df),
#         # name = 'p-value',
#         heatmap_legend_param = list(title='p-value',at=c(0.049999999,0.0500001),labels=c('pval<0.05','pval>=0.05')),
#         column_title = 'Survival Curve P-values',
#         col=color_function,
#         rect_gp=gpar(col='black',lwd=2),
#         cell_fun = cell_function
#         )
# )
# dev.off()
# # cluster heatmap no text
# color_border <- brewer.pal(n=9,"Reds")[c(8,1)]
# color_function <- colorRamp2(c(0,0.049999999,0.05),c(color_border,'grey'))
# png(paste0(project,'_nmf_survival_pvalues_cluster_no_text_sg.png'),units='in',width=7,height=7,res=500)
# h <- draw(Heatmap(as.matrix(pval_df),
#                   # name = 'p-value',
#                   heatmap_legend_param = list(title='p-value',at=c(0:5)/100,labels=c(c(0:4)/100,'>=0.05')),
#                   column_title = 'Survival Curve P-values',
#                   col=color_function,
#                   rect_gp=gpar(col='black',lwd=2),
#                   )
# )
# write.csv(pval_df[row_order(h),column_order(h)], file = paste0(project,'_nmf_survival_pvalues_sg.csv'),row.names=T, quote=F)
# 
# dev.off()


#######################
# cluster_df heatmaps #
#######################
# color-coded pval_df based on significance between surivival curves for each endpoint
cluster_color_df <- ifelse(cluster_df<0.05,'pval<0.05','pval>=0.05')
# function to label heatmap cells
cell_function <- function(j,i,x,y,width,height,fill){
  grid.text(label = round(cluster_df[i,j],digits=7),x = x,y = y)
}
# no cluster heatmap
png(paste0(project,'_nmf_cox_cluster_pvalues_no_cluster_sg.png'),units='in',width=7,height=7,res=500)
h <- draw(Heatmap(cluster_color_df,
        name = 'p-value',
        column_title = 'Cox Regression Cluster P-values',
        col=structure(c('palegreen','firebrick1'),names=c('pval<0.05','pval>=0.05')), 
        rect_gp=gpar(col='black',lwd=2),
        cell_fun = cell_function))
dev.off()

# cluster heatmap (only for TCGA) with and without text
if (project=='TCGA'){
  color_function <- colorRamp2(c(0.049999999,0.0500001),c('palegreen','firebrick1'))
  # order cluster_df according to clustering from pval_df heatmap
  cluster_df <- cluster_df[c('ACC','KIRC','LGG','LIHC','LUAD','CRC','KIRP','PRAD','STAD','UCEC'),]
  cluster_color_df <- ifelse(cluster_df<0.05,'pval<0.05','pval>=0.05')
  png(paste0(project,'_nmf_cox_cluster_pvalues_cluster_sg.png'),units='in',width=7,height=7,res=500)
  Heatmap(cluster_color_df,
          name = 'p-value',
          column_title = 'Cox Regression Cluster P-values',
          cluster_rows = T,
          col=structure(c('palegreen','firebrick1'),names=c('pval<0.05','pval>=0.05')), 
          rect_gp=gpar(col='black',lwd=2),
          cell_fun = cell_function)
  dev.off()
  # cluster heatmap no text
  color_border <- brewer.pal(n=9,"Reds")[c(8,1)]
  color_function <- colorRamp2(c(0,0.049999999,0.05),c(color_border,'grey'))
  png(paste0(project,'_nmf_cox_cluster_pvalues_cluster_no_text_sg.png'),units='in',width=7,height=7,res=500)
  h <- draw(Heatmap(as.matrix(cluster_df), cluster_rows=F, cluster_columns = F,
                    # name = 'p-value',
                    heatmap_legend_param = list(title='p-value',at=c(0:5)/100,labels=c(c(0:4)/100,'>=0.05')),
                    column_title = 'Cox Regression Cluster P-values',
                    col=color_function,
                    rect_gp=gpar(col='black',lwd=2),
                    )
  )
  dev.off()
  write.csv(cluster_df[row_order(h),column_order(h)], file = paste0(project,'_nmf_cox_cluster_pvalues_sg.csv'),row.names=T, quote=F)
}


#######################
# age_df heatmaps #
#######################
# color-coded pval_df based on significance between surivival curves for each endpoint
age_color_df <- ifelse(age_df<0.05,'pval<0.05','pval>=0.05')
# function to label heatmap cells
cell_function <- function(j,i,x,y,width,height,fill){
  grid.text(label = round(age_df[i,j],digits=8),x = x,y = y)
}
# no cluster heatmap
png(paste0(project,'_nmf_cox_age_pvalues_no_cluster_sg.png'),units='in',width=7,height=7,res=500)
Heatmap(age_color_df,
        name = 'p-value',
        column_title = 'Cox Regression Age P-values',
        col=structure(c('palegreen','firebrick1'),names=c('pval<0.05','pval>=0.05')), 
        rect_gp=gpar(col='black',lwd=2),
        cell_fun = cell_function)
dev.off()

# cluster heatmap (only for TCGA)
if(project=='TCGA'){
  color_function <- colorRamp2(c(0.049999999,0.0500001),c('palegreen','firebrick1'))
  # order cluster_df according to clustering from pval_df heatmap
  age_df <- age_df[c('PAAD','ESCA','OV','LUSC','GBM','CESC','STAD','BRCA','PRAD','SARC','PCPG',
                             'BLCA','TGCT','CRC','SKCM','HNSC','THCA','KIRP','LIHC','LUAD','UCEC','LGG','KIRC','ACC'),]
  age_color_df <- ifelse(age_df<0.05,'pval<0.05','pval>=0.05')
  png('nmf_cox_age_pvalues_cluster_sg.png',units='in',width=7,height=7,res=500)
  Heatmap(age_color_df,
          name = 'p-value',
          column_title = 'Cox Regression Age P-values',
          col=structure(c('palegreen','firebrick1'),names=c('pval<0.05','pval>=0.05')), 
          rect_gp=gpar(col='black',lwd=2),
          cell_fun = cell_function)
  dev.off()
  # cluster heatmap no text
  color_border <- brewer.pal(n=9,"Reds")[c(8,1)]
  color_function <- colorRamp2(c(0,0.049999999,0.05),c(color_border,'grey'))
  png(paste0(project,'_nmf_cox_age_pvalues_cluster_no_text_sg.png'),units='in',width=7,height=7,res=500)
  h <- draw(Heatmap(as.matrix(age_df), cluster_rows=F, cluster_columns = F,
                    # name = 'p-value',
                    heatmap_legend_param = list(title='p-value',at=c(0:5)/100,labels=c(c(0:4)/100,'>=0.05')),
                    column_title = 'Cox Regression Age P-values',
                    col=color_function,
                    rect_gp=gpar(col='black',lwd=2),
  )
  )
  dev.off()
}

#######################
# gender_df heatmaps #
#######################
# color-coded gender_df based on significance between surivival curves for each endpoint
gender_color_df <- ifelse(gender_df<0.05,'pval<0.05','pval>=0.05')
# function to label heatmap cells
cell_function <- function(j,i,x,y,width,height,fill){
  grid.text(label = round(gender_df[i,j],digits=8),x = x,y = y)
}
# no cluster heatmap
png(paste0(project,'_nmf_cox_gender_pvalues_no_cluster.png'),units='in',width=7,height=7,res=500)
Heatmap(gender_color_df,
        name = 'p-value',
        column_title = 'Cox Regression Gender P-values',
        col=structure(c('palegreen','firebrick1'),names=c('pval<0.05','pval>=0.05')), 
        rect_gp=gpar(col='black',lwd=2),
        cell_fun = cell_function)
dev.off()

# cluster heatmap (only for TCGA)
if(project=='TCGA'){
  color_function <- colorRamp2(c(0.049999999,0.0500001),c('palegreen','firebrick1'))
  # order cluster_df according to clustering from pval_df heatmap
  gender_df <- gender_df[c('PAAD','ESCA','OV','LUSC','GBM','CESC','STAD','BRCA','PRAD','SARC','PCPG',
                     'BLCA','TGCT','CRC','SKCM','HNSC','THCA','KIRP','LIHC','LUAD','UCEC','LGG','KIRC','ACC'),]
  gender_color_df <- ifelse(gender_df<0.05,'pval<0.05','pval>=0.05')
  png('nmf_cox_gender_pvalues_cluster.png',units='in',width=7,height=7,res=500)
  Heatmap(gender_color_df,
          name = 'p-value',
          column_title = 'Cox Regression Gender P-values',
          col=structure(c('palegreen','firebrick1'),names=c('pval<0.05','pval>=0.05')), 
          rect_gp=gpar(col='black',lwd=2),
          cell_fun = cell_function)
  dev.off()
}

# for TARGET cancers, plot all cox regression pvalues in one heatmap
cox_df <- cbind(cluster_df, age_df, gender_df)
colnames(cox_df) <- c('Cluster','Age','Gender')
write.csv(cox_df,paste0(project,'_nmf_cox_all_pvalues.csv'))
# color-coded gender_df based on significance between surivival curves for each endpoint
cox_color_df <- ifelse(cox_df<0.05,'pval<0.05','pval>=0.05')
# function to label heatmap cells
cell_function <- function(j,i,x,y,width,height,fill){
  grid.text(label = round(cox_df[i,j],digits=8),x = x,y = y)
}

# no cluster heatmap
png(paste0(project,'_nmf_cox_all_pvalues_no_cluster.png'),units='in',width=7,height=7,res=500)
Heatmap(cox_color_df,
        name = 'p-value',
        column_title = 'Cox Regression P-values',
        col=structure(c('palegreen','firebrick1'),names=c('pval<0.05','pval>=0.05')), 
        rect_gp=gpar(col='black',lwd=2),
        cell_fun = cell_function)
dev.off()
