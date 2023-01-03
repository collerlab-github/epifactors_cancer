setwd('../../Summary_tables/')

library(xlsx)
library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

################
#####INPUTS#####
################
project <- 'TCGA' # TCGA or TARGET
if (project=='TCGA'){
  cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                   'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                   'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
  end <- 'PFI'
} else if (project=='TARGET'){
  cancer_list <- c('AML','NBL','OS','WT')
  end <- 'OS'
} else{
  quit(save='no')
}

# epigenes
epigene_list <- read.csv('../metadata/epigenes.csv', sep='\t')[,1]

# gather significant cox hazard ratios for all cancers
cluster_HR <- as.data.frame(matrix(NA, ncol=length(epigene_list), nrow=length(cancer_list)))
rownames(cluster_HR) <- cancer_list
colnames(cluster_HR) <- epigene_list
age_HR <- cluster_HR
gender_HR <- cluster_HR

# function to obtain hazard ratios of genes with significant cox covariate
get_sig_genes <- function(pval_file,HR_file,covariate,c){
  endpoint_pvals <- read.csv(pval_file, header=T, row.names=1)
  endpoint_pvals <- endpoint_pvals[which(endpoint_pvals[,2]<0.05),2,drop=F]
  endpoint_pvals <- endpoint_pvals[order(endpoint_pvals[,1]),1,F]
  sig_HR <- read.csv(HR_file,header=T,row.names=1)
  sig_HR <- sig_HR[rownames(endpoint_pvals),covariate,drop=F]
  sig_HR_file <- paste0('../',project,'_',c,'/05_gene_clinical_prediction_analysis/',c,'_significant_cox_hazard_ratios_',covariate,'_differential_genes_',end,'.csv')
  write.csv(cbind(endpoint_pvals,sig_HR),sig_HR_file,row.names=T,quote=F)
  return(sig_HR)
}

for (c in cancer_list) {
  # fill in cox hazard ratio data frames for genes with significant covariates in cox regression
  # files for cox pvalues
  cluster_file <- paste0('../',project,'_',c,'/05_gene_clinical_prediction_analysis/',c,'_significant_cox_cluster_differential_genes_',end,'.csv')
  age_file <- paste0('../',project,'_',c,'/05_gene_clinical_prediction_analysis/',c,'_significant_cox_age_differential_genes_',end,'.csv')
  gender_file <- paste0('../',project,'_',c,'/05_gene_clinical_prediction_analysis/',c,'_significant_cox_gender_differential_genes_',end,'.csv')
  HR_file <- paste0('../',project,'_',c,'/05_gene_clinical_prediction_analysis/',c,'_significant_cox_hazard_ratios_differential_genes_',end,'.csv')
  # read and fill in cancer significant endpoint genes; use adjusted p-value
  sig_cluster_genes <- get_sig_genes(cluster_file,HR_file,'cluster',c) # cox cluster hazard ratios
  cluster_HR[c,rownames(sig_cluster_genes)] <- sig_cluster_genes$cluster
  
  sig_age_genes <- get_sig_genes(age_file,HR_file,'age',c) # cox age hazard ratios
  age_HR[c,rownames(sig_age_genes)] <- sig_age_genes$age
  
  sig_gender_genes <- get_sig_genes(gender_file,HR_file,'gender',c) # cox gender hazard ratios
  gender_HR[c,rownames(sig_gender_genes)] <- sig_gender_genes$gender
}

write.csv(cluster_HR, file = paste0(project,'_single_gene_',end,'_cox_hazard_ratio_cluster.csv'),row.names=T, quote=F)
write.csv(age_HR, file = paste0(project,'_single_gene_',end,'_cox_hazard_ratio_age.csv'),row.names=T, quote=F)
write.csv(gender_HR, file = paste0(project,'_single_gene_',end,'_cox_hazard_ratio_gender.csv'),row.names=T, quote=F)
