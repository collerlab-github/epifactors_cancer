library(dplyr)
library(ggplot2)


# cancers <- c('KIRC','LUAD')
# KIRC: PFI
# LUAD: OS
cancer <- 'KIRC'
end <- 'PFI'
cancer_df <- data.frame(matrix(nrow = 1, ncol=6), row.names=cancer)

# for (cancer in cancers){
#   val_df <- read.csv(paste0('../../validation_',cancer,'/microarray_data/',cancer,'_epigenes_log2_quantile_norm_counts.csv'),row.names=1)
#   print(cancer)
#   print(ncol(val_df))
#   epigene_file <- list.files(path=paste0('../../TCGA_',cancer,'/RNA-seq_datasets/'),pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'), full.names = T)
#   analysis_epigenes <- ncol(read.csv(epigene_file))
#   print(analysis_epigenes)
#   cancer_surv <- read.csv(file=paste0('../../validation_',cancer,'/microarray_data/patient_clinical_info.csv'),row.names=1,header=T)
#   missing_surv_patients <- rownames(cancer_surv[is.na(cancer_surv[,4]),])
#   print(paste0('Patients with survival data: ', nrow(cancer_surv) - length(missing_surv_patients)))
#   
# }


# tcga validation comparison
if(end=='OS'){
  epi_end <- read.csv(paste0('../../TCGA_',cancer,'/05a_TCGA_clinical_prediction_analysis/',cancer,'_significant_cox_cluster_differential_genes_',end,'.csv'), row.names=1)
  
}else{
  epi_end <- read.csv(paste0('../../TCGA_',cancer,'/05_gene_clinical_prediction_analysis/',cancer,'_significant_cox_cluster_differential_genes_',end,'.csv'), row.names=1)
  
}
epigene_file <- list.files(path=paste0('../../TCGA_',cancer,'/RNA-seq_datasets/'),pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'), full.names = T)
analysis_epigenes <- nrow(read.csv(epigene_file))

rand <- read.csv(list.files(paste0('../../TCGA_val_',cancer,'/prognostic/'), pattern=paste0(cancer,paste0('_significant_cox_cluster_differential_genes_',end)), full.names=T), row.names=1)
sig_epi_end <- rownames(epi_end[which(epi_end[,2]<=0.05),2,F])
sig_rand <- rownames(rand[which(rand[,2]<=0.05),2,F])

cancer_df[cancer, ] <- c(cancer, analysis_epigenes, nrow(rand), length(sig_epi_end), 
                         length(sig_rand), length(intersect(sig_epi_end,sig_rand)))


colnames(cancer_df) <- c('cancer', 'num epigenes used for TCGA prognostic analysis', 'num epigenes in validation',
                         paste0('TCGA num prognostic ',end), 'validation num prognostic',paste0('TCGA ',end,' validation overlap'))

write.csv(cancer_df, paste0('../../TCGA_val_',cancer,'/prognostic/TCGA_',cancer,'_validation_prognostic_comparison',end,'.csv'),quote=F)
rand[,'prognostic_in_primary_dataset'] <- ifelse(rownames(rand)%in%intersect(sig_epi_end,rownames(rand)),'yes','no')
write.csv(rand,paste0('../../TCGA_val_',cancer,'/prognostic/',cancer,'_significant_cox_cluster_differential_genes_',end,'_revised.csv'), quote=F)


# sig_df <- epi_end[which(epi_end[,2]<=0.05),2,drop=F]
# sig_df <- merge(sig_df,rand,by=0)
# rownames(sig_df) <- sig_df$Row.names
# sig_df <- sig_df[,-1]
# colnames(sig_df) <- c('TCGA_pval','Val_pval','Val_adjp')
# tcga_gene_order <- rownames(sig_df %>% arrange(TCGA_pval))
# val_gene_order <- rownames(sig_df %>% arrange(Val_pval))
# 
# ggplot(sig_df, aes(x=TCGA_pval, y=Val_pval))+
#   geom_point()+
#   ggtitle(paste0(cancer,' prognostic pvalue comparison'))
# 
# cum_prop_df <- data.frame()
# for (i in 1:nrow(sig_df)){
#   cum_prop_df <- rbind(cum_prop_df,
#                        c(i/nrow(sig_df),length(intersect(tcga_gene_order[1:i],val_gene_order[1:i]))/i)
#   )
# }
# colnames(cum_prop_df) <- c('num_genes','intersection')
# ggplot(cum_prop_df, aes(x=num_genes, y=intersection))+
#   geom_line()
# 
