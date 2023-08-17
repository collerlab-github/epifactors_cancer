library(dplyr)
library(ggplot2)


# cancers <- c('ACC','LGG','LUAD')
cancers <- c('LGG')
cancer_df <- data.frame(matrix(nrow = length(cancers), ncol=9), row.names=cancers)

for (cancer in cancers){
  val_df <- read.csv(paste0('../../validation_',cancer,'/microarray_data/',cancer,'_epigenes_log2_quantile_norm_counts.csv'),row.names=1)
  print(cancer)
  print(ncol(val_df))
  epigene_file <- list.files(path=paste0('../../TCGA_',cancer,'/RNA-seq_datasets/'),pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'), full.names = T)
  analysis_epigenes <- ncol(read.csv(epigene_file))
  print(analysis_epigenes)
  cancer_surv <- read.csv(file=paste0('../../validation_',cancer,'/microarray_data/patient_clinical_info.csv'),row.names=1,header=T)
  missing_surv_patients <- rownames(cancer_surv[is.na(cancer_surv[,4]),])
  print(paste0('Patients with survival data: ', nrow(cancer_surv) - length(missing_surv_patients)))
  
}

for (cancer in cancers){
  # validation set stats
  val_df <- read.csv(paste0('../../validation_',cancer,'/microarray_data/',cancer,'_epigenes_log2_quantile_norm_counts.csv'),row.names=1)
  num_epigene_probes <- nrow(val_df) # number of epigenes in microarray dataset
  # tcga validation comparison
  epi_PFI <- read.csv(paste0('../../TCGA_',cancer,'/05_gene_clinical_prediction_analysis/',cancer,'_significant_cox_cluster_differential_genes_PFI.csv'), row.names=1)
  epi_OS <- read.csv(paste0('../../TCGA_',cancer,'/05a_TCGA_clinical_prediction_analysis/',cancer,'_significant_cox_cluster_differential_genes_OS.csv'), row.names=1)
  epigene_file <- list.files(path=paste0('../../TCGA_',cancer,'/RNA-seq_datasets/'),pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'), full.names = T)
  analysis_epigenes <- nrow(read.csv(epigene_file))
  
  rand <- read.csv(list.files(paste0('../../validation_',cancer,'/prognostic/'), pattern=paste0(cancer,'_significant_cox_cluster_differential_genes_'), full.names=T), row.names=1)
  sig_epi_PFI <- rownames(epi_PFI[which(epi_PFI[,2]<=0.05),2,F])
  sig_epi_OS <- rownames(epi_OS[which(epi_OS[,2]<=0.05),2,F])
  sig_rand <- rownames(rand[which(rand[,2]<=0.05),2,F])
  
  cancer_df[cancer, ] <- c(cancer, num_epigene_probes, analysis_epigenes, nrow(rand), length(sig_epi_PFI)/analysis_epigenes, 
                           length(sig_epi_OS)/analysis_epigenes, length(sig_rand)/nrow(rand),
                           length(intersect(sig_epi_PFI,sig_rand)), length(intersect(sig_epi_OS,sig_rand)))
  # tcga_epi <- read.csv(paste0('../../TCGA_',cancer,'/05_gene_clinical_prediction_analysis/',cancer,'_significant_cox_cluster_differential_genes_PFI.csv'), row.names=1)
  # val_epi <- read.csv('../../validation_ACC/prognostic/ACC_significant_cox_cluster_differential_genes_OS.csv',row.names=1)
  # sig_df <- tcga_epi[which(tcga_epi[,2]<=0.05),2,drop=F]
  # sig_df <- merge(sig_df,val_epi,by=0)
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
  if(cancer=='LGG'){
    end <- 'PFI'
    rand[,'prognostic_in_primary_dataset'] <- ifelse(rownames(rand)%in%intersect(sig_epi_PFI,rownames(rand)),'yes','no')
    write.csv(rand,paste0('../../validation_',cancer,'/prognostic/',cancer,'_significant_cox_cluster_differential_genes_',end,'_revised.csv'), quote=F)
    
  }
}

colnames(cancer_df) <- c('cancer','num epigenes in validation', 'num epigenes used for TCGA prognostic analysis', 'num epigenes in TCGA analysis and validation',
                         'TCGA fraction prognostic PFI', 'TCGA fraction prognostic OS','validation fraction prognostic','TCGA PFI validation overlap','TCGA OS validation overlap')

dir.create('../../Summary_tables/validation', showWarnings = T)
write.csv(cancer_df, '../../Summary_tables/validation/TCGA_validation_prognostic_comparison.csv',quote=F)



