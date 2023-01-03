
library(dplyr)
library(ggplot2)
# cancer types

cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                 'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                 'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
rank <- 2

# read in clinical data for all patients
clinical_df <- read.delim('../../metadata/TCGA_clinical_data/clinical.tsv', header=T)

clinical_df <- clinical_df[!duplicated(clinical_df$case_submitter_id),]
clinical_df[1:5,1:5]
clinical_df$case_submitter_id <- gsub('-','\\.',clinical_df$case_submitter_id)
# clinical_subset <- clinical_df %>% select(case_submitter_id, project_id, ethnicity, metastasis_at_diagnosis, tumor_stage, tumor_grade)
nmf_df <- data.frame()
top_genes_df <- data.frame()
for (c in cancer_list) {
  print(c)
  # top genes df
  top_df <- data.frame('Gene'=c(rownames(read.csv(paste0('../../TCGA_',c,'/03_nmf/Rank_2/nmf_lee_rank2_feature1_genes.csv'),row.names=1)),
                                rownames(read.csv(paste0('../../TCGA_',c,'/03_nmf/Rank_2/nmf_lee_rank2_feature2_genes.csv'),row.names=1))),
                       'Cancer_Type'=c)
  top_genes_df <- rbind(top_genes_df,top_df)
  df <- read.csv(paste0('../../TCGA_',c,'/03_nmf/Rank_2/patient_clinical_info.csv'))
  df <- df[-nrow(df),]
  df$cancer_type=paste0('TCGA_',c)
  nmf_df <- rbind(nmf_df,df)
  # correlation between age distributions of each cluster
  clust1 <- (df %>% filter(cluster==1))[,'age',drop=T]
  clust2 <- (df %>% filter(cluster==2))[,'age',drop=T]
  # test normality for each cluster
  shap1 <- shapiro.test(clust1)$p.value
  shap2 <- shapiro.test(clust2)$p.value
  # print(shap1)
  # print(shap2)
  # use t test
  pval = t.test(clust1,clust2)$p.value
  # use mann whitney test
  pval_mw = wilcox.test(clust1,clust2, exact = F)$p.value

  print(pval)
  print(pval_mw)
  ggplot(df, aes(x=as.character(cluster), y=as.numeric(age), group=cluster))+
  geom_boxplot()+
  annotate(geom="text", x=1.5, y=max(df$age,na.rm = T)+5, label=paste0("t-test p=",round(pval,5)), color="Black", size=5)+
  annotate(geom="text", x=1.5, y=max(df$age,na.rm = T)+10, label=paste0("mann-whitney p=",round(pval_mw,5)), color="Black", size=5)+
  labs(title= paste0(c,' NMF Rank 2 Clusters by Age'), x = "Cluster", y="Age")
  ggsave(paste0('../../TCGA_', c, '/03_nmf/Rank_2/age_correlation.png'), width=7, height=7)
}
patients <- nmf_df$X
rownames(nmf_df) <- patients
nmf_df <- nmf_df[,-1]
write.csv(nmf_df,'../../Summary_tables/nmf_patient_info_all_revised.csv')
write.csv(top_genes_df,'../../Summary_tables/TCGA_nmf_top_genes_all.csv',row.names=F)

# inner join nmf_df patients with full clinical_df patients

patients <- sort(patients)
nmf_df$X <- gsub('\\.\\d+A','',rownames(nmf_df))
nmf_df1 <- merge(nmf_df,clinical_df, by.x='X',by.y='case_submitter_id',all.x=T)
rownames(nmf_df1) <- patients
colnames(nmf_df1)[1] <- 'case_submitter_id'
write.csv(nmf_df1,'../../Summary_tables/nmf_patient_info_revised_with_TCGA_clinical.csv', row.names=T)
