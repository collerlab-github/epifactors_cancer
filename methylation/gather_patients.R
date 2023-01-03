library(dplyr)
library(TCGAbiolinks)
library(ggplot2)
library(plyr)
library(reshape2)
library(stringr)
# current TCGA patients in our study
clinical_data <- read.csv('../../Summary_tables/nmf_patient_info_all.csv',row.names=1,header=T)
clinical_data$project_id[which(is.na(clinical_data$project_id))] <- 'NA'
sample_id <- rownames(clinical_data)
# TCGA project id's
cancer_list <- paste0('TCGA-',c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','COAD','READ',
                 'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                 'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD'))

# tcga query for patients 
q <- GDCquery(project=cancer_list,legacy = F,data.category = "DNA Methylation")
# load reults
q_res <- getResults(q) 
# table(sample_id %in% q_res$sample.submitter_id)
# which submitter id's included in the methylation data
q_res$sample.submitter_id <- gsub('-','\\.',q_res$sample.submitter_id)
# patients with metylation data
meth_id <- sample_id[which(sample_id %in% q_res$sample.submitter_id)]
write.csv(q_res[which(q_res$sample.submitter_id %in% meth_id),], 'methylation_summary_all.csv',row.names=F)
# find array platform used for each patient
platforms <- q_res[which(q_res$sample.submitter_id %in% meth_id),c('sample.submitter_id','platform')]
platforms$project_id <- clinical_data[platforms$sample.submitter_id,'project_id']
# platforms$project_id[which(is.na(platforms$project_id))] <- 'NA'
sort_paste <- function(x,iter=1) {
  if (iter==2){
    x <- str_split(x,';')[[1]]
  }
  return(paste(unique(sort(x)),collapse=';'))
}
# array platforms for each cancer
cancer_platform <- aggregate(platform~project_id, data=platforms, FUN = sort_paste)
# subset for patients with methylation data
num_patients_df <- data.frame(count(clinical_data,'project_id'),count(platforms,'project_id'))[,-3]
num_patients_df$project_id[which(is.na(num_patients_df$project_id))] <- 'NA'
colnames(num_patients_df) <- c('project_id','all','methylation')
num_patients_df <- merge(num_patients_df, cancer_platform, by.x=1, by.y=1)

write.csv(num_patients_df,'../../metadata/number_of_methylation_patients.csv',row.names=F)
num_patients_df <- melt(num_patients_df,id.vars='project_id')
colnames(num_patients_df)[2:3] <- c('Group','Frequency')
ggplot(num_patients_df,aes(x=project_id,y=Frequency,fill=Group)) +
  geom_bar(stat="identity",position="dodge")+
  theme(axis.text.x = element_text(angle=45,size=10,hjust=1)) +
  ggtitle('Number of Methylation Patients')
  # geom_text(aes(label=Frequency), vjust=-1, hjust=0)
ggsave('number_of_methylation_patients.png',width=7,height=7)
# getManifest(q_res%>%filter(sample.submitter_id%in%meth_id),save = TRUE) 

