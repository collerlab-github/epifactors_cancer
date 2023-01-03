library(argparse)

# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("project",nargs=1, help="type of project (TCGA or TARGET)")
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
parser$add_argument("gene", nargs=1, help='gene by which to stratify patients (default is NA)')
# parse the arguments
args <- parser$parse_args()
print(args)
project <- args$project
cancer <- args$cancer
if (project %in% c('TCGA','TARGET') ==F){
  print('project must be TCGA or TARGET')
  quit(save='no')
}
if(project=='TCGA'){
  end <- 'PFI'
} else{
  end <- 'OS'
}
# cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
#                  'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
#                  'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
# cancer <- 'PAAD'
setwd(paste0("../../",project,"_", cancer, "/05_gene_clinical_prediction_analysis"))


library(readxl)
library(survminer)
library(survival)
library(dplyr)

# read in epigene data
print('Reading patient epigene expression data')
epigenes_filename <- list.files(path='../RNA-seq_datasets/',pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'))
if (is.na(args$gene)){
  gene <- read.csv(paste0(cancer,'_significant_pval_differential_genes_',end,'.csv'))[1,1] # most significant PFI-prognostic gene
} else{
  gene <- args$gene # user-specified PFI-prognostic gene
}
sample_cluster <- as.data.frame(t(read.csv(paste0('../RNA-seq_datasets/',epigenes_filename),row.names=1,header=T)))[,gene,F]
patient_info <- read.csv(file=paste0('../03_nmf/Rank_2','/patient_clinical_info.csv'),row.names=1,header=T)



# read in survival data
print('Reading Survival Data')
if (project=='TCGA'){
  # change sample names to have '-' instead of '.'
  rownames(sample_cluster) <- gsub('\\.','-',rownames(sample_cluster))
  rownames(sample_cluster) <- gsub('-01A','',rownames(sample_cluster))
  rownames(patient_info) <- gsub('\\.','-',rownames(patient_info))
  rownames(patient_info) <- gsub('-01A','',rownames(patient_info))
  # sample age and gender
  sample_age <- patient_info[,'age',F]
  sample_gender <- patient_info[,'gender',F]
  # read survival data
  TCGA_surv <- cbind(read_excel('../../metadata/TCGA_survival_data.xlsx',sheet="TCGA-CDR",
                                col_names=T, range= cell_cols(c("B","H"))),
                     read_excel('../../metadata/TCGA_survival_data.xlsx',sheet="TCGA-CDR",
                                col_names=T, range= cell_cols(c("Z","AG"))))
  # check for patients unaccounted for in survival data:
  missing_surv_patients <- rownames(sample_cluster)[which(rownames(sample_cluster) %in% TCGA_surv$bcr_patient_barcode == F)]
  print(paste0('Patients without survival data: ', missing_surv_patients))
  
  # save patient names with added cbioportal tag
  print('Adding cbioportal tag to patient names')
  cbioportal <- read.csv('../../metadata/TCGA_cbioportal_name.csv', row.names=1, header=T)
  study_id <- cbioportal[cancer,'cbioportal.name']
  patients <- paste0(study_id,':',rownames(sample_cluster)[which(rownames(sample_cluster) %in% TCGA_surv$bcr_patient_barcode)])
  
  # extract samples from survival data
  cancer_surv <- TCGA_surv[which(TCGA_surv$bcr_patient_barcode %in% rownames(sample_cluster)),]
  # add age and gender
  cancer_surv$age <- sample_age[cancer_surv$bcr_patient_barcode,]
  cancer_surv$gender <- sample_gender[cancer_surv$bcr_patient_barcode,]
  # add gene_expression data
  cancer_surv <- cbind(cancer_surv, sample_cluster[cancer_surv$bcr_patient_barcode,,F])
  rownames(cancer_surv) <- cancer_surv$bcr_patient_barcode
  # specify endpoint
  end <- 'PFI'
  print(paste0('Endpoint: ', end))
  
} else{ # project is TARGET
  # save patient names with added cbioportal tag
  print('Adding cbioportal tag to patient names')
  cbioportal <- read.csv('../../TARGET_metadata/TARGET_cbioportal_name.csv', row.names=1, header=T)
  study_id <- cbioportal[cancer,'cbioportal.name']
  patients <- paste0(study_id,':',rownames(sample_cluster)[which(rownames(sample_cluster) %in% rownames(patient_info))])
  
  # survival data is in patient info dataframe already
  # extract samples from survival data
  cancer_surv <- patient_info[which(rownames(patient_info)%in%rownames(sample_cluster)),c("OS.time","OS","age","gender")]
  # add gene expression data
  cancer_surv <- cbind(cancer_surv,sample_cluster[rownames(cancer_surv),,F])
  # specify endpoint
  end <- 'OS'
  print(paste0('Endpoint: ', end))
  
}
rm(patients)

# seperate patients into two expression groups for each gene
print('Performing survminer patient categorization by expression')
cutpoint_df <- surv_cutpoint(cancer_surv, time=paste0(end,".time"), event=end, variables = colnames(sample_cluster))
category_df <- surv_categorize(cutpoint_df)
# combine cancer_surv age, gender and patient information with the expression group data
category_df <- cbind(cancer_surv[rownames(category_df),c("age","gender")],category_df)

write.csv(category_df,paste0(cancer,'_top_gene_',gene,'_pfi_groups.csv'))


print('Finding significance between expression groups for each gene')
pval_df <- data.frame(row.names = colnames(category_df)[5:length(colnames(category_df))],pval=rep(0,length(colnames(category_df))-4)) 
# initialize cox regression pvalue df and hazard ratio effect size
cox_df <- data.frame(row.names = colnames(category_df)[5:length(colnames(category_df))],cluster=rep(0,length(colnames(category_df))-4), age=0, gender=0)
cox_HR <- data.frame(row.names = colnames(category_df)[5:length(colnames(category_df))],cluster=rep(0,length(colnames(category_df))-4), CI_low=0, CI_high=0) 
# function that finds pvalue of difference btwn high and low expression categories
# and performs cox regression on age and gender
pval_func <- function(gene, category_df){
  # print(gene)
  # create survival object (survival_time, survival_status)
  s <- Surv(category_df[,paste0(end,'.time')],category_df[,end])
  
  # find p-value difference btwn gene of interest's high/low expression groups
  diff <- survdiff(s~category_df[,gene],data=category_df)
  pval <- pchisq(diff$chisq,df=length(diff$n)-1,lower.tail=F)
  # store p-value in dataframe
  pval_df[gene,'pval'] <<- pval
  
  
  # see if any expression group has no information for age and gender
  row_na <- apply(category_df[,c('age','gender',gene)],1, function(x) sum(is.na(x)))
  low_no_cox <- sum(row_na[which(category_df[,gene]=='low')]) == 2*length(row_na[which(category_df[,gene]=='low')])
  high_no_cox <- sum(row_na[which(category_df[,gene]=='high')]) == 2*length(row_na[which(category_df[,gene]=='high')])
  # print(low_no_cox)
  # print(high_no_cox)
  
  # # compute cox regression and add pvalues to df
  if (low_no_cox | high_no_cox){
    print('one expression group has no covariate data')
    res.cox <- coxph(s ~ category_df[,gene], data =  category_df)
    coxpval <- c(coef(summary(res.cox))[,5],'NA','NA')
    cox_df[gene,] <<- coxpval
    coxHR <- c(coef(summary(res.cox))[1,2],exp(confint(res.cox))[1,])
    cox_HR[gene,] <<- coxHR
  }
  else if(length(unique(category_df$gender))>1){  # run cox regression with all covariates
    # res.cox <- tryCatch({ # try to run cox on all covariates, if it doesn't work, print error and return null
    #   coxph(s ~ category_df[,gene] + age + gender, data =  category_df)
    #   },
    #   error=function(cond) {
    #     message('covariates age and gender do not seem to work together')
    #     message('Here is the original message:')
    #     message(cond)
    #     return(NULL)
    #   })
    res.cox <- coxph(s ~ category_df[,gene] + age + gender, data =  category_df)
    coxpval <- coef(summary(res.cox))[,5]
    cox_df[gene,] <<- coxpval
    coxHR <- c(coef(summary(res.cox))[1,2],exp(confint(res.cox))[1,])
    cox_HR[gene,] <<- coxHR
  }
  else{ # run cox regression on just cluster and age
    res.cox <- coxph(s ~ category_df[,gene] + age, data =  category_df)
    coxpval <- c(coef(summary(res.cox))[,5],'NA')
    cox_df[gene,] <<- coxpval
    coxHR <- c(coef(summary(res.cox))[1,2],exp(confint(res.cox))[1,])
    cox_HR[gene,] <<- coxHR
  }
}
l <- lapply(colnames(category_df)[5:length(colnames(category_df))], FUN=pval_func, category_df=category_df)
write.csv(cox_HR,paste0(cancer,'_top_gene_',gene,'_cox_HR.csv'))

