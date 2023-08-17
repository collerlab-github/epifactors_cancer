library(argparse)

# create parser object to add command line arguments
# parser <- ArgumentParser()
# # specify desired options
# # by default ArgumentParser will add a help option
# parser$add_argument("project",nargs=1, help="type of project (TCGA or TARGET)")
# parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
# parser$add_argument("gene", nargs=1, help='gene by which to stratify patients (default is NA)')
# # parse the arguments
# args <- parser$parse_args()
# print(args)
# project <- args$project
# cancer <- args$cancer
project <- 'TCGA'

# cancer_list <- c('BRCA','THCA','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
#                  'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
#                  'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
cancer <- 'PAAD'
dir.create(paste0("../../",project,"_", cancer, "/EMTome"), showWarnings = F)
setwd(paste0("../../",project,"_", cancer, "/EMTome"))


library(readxl)
library(survminer)
library(survival)
library(dplyr)

# read in epigene data
print('Reading patient epigene expression data')
# epigenes_filename <- list.files(path='../RNA-seq_datasets/',pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'))
sample_cluster <- as.data.frame(read.csv('../../metadata/EMTome_EMT_scores.csv', row.names=1))
sample_cluster <- sample_cluster[which(sample_cluster$cancer == cancer),"emt",F]
patient_info <- read.csv(file=paste0('../03_nmf/Rank_2','/patient_clinical_info.csv'),row.names=1,header=T)

if (sum(is.na(sample_cluster)) / nrow(sample_cluster) > 0.1){
  print(paste0('Too many NA: ', table(is.na(sample_cluster))['TRUE']))
}
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
  # print('Adding cbioportal tag to patient names')
  # cbioportal <- read.csv('../../metadata/TCGA_cbioportal_name.csv', row.names=1, header=T)
  # study_id <- cbioportal[cancer,'cbioportal.name']
  # patients <- paste0(study_id,':',rownames(sample_cluster)[which(rownames(sample_cluster) %in% TCGA_surv$bcr_patient_barcode)])
  # 
  
  # extract samples from survival data
  cancer_surv <- TCGA_surv[which(TCGA_surv$bcr_patient_barcode %in% rownames(sample_cluster)),]
  # add age and gender
  cancer_surv$age <- sample_age[cancer_surv$bcr_patient_barcode,]
  cancer_surv$gender <- sample_gender[cancer_surv$bcr_patient_barcode,]
  # add gene_expression data
  cancer_surv <- cbind(cancer_surv, sample_cluster[cancer_surv$bcr_patient_barcode,,F])
  rownames(cancer_surv) <- cancer_surv$bcr_patient_barcode
}

# seperate patients into two expression groups for each gene and endpoint

for (end in c('OS','DSS','PFI')){
  print('Performing survminer patient categorization by expression')
  cutpoint_df <- surv_cutpoint(cancer_surv, time=paste0(end,".time"), event=end, variables = colnames(sample_cluster))
  category_df <- surv_categorize(cutpoint_df)
  # combine cancer_surv age, gender and patient information with the expression group data
  category_df <- cbind(cancer_surv[rownames(category_df),c("age","gender")],category_df)
  
  write.csv(category_df,paste0(cancer,'_emt_',end,'_groups.csv'))
  
}


print('Finding significance between expression groups for each gene')
pval_df <- data.frame(row.names = colnames(category_df)[5:length(colnames(category_df))],pval=rep(0,length(colnames(category_df))-4)) 
# initialize cox regression pvalue df and hazard ratio effect size
endpoint <- c('OS','DSS','PFI')
pval_df <- data.frame(row.names = colnames(category_df)[5:length(colnames(category_df))], OS=rep(0,length(colnames(category_df))-4), DSS = rep(0,length(colnames(category_df))-4), PFI=rep(0,length(colnames(category_df))-4))
cox_df <- data.frame(row.names=c('cluster','age','gender'), OS = rep(0,3), DSS=rep(0,3), PFI = rep(0,3))
cox_HR <- data.frame(row.names = c('cluster','CI_low','CI_high'),OS = rep(0,3), DSS=rep(0,3), PFI = rep(0,3)) 
# function that finds pvalue of difference btwn high and low expression categories
# and performs cox regression on age and gender
pval_func <- function(gene, category_df, end){
  # print(gene)
  # create survival object (survival_time, survival_status)
  s <- Surv(category_df[,paste0(end,'.time')],category_df[,end])
  
  # find p-value difference btwn gene of interest's high/low expression groups
  diff <- survdiff(s~category_df[,gene],data=category_df)
  pval <- pchisq(diff$chisq,df=length(diff$n)-1,lower.tail=F)
  # store p-value in dataframe
  pval_df[gene,end] <<- pval
  
  
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
    cox_df[,end] <<- coxpval
    coxHR <- c(coef(summary(res.cox))[1,2],exp(confint(res.cox))[1,])
    cox_HR[,end] <<- coxHR
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
    cox_df[,end] <<- coxpval
    coxHR <- c(coef(summary(res.cox))[1,2],exp(confint(res.cox))[1,])
    cox_HR[,end] <<- coxHR
  }
  else{ # run cox regression on just cluster and age
    res.cox <- coxph(s ~ category_df[,gene] + age, data =  category_df)
    coxpval <- c(coef(summary(res.cox))[,5],NA)
    cox_df[,end] <<- coxpval
    coxHR <- c(coef(summary(res.cox))[1,2],exp(confint(res.cox))[1,])
    cox_HR[,end] <<- coxHR
  }
}


emt_direction <- data.frame(row.names = colnames(category_df)[5:length(colnames(category_df))], OS=rep(NA,length(colnames(category_df))-4), DSS = rep(NA,length(colnames(category_df))-4), PFI=rep(NA,length(colnames(category_df))-4))
# function that fits patients by gene and plots endpoint data
fit_func <- function(gene, category_df, pval, end){
  g <<- gene
  s <<- Surv(category_df[,paste0(end,'.time')],category_df[,end])
  # make survival fit object grouped by the gene column of category_df (high/low expression groups)
  fit <<- survfit(s~category_df[,g], data=category_df)
  # specify which expression group is associated with worse clinical outcome
  # compare average survival prob in each group
  high_n <- fit$strata[1]
  low_n <- fit$strata[2]
  high_surv_prob <- mean(fit$surv[1:fit$strata[1]]) # first entry in fit$strata is the number of high expression patients
  low_surv_prob <- mean(fit$surv[(fit$strata[1]+1):sum(fit$strata)]) # first+second = index of the last low expression patients
  
  # append gene to survival-promoting or inhibiting group
  if (high_surv_prob> low_surv_prob){
    # gene associated with promoting survival outcome (high expression is good)
    emt_direction[gene,end] <<- 'high_pfi'
  }
  else{
    # gene associated with inhibiting survival outcome (high expression is bad)
    emt_direction[gene,end] <<- 'low_pfi'
  }
  
  # plot df
  # plot curves
  group_lab <- c(paste0('High Expr (n=',high_n,')'),paste0('Low Expr (n=',low_n,')'))
  # specify y label
  if (end=='OS'){
    ylabel <- 'Survival (overall) (%)'
  } else if (end=='DSS'){
    ylabel <- 'Survival (disease specific) (%)'
  } else if (end=='PFI'){
    ylabel <- 'Progression Free (%)'
  }
  survp <- ggsurvplot(fit, title=paste0(cancer,' ',end,' by ', g), font.main=30,
                      font.x=20, font.y=20, font.tickslab=14, font.legend=20,
                      xlab='Time (days)', ylab=ylabel,
                      legend=c(0.8,1), legend.labs = group_lab, legend.title='',
                      text=element_text(size=20)
  )
  survp$plot <- survp$plot + 
    ggplot2::annotate("text", 
                      # x and y coordinates of the text
                      x = max(na.omit(category_df[,paste0(end,'.time')]))-(0.3*max(na.omit(category_df[,paste0(end,'.time')]))), y = 0.8, 
                      label = paste0('p=', signif(pval,digits=5)), size = 10)
  ggsave(file = paste0(cancer, '_', end, "_survival_by_", gene, "_expression.png"), survp$plot, width=7, height=7)
}

for ( end in c('OS','DSS','PFI')){
  category_df <- read.csv(paste0(cancer,'_emt_',end,'_groups.csv'), row.names=1)
  l <- lapply(colnames(category_df)[5:length(colnames(category_df))], FUN=pval_func, category_df=category_df, end=end)
  print('performing survival analysis based on the expression groups')
  fit_func('emt',category_df = category_df, pval = as.numeric(cox_df['cluster', end]), end=end)
}

write.csv(cox_HR,paste0(cancer,'_emt_cox_HR.csv'))
write.csv(cox_df, paste0(cancer,'_emt_cox_df.csv'))
write.csv(emt_direction, paste0(cancer,'_emt_direction.csv'))




