library(argparse)

# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("project",nargs=1, help="type of project (TCGA or TARGET)")
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
# parse the arguments
args <- parser$parse_args()
print(args)
project <- args$project
cancer <- args$cancer
if (project %in% c('TCGA','TARGET') ==F){
  print('project must be TCGA or TARGET')
  quit(save='no')
}

setwd(paste0("../../",project,"_", cancer, "/05_gene_clinical_prediction_analysis"))


library(readxl)
library(survminer)
library(survival)
library(dplyr)

# read in epigene data
print('Reading patient epigene expression data')
epigenes_filename <- list.files(path='../RNA-seq_datasets/',pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'))
sample_cluster <- as.data.frame(t(read.csv(paste0('../RNA-seq_datasets/',epigenes_filename),row.names=1,header=T)))
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
  cancer_surv <- cbind(cancer_surv, sample_cluster[cancer_surv$bcr_patient_barcode,])
  rownames(cancer_surv) <- cancer_surv$bcr_patient_barcode
  # specify endpoint
  end <- 'PFI'
  print(paste0('Endpoint: ', end))
  
} else{ # project is TARGET
  # save patient names with added cbioportal tag
  print('Adding cbioportal tag to patient names')
  cbioportal <- read.csv('../../TARGET_metadata/TARGET_cbioportal_name.csv', row.names=1, header=T)
  study_id <- cbioportal[cancer,'cbioportal.name']
  patient_names <- rownames(sample_cluster)[which(rownames(sample_cluster) %in% rownames(patient_info))]
  patient_names <- gsub('\\.','-',patient_names)
  patient_names <- gsub('(\\d\\d)(A)','\\1',patient_names)
  patients <- paste0(study_id,':',patient_names)
  
  # survival data is in patient info dataframe already
  # extract samples from survival data
  cancer_surv <- patient_info[which(rownames(patient_info)%in%rownames(sample_cluster)),c("OS.time","OS","age","gender")]
  # add gene expression data
  cancer_surv <- cbind(cancer_surv,sample_cluster[rownames(cancer_surv),])
  # specify endpoint
  end <- 'OS'
  print(paste0('Endpoint: ', end))
  
}
# write patient names to file
write.csv(patients, file='survival_patients.csv',quote=F,row.names = F)
rm(patients)


# seperate patients into two expression groups for each gene
print('Performing survminer patient categorization by expression')
cutpoint_df <- surv_cutpoint(cancer_surv, time=paste0(end,".time"), event=end, variables = colnames(sample_cluster))
category_df <- surv_categorize(cutpoint_df)
# combine cancer_surv age, gender and patient information with the expression group data
category_df <- cbind(cancer_surv[rownames(category_df),c("age","gender")],category_df)


# Initialize p-value dataframes that defines significance between expression groups for each gene
print('Finding significance between expression groups for each gene')
pval_df <- data.frame(row.names = colnames(category_df)[5:length(colnames(category_df))],pval=rep(0,length(colnames(category_df))-4)) 
# initialize cox regression pvalue df and hazard ratio effect size
cox_df <- data.frame(row.names = colnames(category_df)[5:length(colnames(category_df))],cluster=rep(0,length(colnames(category_df))-4), age=0, gender=0)
cox_HR <- data.frame(row.names = colnames(category_df)[5:length(colnames(category_df))],cluster=rep(0,length(colnames(category_df))-4), age=0, gender=0) 
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
    coxHR <- c(coef(summary(res.cox))[,2],'NA','NA')
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
    coxHR <- coef(summary(res.cox))[,2]
    cox_HR[gene,] <<- coxHR
  }
  else{ # run cox regression on just cluster and age
    res.cox <- coxph(s ~ category_df[,gene] + age, data =  category_df)
    coxpval <- c(coef(summary(res.cox))[,5],'NA')
    cox_df[gene,] <<- coxpval
    coxHR <- c(coef(summary(res.cox))[,2],'NA')
    cox_HR[gene,] <<- coxHR
  }
}

# collect pval difference between the 2 groups of each gene
l <- lapply(colnames(category_df)[5:length(colnames(category_df))], FUN=pval_func, category_df=category_df)

# sort by ascending pvalue
pval_df <- pval_df[order(pval_df$pval), , drop=F]

# add benjamini hochberg adjusted p-value
pval_df$adjpval <- p.adjust(pval_df$pval, method='BH')
cox_df$adjpval_cluster <- p.adjust(cox_df$cluster, method='BH')
cox_df$adjpval_age <- p.adjust(cox_df$age, method='BH')
cox_df$adjpval_gender <- p.adjust(cox_df$gender, method='BH')
# round pvals to 2 decimals
# pval_df$adjpval <- round(pval_df$adjpval, digits=7)

# define gene lists to group based on effect on clinical outcome
surv_promoting_genes <- c()
surv_inhibiting_genes <- c()


# function that fits patients by gene and plots endpoint data
fit_func <- function(gene, category_df, pval_df){
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
    surv_promoting_genes <<- c(surv_promoting_genes,g)
  }
  else{
    # gene associated with inhibiting survival outcome (high expression is bad)
    surv_inhibiting_genes <<- c(surv_inhibiting_genes,g)
  }
  adjp <- pval_df[gene, 'adjpval']
  p <- pval_df[gene, 'pval']
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
                      label = paste0('adj_p=', signif(adjp,digits=5)), size = 10)
  ggsave(file = paste0(cancer, '_', end, "_survival_by_", gene, "_expression.png"), survp$plot, width=7, height=7)
}

# plot significant p-value survival curves
print('Saving significant epigenes')
sig_pval <- pval_df[which(pval_df$adjpval<0.05),]
sig_pval <- sig_pval[order(sig_pval$adjpval),]
write.csv(sig_pval, file = paste0(cancer,"_significant_pval_differential_genes_",end,".csv"), quote=F, row.names=T)
  # save respective cox pvalues for age and gender for genes with significant survival differences
sig_cluster <- cox_df[rownames(sig_pval), c('cluster','adjpval_cluster')]
write.csv(sig_cluster, file = paste0(cancer,"_significant_cox_cluster_differential_genes_",end,".csv"), quote=F, row.names=T)

sig_age <- cox_df[rownames(sig_pval),c("age","adjpval_age")]
write.csv(sig_age, file = paste0(cancer,"_significant_cox_age_differential_genes_",end,".csv"), quote=F, row.names=T)

sig_gender <- cox_df[rownames(sig_pval),c("gender","adjpval_gender")]
write.csv(sig_gender, file = paste0(cancer,"_significant_cox_gender_differential_genes_",end,".csv"), quote=F, row.names=T)
  # save cox hazard ratios
sig_HR <- cox_HR[rownames(sig_pval),]
write.csv(sig_HR, file = paste0(cancer,"_significant_cox_hazard_ratios_differential_genes_",end,".csv"), quote=F, row.names=T)

# quit if no significant predictive genes
if (nrow(sig_pval)==0){
  print('No significant epigenes...quitting')
  quit(save='no')
}


print('performing survival analysis based on the expression groups')
l <- lapply(rownames(sig_pval),FUN=fit_func, category_df=category_df, pval_df=pval_df)
# write significant epigenes
sig_epigenes_df <- as.data.frame(t(sample_cluster[,rownames(sig_pval)]))
write.csv(sig_epigenes_df, file=paste0('../RNA-seq_datasets/',cancer,'_',tolower(end),'_epigenes_log2_counts.csv'))

# find enriched complexes and functional groups for significant genes
print('Finding survival promoting and inhitibing epigenes and complexes their involved in')
complexes <- read.delim('../../metadata/epigenes.csv')
complexes$Function <- gsub(',', ';', complexes$Function)
complexes$Protein.complex <- gsub(',', ';', complexes$Protein.complex)

enrichments <- function(genes, complex_df, outcome, gene_effect){
  if(gene_effect %in% c('promoting','inhibiting')==F){
    return('gene_effect must either be "promoting" or "inhibiting"')
  }
  # find main protein complexes and functional types
  # subset epigene complex df with top_genes and order by protein complex
  complex_subset <- complex_df %>% filter(HGNC.approved.symbol%in%genes) %>% select(HGNC.approved.symbol,Function, Protein.complex)
  complex_subset_prot_order <- arrange(complex_subset, desc(Protein.complex))
  true_prot.complex <- complex_subset_prot_order %>% filter(Protein.complex!='#') 
  if(nrow(true_prot.complex)!=0){
    complex_subset_prot_order[1:nrow(true_prot.complex),] <- arrange(true_prot.complex,Protein.complex)
  }
  # find enriched functinal groups
  functional_groups <- aggregate(HGNC.approved.symbol~Function, data=complex_subset,FUN=paste,collapse='; ')
  
  # write tables
  write.csv(complex_subset_prot_order, file=paste0(outcome, '_', gene_effect,'_genes_protein_complexes.csv'),quote=F,row.names=F)
  write.csv(functional_groups, file=paste0(outcome, '_', gene_effect, '_genes_functional_groups.csv'),quote=F,row.names=F)
  
}

enrichments(surv_inhibiting_genes, complexes, tolower(end), 'inhibiting')
enrichments(surv_promoting_genes, complexes, tolower(end), 'promoting')

