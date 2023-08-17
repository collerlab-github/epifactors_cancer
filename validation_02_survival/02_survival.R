library(argparse)

# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("project", nargs=1, help="type of project (e.g TCGA or TARGET")
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
parser$add_argument("rank", nargs=1, help='rank of NMF')
# parse the arguments
args <- parser$parse_args()
print(args)
# args <- list(project='TCGA_val',cancer='LUAD',rank=2) KIRC or LUAD
project <- args$project # TCGA_val
cancer <- args$cancer # LUAD
rank <- as.numeric(args$rank)

setwd(paste0("../../", project, "_", cancer, "/03_nmf_survival/"))
# create folder for rank

library(survival)
library(readxl)
library(survminer)
library(dplyr)

# if you want pairwise comparison, set to false
groupwise=F
# read in patient clinical info and cluster membership
print('reading in NMF cluster')
# sample_cluster <- read.csv(paste0(cancer,'_NMFgenes_GSVA_scores.csv'), row.names=1, header=T)
patient_info <- read.csv(file=paste0('patient_clinical_info.csv'),row.names=1,header=T)
sample_cluster <- patient_info[,1,F]
# read in project survival data
print('Reading Survival Data')
rownames(sample_cluster) <- gsub('-\\d\\d\\w$','',rownames(sample_cluster))
rownames(patient_info) <- gsub('-\\d\\d\\w$','',rownames(patient_info))
# sample age and gender
sample_age <- patient_info[,'age',F]
sample_gender <- patient_info[,'gender',F]
# read survival data
TCGA_surv <- cbind(read_excel('../../metadata/TCGA_survival_data.xlsx',sheet="TCGA-CDR",
                              col_names=T, range= cell_cols(c("B","C"))),
                   read_excel('../../metadata/TCGA_survival_data.xlsx',sheet="TCGA-CDR",
                              col_names=T, range= cell_cols(c("Z","AG"))))
# check for patients unaccounted for in survival data:
missing_surv_patients <- rownames(sample_cluster)[which(rownames(sample_cluster) %in% TCGA_surv$bcr_patient_barcode == F)]
print(paste0('Patients without survival data: ', ifelse(length(missing_surv_patients)==0,'none',missing_surv_patients)))
# extract samples from survival data
cancer_surv <- TCGA_surv[which(TCGA_surv$bcr_patient_barcode %in% rownames(sample_cluster)),]
cancer_surv$cluster <- sample_cluster[cancer_surv$bcr_patient_barcode,]
# add age
cancer_surv$age <- sample_age[cancer_surv$bcr_patient_barcode,]
cancer_surv$gender <- sample_gender[cancer_surv$bcr_patient_barcode,]
rownames(cancer_surv) <- cancer_surv$bcr_patient_barcode
# write to file
write.csv(cancer_surv[,c("OS","OS.time","DSS","DSS.time","PFI","PFI.time","cluster","age","gender")], paste0(project,'_',cancer,'_survival_info.csv'),row.names=T,quote=F)
# filter out patients with no age or gender info
cancer_surv <- cancer_surv %>% filter(!is.na(age),!is.na(gender))

# for rank 3, do groupwise comparisons (e.g 1+2 vs. 3, 1+3 vs. 2, 2+3 vs. 1)
if (rank==3 && groupwise==T){
  cancer_surv$'1+2' <- ifelse(cancer_surv$cluster%in%c(1,2), '1+2','3')
  cancer_surv$'1+3' <- ifelse(cancer_surv$cluster%in%c(1,3), '1+3','2')
  cancer_surv$'2+3' <- ifelse(cancer_surv$cluster%in%c(2,3),'2+3','1')
  cluster_comb <- list(c('1+2','3'),c('1+3','2'),c('2+3',1))
} else{
  cluster_comb <- combn(x=unique(cancer_surv$cluster), m=2, simplify=F)
}

# create survival plot (kaplan meier) for each survival endpoint

endpoint <- c('OS','DSS','PFI')
pval_df <- data.frame(OS=rep(0,length(cluster_comb)), DSS = rep(0,length(cluster_comb)), PFI=rep(0,length(cluster_comb)))
cox_df <- data.frame(OS = rep(0,3), DSS=rep(0,3), PFI = rep(0,3))


# format pvalue and cox df's
sorted_cluster_comb <- lapply(cluster_comb,sort)
rownames(pval_df) <- lapply(sorted_cluster_comb,paste,collapse=':')
rownames(pval_df) <- sort(rownames(pval_df))
rownames(cox_df) <- c('cluster','age','gender')

# plot survival curves for each endpoint
print('Performing survival analysis')
plot_surv <- function(clust, cancer_surv){
  print(paste(sort(clust),collapse=':')) # formats cluster to rownames for pval_df
  for (e in endpoint){
    t <- paste0(e,'.time')
    # subset the respective endpoint
    if (groupwise==T){
      e_data <- cancer_surv[,c(e,t,'cluster','1+2','1+3','2+3')]
    } else{
      e_data <- cancer_surv[,c(e,t,'cluster','age','gender')]
    }
    # # plot data for each pairwise cluster combination
    # lapply(cluster_comb, FUN=plot_surv, e_data=e_data)
    
    if (groupwise==T){
      e_data$cluster <- e_data[,clust[1]]
    }
    # subset the respective endpoint
    e_clust_data <<- e_data[which(e_data$cluster %in% clust),]
    # survival object
    survival_time <<- e_clust_data[,2]
    survival_status <<- e_clust_data[,1]
    s <<- Surv(survival_time, survival_status)
    fit <- survfit(s~cluster,data=e_clust_data) # fits survival curve to data
    
    # # compute cox regression and add pvalues to df
    if(length(unique(e_clust_data$gender))>1){
      res.cox <- coxph(s ~ cluster + age + gender, data =  e_clust_data)
      coxpval <- coef(summary(res.cox))[,5]
      cox_df[,e] <<- coxpval
    }
    else{
      res.cox <- coxph(s ~ cluster + age, data =  e_clust_data)
      coxpval <- c(coef(summary(res.cox))[,5],'NA')
      cox_df[,e] <<- coxpval
    }
    
    # write.csv(coxpval, paste0('Rank_',rank,'/',cancer,"_rank",rank,"_cluster", paste(clust, collapse="_"), "_cox_pval_",e,".csv"))
    # compute p-value of survival curve difference and write to pval_df
    diff <- survdiff(s~cluster,data=e_clust_data)
    pval <- pchisq(diff$chisq,df=length(diff$n)-1,lower.tail=F)
    pval_df[paste(sort(clust),collapse=':'),e] <<- pval
    # pval <- round(pchisq(diff$chisq,df=length(diff$n)-1,lower.tail=F), digits=5)
    coxpval <- formatC(cox_df['cluster',e],format="e", digits=2)
    # number of patients in each cluster
    group_labs <- paste0(sort(clust), ' (n=', c(fit$n[1],fit$n[2]),')')
    # plot curves
    # specify y label
    if (e=='OS'){
      ylabel <- 'Survival (overall) (%)'
    } else if (e=='DSS'){
      ylabel <- 'Survival (disease specific) (%)'
    } else if (e=='PFI'){
      ylabel <- 'Progression Free (%)'
    }
    print('plotting')
    survp <- ggsurvplot(fit, title=paste0('Low Purity ',cancer,' ', e, ' Curve'), font.main=30,
                         font.x=20, font.y=20, font.tickslab=14, font.legend=20,
                         xlab='Time (days)', ylab=ylabel,
                         legend=c(0.8,1), legend.labs = group_labs, legend.title="",
                         text=element_text(size=20)
                         )
    survp$plot <- survp$plot + 
      ggplot2::annotate("text", 
                        x = max(na.omit(survival_time))-(0.2*max(na.omit(survival_time))), y = 0.8, # x and y coordinates of the text
                        label = paste0("p=",coxpval), size = 7)
    ggsave(file = paste0(project,'_',cancer,"_rank",rank,"_", paste(clust, collapse="_"), "_",e,".png"))
  }
}
# plot curves for each pair of clusters
l <- lapply(cluster_comb, FUN=plot_surv, cancer_surv=cancer_surv)
write.csv(pval_df,file = paste0(project,'_',cancer,'_endpoint_pval.csv'),row.names=T, quote=F)
write.csv(cox_df, file = paste0(project,'_',cancer,'_cox_pval.csv'),row.names=T, quote=F)

