library(dplyr)
library(scales)
# library(forestplot)
if ('ggforestplot' %in% installed.packages()==F){
  devtools::install_github("NightingaleHealth/ggforestplot")
}
library(ggforestplot)
library(ggplot2)

project <- 'TCGA'

cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                 'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                 'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
HR_df <- data.frame(matrix(nrow=1,ncol=6))
for (cancer in cancer_list) {
  HR_file <- list.files(paste0('../../',project,'_',cancer,'/05_gene_clinical_prediction_analysis/'),pattern='cox_HR.csv')
  if(cancer=='LIHC'){
    HR_file <- HR_file[4]
  } else if(cancer=='ACC'){
    HR_file <- HR_file[2]
  }
  cluster_HR <- read.csv(paste0('../../',project,'_',cancer,'/05_gene_clinical_prediction_analysis/',HR_file))
  pval <- read.csv(paste0('../../',project,'_',cancer,'/05_gene_clinical_prediction_analysis/',cancer,'_significant_pval_differential_genes_PFI.csv'))[1,3]
  cluster_HR$pval <- pval
  HR_df <- rbind(HR_df,c(cancer,t(cluster_HR[1,])))
}
colnames(HR_df) <- c('Study','Gene','HR','CI_lower','CI_upper','pvalue')
HR_df$HR <- as.numeric(HR_df$HR)
HR_df$CI_lower <- as.numeric(HR_df$CI_lower)
HR_df$CI_upper <- as.numeric(HR_df$CI_upper)
HR_df <- HR_df[-1,]
# true coefficient is the log of the Hazard ratio
HR_df$beta <- log(HR_df$HR)
HR_df$se <- -(log(HR_df$CI_lower) - HR_df$beta)/1.96
HR_df$Study_Gene <- paste0(HR_df$Study,' (',HR_df$Gene,')')

png('../../single_gene_figures/TCGA_single_gene_forest_plot.png',width=1000,height=800)

f <- forestplot(HR_df, name = Study_Gene, estimate=beta, se = se, logodds=T,
           xlab='Hazard Ratio with 95% CI', title = "Hazard Ratio of Top PFI-Prognostic Genes",
           xlim = c(0.01,100)) + scale_x_continuous(trans = "log10",labels=comma)
# adjust format
f$theme$title <- element_text(family='Helvetica',size=27)
f$theme$axis.text.x$family <- 'Helvetica'
f$theme$axis.text.x$size <- 25
f$theme$axis.title.x$size <- 25
f$theme$axis.text.y$family <- 'Helvetica'
f$theme$axis.text.y$size <- 25
f$theme$plot.margin[2] <- unit(as.numeric(f$theme$plot.margin[2])+100,"points")
f
dev.off()


