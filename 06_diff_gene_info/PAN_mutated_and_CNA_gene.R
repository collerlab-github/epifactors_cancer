args <- commandArgs(trailingOnly=TRUE)
if (length(args)<=1){
  stop('Usage: Rscript 01_raw_counts_dataframe.R [cancer1] [cancer2] ...\nAt least 2 cancers must be supplied (i.e.)')
} else{
  print(args)
}

setwd('../../TCGA_PAN/gene_predictive_clinical_prediction_analysis/')

library(stringr)
# get mutated and cnv info in each cancer for overlapping genes

# get overlapping genes
cancer_genes <- read.csv(paste0(paste(args,collapse='_'), '_overlapping_genes.csv'),row.names=1)
common_genes <- str_split(string=cancer_genes[nrow(cancer_genes),'gene'],pattern = '; ')[[1]]
# biomart dataframe for genes
biomart <- read.csv('../../metadata/epigenes_biomart.csv')
biomart$Gene.description <- gsub(' \\[.*\\]','',biomart$Gene.description)
biomart$Gene.description <- gsub(',', ';', biomart$Gene.description)
bm_df <- biomart %>% filter(Gene.name%in%common_genes)
write.csv(bm_df, paste0(paste(args,collapse='_'),'_biomart.csv'),quote=F)


# define dataframe for the information
mut_df <- data.frame(row.names=common_genes)
cna_df <- data.frame(row.names=common_genes)

mc_info <- function(cancer){
  # read mutation data for overlapping genes
  mut_data <- read.csv(paste0('../../TCGA_',cancer,'/06_diff_gene_info/PFI_mutated_genes.csv'))
  mut_genes <- mut_data %>% filter(Gene %in% rownames(mc_df))
  rownames(mut_genes) <- mut_genes$Gene
  # get fraction of patients with mutation in each gene
  mut_genes <- mut_genes %>% mutate(frac = paste0(Num_Mutant_Patients,'/',Total_Num_Patients, ' patients'))
  # add fraction of gene mutation for cancer patients, NA if gene not mutated
  mut_df$cancer <<- ifelse(rownames(mut_df) %in% mut_genes$Gene, mut_genes[rownames(mut_df),'frac'],'NA')
  colnames(mut_df)[ncol(mut_df)] <<- cancer 
  
  # read cna data for overlapping genes
  # read mutation data for overlapping genes
  cna_data <- read.csv(paste0('../../TCGA_',cancer,'/06_diff_gene_info/PFI_cna_genes.csv'))
  cna_genes <- cna_data %>% filter(Gene %in% rownames(mc_df))
  # get fraction of patients with cna in each gene
  cna_genes <- cna_genes %>% mutate(frac = paste0(CNA,' ',Num_CNA_Patients, '/', Total_Num_Patients))
  # aggregate gene's cna types
  cna_genes1 <- aggregate(frac~Gene, data=cna_genes, FUN=paste, collapse='; ', 'patients')
  rownames(cna_genes1) <- cna_genes1$Gene
  # add fraction of gene cna for cancer patients, NA if gene not there
  cna_df$cancer <<- ifelse(rownames(cna_df) %in% cna_genes1$Gene, cna_genes1[rownames(cna_df),'frac'],NA)
  colnames(cna_df)[ncol(cna_df)] <<- cancer 
}

l <- lapply(args, FUN=mc_info)

write.csv(cna_df, paste0(paste(args, collapse='_'),'_cna_genes.csv'), quote=F)
write.csv(mut_df, paste0(paste(args, collapse='_'),'_mutated_genes.csv'), quote=F)
