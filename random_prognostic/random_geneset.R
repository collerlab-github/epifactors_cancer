library(readxl)
library(survminer)
library(survival)
library(dplyr)


cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                 'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                 'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
project <- 'TCGA'
end <- 'PFI'
dir.create(paste0("../../Summary_tables/random_prognostic"), showWarnings = F)
setwd(paste0("../../Summary_tables/random_prognostic"))


# read in gencode genename index
gencode <- read_excel("../../metadata/gencode.v22.annotation.bed.xlsx",sheet="gencode.v22.annotation.bed", range = cell_cols(4:5),col_names = F)
gencode <- as.data.frame(gencode)
# fix gene names
gencode[,1]<- gsub(';','',gencode[,1])
gencode[,2] <- gsub(';','',gencode[,2])
colnames(gencode) <- c('ENSG','HUGO')
rownames(gencode) <- gencode$ENSG
# convert ensg to hugo
filter_func <- function(df){
  # add corresponding hugo name column to counts data
  df$HUGO <- gencode[match(rownames(df), rownames(gencode)),2]
  # check
  head(df %>%
         select(HUGO))
  # remove duplicate HUGO ids
  dup_ids <- df$HUGO[duplicated(df$HUGO)]
  df[df$HUGO %in% dup_ids,] %>% select(HUGO)
  hugo_df <- df %>% filter(!duplicated(df$HUGO))
  # set rownames to hugo ids
  hugo_to_ensg <- rownames(hugo_df)
  names(hugo_to_ensg) <- hugo_df$HUGO
  rownames(hugo_df) <- hugo_df$HUGO
  hugo_df <- hugo_df[,-ncol(hugo_df)]
  return(hugo_df)
}

# get all common genes
num_epigenes <- rep(0,length(cancer_list))
names(num_epigenes) <- cancer_list
common_genes <- c()
count <- T
for (cancer in cancer_list){
  # read full gene expression data
  print(cancer)
  epigenes_data <- read.csv(list.files(paste0('../../',project,'_',cancer,'/RNA-seq_datasets/'),pattern='log2norm_sd', full.names=T))
  num_epigenes[cancer] <- nrow(epigenes_data)
  print(paste('Number epigenes:',nrow(epigenes_data)))
  print(all(!is.na(epigenes_data)))
  log2_counts <- read.csv(paste0('../../',project,'_',cancer,'/RNA-seq_datasets/',cancer,'_log2normalized_counts.csv'), row.names=1)
  hugo_df <- filter_func(log2_counts)
  # gene_sds <- apply(hugo_df, 1,sd)
  # candidate_genes <- rownames(hugo_df)
  if (count){
    common_genes <- rownames(hugo_df)
    count <- F
  } else{
    common_genes <- intersect(common_genes, rownames(hugo_df)) # ensures random genes will meet SD cutoff of previous cancers
    # # remove sd cutoff
    # common_genes <- common_genes[gene_sds[common_genes] >= 0.5] # ensures random genes will meet SD cutoff of current cancer
  }
  print(all(common_genes %in% rownames(hugo_df)))
}

# use 720
epigenes <- read.csv('../../metadata/epigene_list.txt', sep='\t', header=F)[,1]
set.seed(0)
# sample random genes that are not epigenes
common_genes <- common_genes[common_genes %in% epigenes ==F]
random_genes <- sample(common_genes, length(epigenes), replace=F)
write.csv(random_genes, 'random_genes.csv', row.names=F)

# perform sd cutoff in surminer
