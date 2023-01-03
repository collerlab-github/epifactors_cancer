library(dplyr)
library(readxl)

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")



# define inputs
project <- 'TCGA'
cancer <- 'LIHC'
module <- 'turquoise'
setwd(paste0("../../",project,"_", cancer, "/WGCNA/"))



# gencode file contains bed file coordinates
gencode <- read_excel("../../metadata/gencode.v22.annotation.bed.xlsx",sheet="gencode.v22.annotation.bed", range = cell_cols(4:5),col_names = F)
gencode <- as.data.frame(gencode)

# fix gene names
gencode[,1]<- gsub(';','',gencode[,1])
gencode[,2] <- gsub(';','',gencode[,2])

colnames(gencode) <- c('ENSG','HUGO')
rownames(gencode) <- gencode$ENSG
