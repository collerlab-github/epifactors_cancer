
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<=1){
  stop('Usage: Rscript 01_raw_counts_dataframe.R [cancer1] [cancer2] ...\nAt least 2 cancers must be supplied (i.e.)')
} else{
  print(args)
}

project <- 'TARGET'
setwd(paste0('../../',project,'_PAN'))

# BiocManager::install("DESeq2")
library(DESeq2)
# read in raw counts
raw_counts <- read.csv(paste0("RNA-seq_datasets/",paste0(args, collapse='_'),"_raw_counts.csv"),row.names=1,header=T)
# read in metadata
metadata <- read.csv(paste0('01_filtering/',paste(args,collapse='_'),'_patient_data.csv'),row.names=1)
# rearrange columns of counts data to match metdata
raw_counts <- raw_counts[,rownames(metadata)]
# edit column names
colnames(raw_counts) <- gsub('\\.','-',colnames(raw_counts))

# make DESeq object
dds <-  DESeqDataSetFromMatrix(countData = raw_counts,
                               colData = metadata,
                               design = ~ 1)
# get normalization factors
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalize=TRUE)
write.csv(normalized_counts, file=paste0("RNA-seq_datasets/",paste(args,collapse='_'),"_normalized_counts.csv"),quote=F, row.names=T)


# # variance stabilize transformation (log2) to remove the dependence of the variance on the mean, 
# # particularly the high variance of the logarithm of count data when the mean is low.
# # This is done by normalizing according to the trend of patient variance over patient mean accross all genes in the experiment.
# # http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations
# vsd <- vst(dds)

# log2(x+1) transformation
vsd <- log2(normalized_counts+1)

# write.table(normalized_counts,"RNA-seq_datasets/brca_normalized_counts.txt", sep="\t",quote=F)
write.csv(vsd,file=paste0("RNA-seq_datasets/",paste(args,collapse='_'),"_log2normalized_counts.csv"),quote=F)
