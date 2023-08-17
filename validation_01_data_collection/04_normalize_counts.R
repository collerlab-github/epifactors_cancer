
library(argparse)

# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("project", nargs=1, help="project name (TCGA or TARGET)")
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
# parse the arguments
args <- parser$parse_args()
print(args)
project <- args$project
cancer <- args$cancer

setwd(paste0("../../",project,"_", cancer, "/02_filtering/"))

# BiocManager::install("DESeq2")
library(DESeq2)
library(readxl)
# read in raw counts
raw_counts <- read.csv(paste0("../RNA-seq_datasets/",cancer,"_raw_counts.csv"),row.names=1,header=T)
colnames(raw_counts) <- gsub('\\.','-',colnames(raw_counts))
# read in metadata
metadata <- read.delim(paste0("../01_data_collection/purity_gdc_sample_sheet_",cancer,".txt"), row.names="Sample.ID", header=T,sep="\t")

# rearrange columns of counts data to match metdata
raw_counts <- raw_counts[,rownames(metadata)]

# make DESeq object
dds <-  DESeqDataSetFromMatrix(countData = raw_counts,
                               colData = metadata,
                               design = ~ 1)
# get normalization factors
print('Outputting normalized counts')
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalize=TRUE)

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

normalized_counts <- filter_func(as.data.frame(normalized_counts))

write.csv(normalized_counts, file=paste0("../RNA-seq_datasets/",cancer,"_normalized_counts.csv"),quote=F, row.names=T)


# # variance stabilize transformation (log2) to remove the dependence of the variance on the mean, 
# # particularly the high variance of the logarithm of count data when the mean is low.
# # This is done by normalizing according to the trend of patient variance over patient mean accross all genes in the experiment.
# # http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations
# vsd <- vst(dds)

# log2(x+1) transformation
print('Outputting Log2 transform normalized counts')
vsd <- log2(normalized_counts+1)


# write.table(normalized_counts,"RNA-seq_datasets/brca_normalized_counts.txt", sep="\t",quote=F)
write.csv(vsd,file=paste0("../RNA-seq_datasets/",cancer,"_log2normalized_counts.csv"),quote=F)
