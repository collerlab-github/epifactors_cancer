library(argparse)

# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("project", nargs=1, help="project name (TCGA or TARGET)")
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
parser$add_argument("rank", nargs=1, help='rank of NMF')

# parse the arguments
args <- parser$parse_args()
print(args)
# cancer type and rank
project <- args$project
cancer <- args$cancer
rank <- as.numeric(args$rank)

setwd(paste0("../../",project,"_", cancer, "/02_filtering/"))

# BiocManager::install("DESeq2")
library(DESeq2)
library(readxl)
library(dplyr)

# read in raw counts
raw_counts <- read.csv(paste0("../RNA-seq_datasets/",cancer,"_raw_counts.csv"),row.names=1,header=T)
colnames(raw_counts) <- gsub('\\.','-',colnames(raw_counts))
# # read in metadata
# metadata <- read.delim(paste0("../01_data_collection/purity_gdc_sample_sheet_",cancer,".txt"), row.names="Sample.ID", header=T,sep="\t")

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

  df2 <- df %>% filter(!duplicated(df$HUGO))

  return(df2)
}

# filter epigenes for log2 and normalized counts
# outfile <- paste0('../RNA-seq_datasets/', cancer, '_hugo_raw_counts.csv')
hugo_df <- filter_func(raw_counts)
hugo_to_ensg <- rownames(hugo_df)
names(hugo_to_ensg) <- hugo_df$HUGO
rownames(hugo_df) <- hugo_df$HUGO
hugo_df <- hugo_df[,-ncol(hugo_df)]

colnames(hugo_df) <- gsub('-','\\.',colnames(hugo_df))


# read nmf consensus cluster memberships
if(project == 'TCGA'){
  event <- 'pfi'
} else if(project == 'TARGET'){
  event <- 'os'
} else{
  print('project is not TCGA or TARGET')
}
sample_cluster <- read.csv(paste0("../03_nmf/Rank_",rank,"/nmf_lee_rank",rank,"_cluster_membership.csv"),row.names=1, header=T)
sample_cluster[,1] <- as.factor(sample_cluster[,1])
colnames(sample_cluster) <- c('Cluster')
event_group <- read.csv(paste0('../../metadata/',project,'_nmf_cluster_',event,'_prognostic_groups.csv'),row.names=1,header=T)
sample_cluster[,1] <- as.factor(event_group[cancer,sample_cluster$Cluster])
# match order of samples with metdata
hugo_df <- hugo_df[,rownames(sample_cluster)]


# create DESeq object
dds <- DESeqDataSetFromMatrix(countData=hugo_df, 
                              colData=sample_cluster, 
                              design=~Cluster)
if(event_group[cancer,2] == 'high_pfi'){
  dds$Cluster<- relevel(dds$Cluster, ref = paste0("low_",event))
}

# run DESeq
dds <- DESeq(dds)

res <- results(dds)
print(res)
res_df <- as.data.frame(res)
res_df$ENSG <- hugo_to_ensg[rownames(res_df)]
res_df <- res_df%>%arrange(padj)

write.csv(res_df,paste0(project,'_',cancer,'_DESeq_output.csv'),row.names=T)

# res_df$DE <- "No"
# fc_cutoff <- 0.6
# res_df$DE[res_df$log2FoldChange > fc_cutoff & res_df$padj < 0.05] <- "Up"
# res_df$DE[res_df$log2FoldChange < -fc_cutoff & res_df$padj < 0.05] <- "Down"
# # volcano plot
# library(ggplot2)
# # plot adding up all layers we have seen so far
# ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(padj), col=DE)) +
#   geom_point() + 
#   scale_color_manual(values=c("blue", "black", "red")) +
#   geom_vline(xintercept=c(-0.6, 0.6), col="red", linetype='dashed') +
#   geom_hline(yintercept=-log10(0.05), col="red", linetype='dashed') +
#   ggtitle(paste0(project,' ',cancer,' Fold Change '))
