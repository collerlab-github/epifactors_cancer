# BiocManager::install("DESeq2")
library(DESeq2)
library(readxl)
library(dplyr)

# DESEQ2 vignette: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts
# read in raw counts
hugo_df <- read.csv(paste0("../RNA-seq_datasets/",cancer,"_raw_counts.csv"),row.names=1,header=T)
colnames(raw_counts) <- gsub('\\.','-',colnames(hugo_df))

# read in metadata
sample_cluster <- read.csv(paste0("../03_nmf/Rank_",rank,"/nmf_lee_rank",rank,"_cluster_membership.csv"),row.names=1, header=T)
sample_cluster[,1] <- as.factor(sample_cluster[,1])
colnames(sample_cluster) <- c('Cluster')
event_group <- read.csv(paste0('../../metadata/',project,'_nmf_cluster_',event,'_prognostic_groups.csv'),row.names=1,header=T)
sample_cluster[,1] <- as.factor(event_group[cancer,sample_cluster$Cluster])
# match order of samples with metdata
hugo_df <- hugo_df[,rownames(sample_cluster)]


# create DESeq object
# coldata: metadata dataframe
# design: Column name containing the conditions you want to compare (ex treatment vs control)
dds <- DESeqDataSetFromMatrix(countData=hugo_df, 
                              colData=sample_cluster, 
                              design=~Cluster)
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
# volcano plot
library(ggplot2)
# plot adding up all layers we have seen so far
ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(padj), col=DE)) +
  geom_point() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red", linetype='dashed') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='dashed') +
  ggtitle(paste0(project,' ',cancer,' Fold Change '))
