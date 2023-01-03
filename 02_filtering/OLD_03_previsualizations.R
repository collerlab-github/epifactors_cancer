library(argparse)

# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
parser$add_argument("sd_cutoff", nargs=1, help="log2 gene expression standard deviation cutoff")
# parse the arguments
args <- parser$parse_args()
print(args)
cancer <- args$cancer
sd_cutoff <- as.numeric(args$sd_cutoff)
setwd(paste0("../../TCGA_", cancer, "/02_filtering/"))

library(pheatmap)
library(ggplot2)
library(genefilter)
library(readxl)
library(dplyr)



genes <- read.csv(paste0('../RNA-seq_datasets/', cancer, "_log2_75th_percentile>5_counts.csv"),row.names=1, header=T)
normalized_genes <- read.csv(paste0('../RNA-seq_datasets/', cancer, "_75th_percentile>5_counts.csv"), row.names=1, header=T)

# calculate row SD 
gene_sd <- rowSds(as.matrix(genes))
names(gene_sd) <- rownames(genes)


# filter genes with sd cutoff (for log2 genes, and normalized genes)
vargenes <- gene_sd[which(gene_sd>sd_cutoff)]

length(vargenes)

vargenes_df <- genes[which(rownames(genes)%in%names(vargenes)),]

print(paste0('number of sd>', sd_cutoff,' genes:',dim(vargenes_df)[1]))

if (cancer %in% c('BRCA','PRAD','OV','UCEC','CESC','TGCT')==F){
  gencode <- read_excel("../../metadata/gencode.v22.annotation.bed.xlsx",sheet="gencode.v22.annotation.bed",col_names = F)
  gencode <- as.data.frame(gencode)[,c(1,5)]
  # fix gene names
  gencode[,2]<- gsub(';','',gencode[,2])
  # set column names
  colnames(gencode) <- c('Chromosome','HUGO')
  # filter out any sd cutoff gene on y chromosome
  non_ychrom_vargenes <- unique(gencode %>% filter(HUGO %in% rownames(vargenes_df), Chromosome != 'chrY'))$HUGO
  vargenes_df <- vargenes_df[non_ychrom_vargenes,]
}

print(paste0('number of sd>', sd_cutoff, ' genes after y chromosome filter:', dim(vargenes_df)[1]))




# save table
print('saving filtered dataset')
write.csv(vargenes_df,file=paste0("../RNA-seq_datasets/",cancer,"_epigenes_log2norm_sd",sd_cutoff,"_counts.csv"),quote=F)
norm_vargenes_df <- normalized_genes[rownames(vargenes_df),]
write.csv(norm_vargenes_df, file=paste0("../RNA-seq_datasets/",cancer,"_epigenes_norm_sd",sd_cutoff,"_counts.csv"),quote=F)

# # find median of normalized var genes
# vargenes_df_norm_meds <- sort(round(apply(vargenes_df_norm,MARGIN = 1, FUN=median),digits=3))
# write.csv(vargenes_df_norm_meds, file=paste0("../RNA-seq_datasets/", cancer, "_median_normalized_counts.csv"), quote=F)

# 
# # principal component analysis using PCAtools
# if ('PCAtools' %in% installed.packages() ==F){
#   BiocManager::install('PCAtools')
# }
# library(PCAtools)
# 
# # load DESeq vst transformed counts data (already filtered for epigenes), currently stored in vargenes_df object
# vst <- vargenes_df
# # run pca
# p <- pca(vst, removeVar=0.1)
# 
# # screeplot
# screeplot(p, axisLabSize=18, titleLabSize=22,components = getComponents(p, 1:20), title=paste0(cancer, ' SCREE Plot'))
# ggsave(paste0(cancer, '_epigenes_screeplot.png'))
# 
# # bi-plot
# # with just first two pc's
# biplot(p, title = paste0(cancer, ' Biplot'))
# ggsave(paste0(cancer, '_pca.png'))
# biplot(p, showLoadings=T, labSize=5, lab=NULL, title=paste0(cancer, ' Biplot with Loadings'))
# ggsave(paste0(cancer, '_biplot_loadings.png'))
# 
# # pairs plot
# pairsplot(p, title=paste0(cancer, ' Pairs Plot'))
# ggsave(paste0(cancer, '_pairsplot.png'))
# 
# # loadings plot
# plotloadings(p, labSize = 3, title=paste0(cancer, ' Loadings Plot'))
# ggsave(paste0(cancer, '_loadingsplot.png'))
# 
# # optimal number of PC's to retain
# # horn method
# horn <- parallelPCA(vst)
# horn$n
# # elbow method
# elbow <- findElbowPoint(p$variance)
# elbow
# 
# screeplot(p, components = getComponents(p, 1:(horn$n +5)), vline = c(horn$n, elbow), title=paste0(cancer,' SCREE Plot')) +
#   geom_label(aes(x=horn$n - 1, y=50, label='Horn\'s', vjust=-1, size=8)) +
#   geom_label(aes(x=elbow+1, y=50, label='Elbow method', vjust=-1, size=8))
# ggsave(paste0(cancer, '_scree_optimal_pcs.png'))
