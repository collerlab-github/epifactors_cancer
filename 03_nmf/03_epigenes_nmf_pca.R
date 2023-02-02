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
setwd(paste0("../../",project,"_", cancer, "/03_nmf/"))

library(PCAtools)
library(stringr)
library(dplyr)


# read in data
print('reading log2 normalized patient expression data')
epigenes_filename <- list.files(path='../RNA-seq_datasets/',pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'))
epigenes <- read.csv(paste0('../RNA-seq_datasets/',epigenes_filename),row.names=1,header=T)
# sd cutoff (use later for naming convention)
sd_cutoff <- str_extract(string=epigenes_filename, pattern='0\\.\\d+')


# order epigenes samples according to nmf consensus cluster memberships
sample_cluster <- read.csv(paste0("Rank_",rank,"/nmf_lee_rank",rank,"_cluster_membership.csv"), row.names=1, header=T)
if(cancer%in% c('LIHC','NBL')){
  sample_cluster <- sample_cluster %>% arrange(desc(cluster))
} else{
  feature_order <- 1:rank
}
colnames(sample_cluster)[1] <- 'Cluster'
sample_cluster[,1] <- factor(sample_cluster[,1], levels=unique(sample_cluster$Cluster))
epigenes <- epigenes[,rownames(sample_cluster)]

# pca
print('Performing PCA')
p <- pca(epigenes,metadata=sample_cluster, removeVar=0.1)


# bi-plot
# with just first two pc's
print('Generating PCA visualizations')
biplot(p, title = paste0(cancer," NMF PCA"),titleLabSize=25, colby='Cluster', lab=NULL, 
       legendPosition = 'right', legendLabSize=18, legendTitleSize = 15, encircle=T, axisLabSize=20)
ggsave(paste0('Rank_',rank,'/NMF_rank',rank,'_pca_1+2_biplot.png'))
biplot(p, showLoadings=T, title = paste0(cancer, " NMF PCA"), titleLabSize=30, colby='Cluster', lab=NULL, 
       legendPosition = 'right', legendLabSize=18, legendTitleSize = 20, axisLabSize=25, encircle=T)
ggsave(paste0('Rank_',rank,'/NMF_rank',rank,'_pca_1+2_biplot_loadings.png'))

