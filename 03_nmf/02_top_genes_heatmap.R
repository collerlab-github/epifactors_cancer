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
project <- args$project
cancer <- args$cancer
rank <- as.numeric(args$rank)
setwd(paste0("../../",project,"_", cancer, "/03_nmf/"))

if ('ComplexHeatmap'%in%installed.packages()==F){
  library(devtools)
  install_github("jokergoo/ComplexHeatmap")
}
library(ComplexHeatmap)
# library(pheatmap)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(dplyr)

# read in data
print('reading log2 normalized patient expression data')
epigenes_filename <- list.files(path='../RNA-seq_datasets/',pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'))
epigenes <- read.csv(paste0('../RNA-seq_datasets/',epigenes_filename),row.names=1,header=T)
# sd cutoff (use later for naming convention)
sd_cutoff <- str_extract(string=epigenes_filename, pattern='0\\.\\d+')

# order epigenes samples according to nmf consensus cluster memberships
sample_cluster <- read.csv(paste0("Rank_",rank,"/nmf_lee_rank",rank,"_cluster_membership.csv"), row.names=1, header=T)
if(cancer %in% c('LIHC','NBL')){
  sample_cluster <- sample_cluster %>% arrange(desc(cluster))
  feature_order <- rank:1
} else{
  feature_order <- 1:rank
}
sample_cluster[,1] <- as.factor(sample_cluster[,1])
epigenes <- epigenes[,rownames(sample_cluster)]


# read in top genes for each nmf feature (metagene/geneset)
top_genes <- c() # top genes list
feat <- c() # list of feature identity of each top gene
# looping through each feature index
for (feature in feature_order){
  print(paste0("Rank ", rank, " feature ", feature))
  # append top genes list
  top_genes <- c(top_genes,rownames(read.csv(paste0("Rank_",rank,"/nmf_lee_rank",rank,"_feature",feature,"_genes.csv"),header=T,row.names=1)))
  # append feat list
  feat <- c(feat,rep(as.character(feature),length(top_genes)-length(feat)))
  print(length(top_genes))
}
########################
##### Draw Heatmap #####
########################
# heatmap of top_genes gene expression
# epigene counts matrix of only top_genes
print('Generating heatmap of top epigene expression')
top_df <- epigenes[top_genes,]
# dataframe for annotating genes in heatmap
feature_df <- data.frame(feature = paste0('Gene Set ',feat),row.names=top_genes)

library(circlize)
plot_heatmap <- function(top_df, rank, sd_cutoff, rclust=F, filename){
  # make png
  png(file=filename, units="in", width=15, height=10, res=500)
  # specify color, scaled expression, column names, row_names
  color_fun <- colorRamp2(c(-2,0,2), c('blue','white','red'))
  top_df_scaled <- t(scale(t(as.matrix(top_df))))
  if(cancer %in% c('LIHC','NBL')){
    column_cluster_order <- rank:1
    row_cluster_order <- sort(unique(feature_df$feature), decreasing = T)
  } else {
    column_cluster_order <- 1:rank
    row_cluster_order <- sort(unique(feature_df$feature))
  }
  column_split_names <- factor(paste0('Cluster ', sample_cluster$cluster), levels=paste0('Cluster ',column_cluster_order))
  row_split_names <- factor(feature_df$feature, levels=row_cluster_order)
  # heatmap
  h <- Heatmap(as.matrix(top_df_scaled), name='Expression', col=color_fun,
               heatmap_legend_param = list(labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                           title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                           legend_height = unit(2.5, "in")),
               show_row_names=TRUE, row_names_gp = gpar(fontsize=9, fontfamily='Helvetica'),
               show_column_names=FALSE,
               row_split=row_split_names, row_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
               column_split=column_split_names, column_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
               cluster_row_slices=F, cluster_rows=rclust,
               cluster_column_slices=F, cluster_columns=T,
               row_gap=unit(2,'mm'), column_gap=unit(2,'mm'), border = TRUE
  )
  # add more titles to heatmap
  ht <- draw(h,
       row_title='Top NMF Epigenes', row_title_gp=gpar(fontsize=25, fontfamily='Helvetica'),
       column_title=paste0(cancer,' Patients'), column_title_gp=gpar(fontsize=25, fontfamily='Helvetica'),
  )
  saveRDS(row_order(ht), gsub('.png','_order.RDS',filename))
  dev.off()
}
# plot top genes heatmap with and without row clustering
plot_heatmap(top_df=top_df, rank=rank, sd_cutoff=sd_cutoff, filename=paste0("Rank_",rank,"/Rank",rank,"_log2_sd",sd_cutoff,"_topgenes_heatmap.png"))
plot_heatmap(top_df=top_df, rank=rank, sd_cutoff=sd_cutoff, rclust=T, filename=paste0("Rank_",rank,"/Rank",rank,"_log2_sd",sd_cutoff,"_topgenes_heatmap_row_clust.png"))


