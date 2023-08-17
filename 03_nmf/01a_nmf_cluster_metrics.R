
project <- 'TARGET' # TCGA or TARGET
if (project=='TCGA'){
  cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                   'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                   'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
  end <- 'PFI'
} else if (project=='TARGET'){
  cancer_list <- c('AML','NBL','OS','WT')
  end <- 'OS'
} else{
  quit(save='no')
}
# setwd(paste0("../../",project,"_", cancer, "/03_nmf/"))

if ("clValid" %in% installed.packages() == F){
  install.packages('clValid')
}
library(clValid)
library(dplyr)
library(RColorBrewer)


# read in metric dataframe
ranks <- c(2:4)

# load expression data
print('reading log2 normalized patient expression data')
epigenes_filename <- list.files(path='../RNA-seq_datasets/',pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'))
epigenes <- read.csv(paste0('../RNA-seq_datasets/',epigenes_filename),row.names=1,header=T)
dist_epi <- dist(t(epigenes), method="euclidean")
idx_df <- data.frame(rank=ranks, dunn=rep(0,3),silhouette=0,connectivity=0, row.names = ranks)
# calculate dunn for each rank
for (r in ranks){
  
  print(r)
  cluster_membership <- read.csv(paste0("Rank_",r,"/nmf_lee_rank",r,"_cluster_membership.csv"), row.names=1, header=T)
  
  # print(cluster_membership)
  # dunn: ratio of distances: smallest distance across clusters : largest distance within cluster (want smaller denom and bigger numer --> larger dunn index)
  dunn_idx <- dunn(clusters = cluster_membership[colnames(epigenes),], distance = dist_epi)
  # silhouette: # average score of how each data point is similar to its own cluster compared to other clusters. [-1,1] range, closer to 1 mean more similar within cluster than outside of cluster (better)
  sil <- summary(silhouette(x = cluster_membership[colnames(epigenes),], dist = dist_epi))$avg.width
  # connectivity
  conn <- connectivity(dist=dist_epi, clusters = cluster_membership[colnames(epigenes),])
  idx_df[idx_df$rank==r,2:4] <- c(dunn_idx, sil, conn)
}


