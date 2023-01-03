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
setwd(paste0("../../",project,"_", cancer, "/03_nmf/"))

if ("NMF" %in% installed.packages() == F){
  install.packages('NMF')
}
library(NMF)
library(RColorBrewer)


# read in data
epigenes_filename <- list.files(path='../RNA-seq_datasets/',pattern=paste0(cancer, '_epigenes_log2norm_sd0\\.\\d+'))
epigenes <- read.csv(paste0('../RNA-seq_datasets/',epigenes_filename),row.names=1,header=T)

# create nmf object
# NMF : nonnegative matrix factorization. All matrices are nonnegative
      # X = WH, where X is the target matrix (nxp) (n genes and p patients), 
      # W is the basis/features matrix (matrix of metagenes/features, dim(nxr), entries (i,j)  = gene i's rank of importance to metagene j, 
      # H is the mixture coefficient matrix (metagene expression profiles dim(rxp), entry (i,j) is metagene i's expression in patient j)
      # Theoretically, X = WH
      # W and H are optimized by minimzing an objective (cost) function that evaluates the estimation accuracy of WH to X
      # The loss function under the lee methods is Euclidean Distance.
      # The output of nmf is a basis matrix W, with the r features and n genes' contribution to the features,
      # and coefficient matrix H, with p patients and r features' expressions in the patients.
      # WH should approximate X, if nmf is performed correctly

# first entry: target matrix - the original dataframe
# second entry: rank - number of clusters ('metagenes')
# third entry: method - form of nmf, default = brunet
# fourth entry: seeding method - default = random
# optional: nrun - number of runs (default is 30)


# do multiple NMF run with different r values to identify best rank
print('performing NMF')
multires <- nmf(epigenes,rank=2:4,method="lee",seed=123456,nrun=100,.options='vp3')
print('saving NMF object')
save(file='multires.Rdata', multires)
# # inspecting data
# # summary of NMF run
# summary(multires)
# # identify different objects inside multires
# attributes(multires)
# # look at fit, each entry in the list multires$fit is a NMVfitX1 object, which is an NMF object for a single rank
# multires$fit



# extract features function to find top NMF features
print('finding top genes in each feature and rank')
find_genes <- function(gene_index, W, res_object,rank){
  feature <<- feature+1
  sorted_genes <- as.data.frame(sort(W[gene_index,feature],decreasing = T))
  colnames(sorted_genes)[1] <- paste0("Feature ",feature)
  write.csv(sorted_genes, file=paste0("Rank_",rank,"/nmf_lee_rank",rank,"_feature",feature,"_genes.csv"),quote=F,row.names=T)
  return(paste0("Rank_",rank,"/nmf_lee_rank",rank,"_feature",feature,"_genes.csv has been written"))
}



# saving all the extract feature outputs
for (r in c(2:4)) {
  print(paste0('Rank: ', r))
  # create folder for rank
  dir.create(paste0('Rank_',r))
  # grab the nmvfit object
  nmf_res <- multires$fit[[r-1]]
  # inspect nmf_res basis matrix
  W <- basis(nmf_res)
  # inspect nmf_res coefficient matrix
  H <- coef(nmf_res)
  # extract features of nmf_res
  e <- extractFeatures(nmf_res)
  # set feature to 0
  # IMPORTANT: DO NOT CHANGE THE NAME OF feature. IT IS A GLOBAL VARIABLE THAT IS MODIFIED BY THE FUNCTION, THE NAME MUST MATCH
  # find top contributin genes for each feature
  feature <- 0
  lapply(e, find_genes, W=W, res_object=nmf_res, rank=r)
  # write cluster membership for each sample
  p <- data.frame(cluster=sort(predict(nmf_res)))
  write.csv(p, file=paste0("Rank_", r, "/nmf_lee_rank",r,"_cluster_membership.csv"), quote=F,row.names=T)
  # plot heatmaps of the consensus matrix
  png(file=paste0("Rank_", r, '/Rank',r,'_consensus_map.png'), width=700, height=700)
  consensusmap(nmf_res, tracks=c('consensus:'),annColors='Dark2', labRow=NA, labCol=NA, Rowv=F, main=paste0('Rank ', r))
  dev.off()
}
# save the measures of the multires. Pay attention to the cophenetic (cophenetic coefficient)
write.csv(multires$measures,file="multi_rank_measures.csv",quote=F,row.names=F)

# plot measures taken by nmf
png(file='multi_rank_measures_plot.png', width=700, height=700)
plot(multires)
dev.off()

