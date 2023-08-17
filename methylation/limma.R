library(limma)
library(dplyr)

cancer <- 'LIHC'
setwd(paste0('../../TCGA_',cancer))
betas <- read.csv(paste0('methylation/TCGA_',cancer,'_methylation_beta_values_bc.csv'), row.names=1)
sample_cluster <- read.csv(paste0("03_nmf/Rank_2/nmf_lee_rank2_cluster_membership.csv"), row.names=1, header=T)
rownames(sample_cluster) <- gsub('01A','01',rownames(sample_cluster))
sample_cluster <- sample_cluster[colnames(betas),,F]
colnames(sample_cluster)[1] <- 'Cluster'

event_group <- read.csv('../metadata/TCGA_nmf_cluster_pfi_prognostic_groups.csv', row.names=1, header=T)
# if (event_group[cancer,1]=='high_pfi'){
#   sample_cluster[,1] <- factor(sample_cluster[,1], levels=c(2,1))
# } else {
#   sample_cluster[,1] <- factor(sample_cluster[,1], levels=c(1,2))
#   
# }
sample_cluster$Cluster <- as.character(event_group[cancer,])[sample_cluster$Cluster]
sample_cluster$Cluster <- factor(sample_cluster$Cluster, levels=c('low_pfi','high_pfi')) 
# this is the factor of interest (NMF group)
nmf_group <- sample_cluster$Cluster
# explore data
# plotMDS(betas, top=1000, gene.selection="common", 
#         col=levels(nmf_group))
# legend("right", legend=levels(nmf_group), text.col=levels(nmf_group),
#        bg="white", cex=0.7)

# Limma
# create contrast matrix

# MAY POTENTIALLY NEED TO ACCOUNT FOR STAGE/GRADE
# # this is the individual effect that we need to account for
# individual <- factor(targets$Sample_Source) 

# use the above to create a design matrix
design <- model.matrix(~0+nmf_group, data=betas)
colnames(design) <- c(levels(nmf_group))

# fit the linear model 
fit <- lmFit(betas, design)
# create a contrast matrix for specific comparisons
# first term is the numerator and second term is the denominator when calculating fold change/differential methylation
contMatrix <- makeContrasts(high_pfi-low_pfi,
                            levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))


# # get the table of results for the first contrast (naive - rTreg)
# ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
#                       c(1:4,12:19,24:ncol(ann450k))]
# DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
# read probe-gene map
probe_map <- read.delim('../TCGA_PAN/xena_methylation/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy.txt', sep='\t', header=T, row.names=1)
DMPs <- topTable(fit2, num=Inf, coef=1)
DMPs <- merge(DMPs, probe_map, by = 0, all.x = T)
DMPs <- DMPs %>% arrange(adj.P.Val)
DMPs <- tibble::column_to_rownames(DMPs, "Row.names")
# write.csv(DMPs, file=paste0("methylation/TCGA_",cancer,"_diff_meth_probes.csv"), row.names=FALSE)
# dim(DMPs %>% filter(adj.P.Val <= 0.05))
write.csv(DMPs %>% filter(adj.P.Val <= 0.05, logFC>0), file=paste0("methylation/TCGA_",cancer,"_diff_meth_probes_signif_up_in_highpfi.csv"), row.names=T)
write.csv(DMPs %>% filter(adj.P.Val <= 0.05, logFC<0), file=paste0("methylation/TCGA_",cancer,"_diff_meth_probes_signif_down_in_highpfi.csv"), row.names=T)
