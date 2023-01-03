setwd("~/Desktop/SCAD_Dseq/MEGENA")

# create object datExpr0 from .Rdata file
load("~/Desktop/SCAD_Dseq/MEGENA/datExpr0.RData")


# get row means
gene_means <- rowMeans(datExpr0)
summary(gene_means)
# filter out half of genes (using median)
library(dplyr)
datExpr1 <- filter(.data=datExpr0, gene_means>13)
summary(gene_means[rownames(datExpr1)])
datExpr1[1:6,1:6]
save(datExpr1, file="datExpr1.RData")

