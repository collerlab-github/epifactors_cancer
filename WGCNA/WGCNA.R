



project <- 'TCGA'
cancer <- 'ACC'

setwd(paste0("../../", project, "_", cancer, "/WGCNA/"))


BiocManager::install("WGCNA")


library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

