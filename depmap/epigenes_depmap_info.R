setwd('../../depmap/')

library(depmap)
library(ExperimentHub)
library(dplyr)

# read TCGA-Depmap info table
depmap <- read.csv('TCGA_depmap_matchings.csv')
# read all epigene list
epigenes <- read.delim('../metadata/epigene_list.txt', header=F)


# get all depmap subtypes among TCGA cancers of study
depmap_subtypes <- depmap$Subtype_Disease_Depmap
# entries that have multiple subtypes separted by ';' need to be split and treated as individual entries
depmap_subtypes <- unlist(strsplit(depmap_subtypes,';'))

depmap_subtypes <- depmap_subtypes[-c(which(depmap_subtypes %in% c(NA,'none')))]


###########
# initialize experimental hub object to access depmap data
###########
eh <- ExperimentHub()
query(eh, 'depmap')

###########
# Gathering dependency scores for epigenes in cell lines that group with the depmap subtypes
###########
# rnai <- eh[['EH5357']] # rnai info for 21Q1 not present

# 1. get cell line and depmap id's with the desired subtypes and cell lineage subtype
# read metadata
metadata <- eh[['EH5362']]
metadata2 <- metadata %>% filter(lineage_subtype %in% depmap$Lineage_Subtype_Depmap)
metadata2 <- metadata %>% filter(subtype_disease%in%depmap_subtypes)
# get cell lines and depmap id of the subtypes of interest
cell_lines <- metadata2%>%select(depmap_id,cell_line)


# read crispr data
crispr <- eh[['EH5358']]
# filter depmap id's
crispr2 <- crispr %>% filter(depmap_id%in%depmap$Depmap_id)
# filter epigenes
crispr2 <- crispr2 %>% filter(gene_name%in%epigenes$V1)
# write table
write.csv(crispr2, 'depmap_TCGA_matched_epigenes_filtered_crispr.csv', row.names=F, quote=F)


# read copy number
cnv <- eh[['EH5359']]
cnv2 <- cnv %>% filter(depmap_id%in%depmap$Depmap_id)
cnv2 <- cnv2 %>% filter(gene_name%in%epigenes$V1)
write.csv(cnv2, 'depmap_TCGA_matched_epigenes_filtered_cnv.csv',row.names=F, quote=F)
