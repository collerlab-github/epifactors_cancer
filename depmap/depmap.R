# BiocManager::install('depmap')

setwd('../../depmap/')
library(depmap)
library(ExperimentHub)
library(dplyr)

###########
# Identify cancer cell line 
# initialize experimental hub object to access depmap data
###########
eh <- ExperimentHub()
query(eh, 'depmap')

# match cancer cell lines to our cancer types
metadata <- eh[['EH5362']]

# explore data
head(metadata)
head(metadata)[1:9]
head(metadata)[10:17]
head(metadata)[18:25]
head(metadata)[26]
unique(metadata[,c('primary_disease','subtype_disease')])
# lineage subtype distinguishes the different cell lines
# apply(unique(metadata[,'primary_disease']),1,print)
# 
# unique(metadata[which(metadata['lineage_subtype']=='exocrine'),'lineage'])
# unique(metadata[which(metadata['primary_disease']=='Sarcoma'),'lineage_subtype'])

# curate table for TCGA cancer type and similar depmap cell line subtype lineage
TCGA_list <- c('BRCA','THCA','OV','OV','LGG','PRAD','SKCM','UCEC','UCEC','KIRC','CRC',
                 'CESC','CESC','LIHC','SARC','SARC','SARC','SARC','SARC','SARC','SARC','SARC','SARC','SARC',
               'HNSC','HNSC','KIRP','GBM','GBM','GBM','PCPG','LUAD','STAD','LUSC','ESCA','ESCA','TGCT','ACC','BLCA','PAAD')
cell_line_list <- c('breast_adenocarcinoma','thyroid_carcinoma','ovary_adenocarcinoma','ovary_carcinoma','glioma',
                    'prostate_adenocarcinoma','melanoma', 'endometrial_adenocarcinoma','endometrial_squamous','renal_cell_carcinoma',
                    'colorectal_adenocarcinoma','cervical_squamous','cervical_adenocarcinoma','hepatocellular_carcinoma',
                    'thyroid_sarcoma','rhabdomyosarcoma','fibrosarcoma','leiomyosarcoma','uterine_sarcoma','pleomorphic_sarcoma',
                    'undifferentiated_sarcoma','synovial_sarcoma','MPNST','malignant_fibrous_histiocytoma','upper_aerodigestive_squamous',
                    'upper_aerodigestive_carcinoma','none','meningioma','medulloblastoma','neuroblastoma','none',
                    'lung_carcinoma', 'gastric_adenocarcinoma', 'NSCLC', 'esophagus_adenocarcinoma','esophagus_squamous',
                    'none', 'adrenal_carcinoma', 'bladder_carcinoma', 'exocrine'
                    )

# TCGA study information
TCGA <- read.csv('../metadata/TCGA_full_study_name.csv',row.names=1)
full_name_TCGA_list <- TCGA[TCGA_list,'Study.Name']

TCGA_depmap_index <- data.frame('TCGA'=TCGA_list, 'TCGA_full_name'=full_name_TCGA_list, 'depmap'=cell_line_list)


# make table
cancer_table <- metadata[which(metadata$lineage_subtype%in%cell_line_list),c('depmap_id','cell_line','primary_disease','subtype_disease','lineage','lineage_subtype')]
# merge TCGA names with depmap cancer table
cancer_table <- merge(cancer_table,TCGA_depmap_index, by.x='lineage_subtype', by.y='depmap')
# add TCGA cancers with no lineage subtype in depmap
TCGA_na_lineage <- TCGA_depmap_index %>% filter(TCGA %in% unique(cancer_table$TCGA)==F)
na_rows <- data.frame(TCGA_na_lineage$depmap, NA, NA, NA, NA, NA, TCGA_na_lineage$TCGA, TCGA_na_lineage$TCGA_full_name)
colnames(na_rows) <- colnames(cancer_table)
cancer_table <- rbind(cancer_table, na_rows)
colnames(cancer_table) <- c('Lineage_Subtype_Depmap','Depmap_id','Cell_Line_Depmap','Primary_Disease_Depmap',
                            'Subtype_Disease_Depmap','Lineage_Depmap','Study_Abbreviation_TCGA','Study_Name_TCGA')
cancer_table <- cancer_table[,c(7,8,2,3,4,5,6,1)] %>% arrange(Study_Abbreviation_TCGA)

cancer_table$Subtype_Disease_Depmap <- gsub(', ','-',cancer_table$Subtype_Disease_Depmap)

write.csv(cancer_table,'TCGA_depmap_matchings.csv',quote=F,row.names=F)
