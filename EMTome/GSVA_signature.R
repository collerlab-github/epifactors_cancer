# install.packages("BiocManager")
BiocManager::install("GSVA")

library(GSVA)
library(readxl)

# read gene expression data
genes <- read.csv('../../TCGA_PAN/RNA-seq_datasets/BRCA_THCA_OV_LGG_PRAD_SKCM_UCEC_KIRC_CRC_CESC_LIHC_SARC_HNSC_KIRP_GBM_PCPG_LUAD_STAD_LUSC_ESCA_TGCT_ACC_BLCA_PAAD_epigenes_log2norm_sd0.5_counts.csv',
                  row.names=1)
# patients <- gsub('(.*)\\w$','\\1',colnames(epigenes))
patients <- gsub('\\.','-',colnames(genes))
rm(genes)
# epigenes <- rownames(epigenes)

# counts already log normalized
pancan <- read.table('../../TCGA_PAN/RNA-seq_datasets/PANCAN_rsem_batch_norm_counts')
colnames(pancan) <- pancan[1,]
names(patients) <- gsub('(.*)\\w$','\\1',patients)

pancan_filtered <- pancan[, which(colnames(pancan) %in% names(patients))]

colnames(pancan_filtered) <- patients[colnames(pancan_filtered)]
pancan_filtered <- na.omit(pancan_filtered)
rnames <- pancan[rownames(pancan_filtered),'sample']
pancan_filtered <- pancan_filtered[!duplicated(rnames),]
rnames <- rnames[!duplicated(rnames)]
pancan_filtered <- apply(pancan_filtered, 2, as.numeric)
rownames(pancan_filtered) <- rnames
pancan_filtered <- na.omit(pancan_filtered)


# read EMT genes
emt_genes <- read_excel('../../metadata/EMTome_genes.xlsx')
emt_genes <- emt_genes$Gene
emt_genes <- emt_genes[emt_genes%in%rownames(pancan_filtered)]

gs <- list('EMT'=emt_genes)

gsva.es <- gsva(pancan_filtered, gs, verbose=FALSE)

write.csv(as.data.frame(t(gsva.es)), '../../metadata/EMTome_EMT_scores.csv')
