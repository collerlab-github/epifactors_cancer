library(dplyr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggsci)

metapcna <- F # T or F
if(metapcna){
  setwd('../../Summary_tables/metapcna')
} else{
  setwd('../../Summary_tables/')
}
# read in top genes
read_genes <- function(project, end){
  if(metapcna){
    project <- paste0(project,'_metapcna')
  }
  analysis_genes <- read.csv(paste0(project,'_single_gene_',end,'_adjpvalues_cox_cluster.csv'))
  # columns are Gene, Cancer Type, and pvalue
  analysis_genes <- melt(analysis_genes, id.vars = 'X', variable.name = 'Gene', value.name = 'pval')
  analysis_genes <- analysis_genes[is.na(analysis_genes$pval)==F, c(2,1,3)]
  colnames(analysis_genes)[2] <- 'Cancer_Type'
  rownames(analysis_genes) <- NULL
  return(analysis_genes)
}

# read and combined prognostic genes from both analyses
TCGA_genes <- read_genes('TCGA','PFI')
TARGET_genes <- read_genes('TARGET','OS')
cancers <- c(unique(TCGA_genes$Cancer_Type),unique(TARGET_genes$Cancer_Type))
combined_genes <- rbind(TCGA_genes,TARGET_genes)
# count number of cancers each gene is prognostic in both cancer types
TCGA_gene_count <- table(TCGA_genes$Gene)
TARGET_gene_count <- table(TARGET_genes$Gene)
combined_gene_count <- table(combined_genes$Gene)
combined_genes$Number_Cancer_Types_TCGA <- TCGA_gene_count[combined_genes$Gene]
combined_genes$Number_Cancer_Types_TARGET <- TARGET_gene_count[combined_genes$Gene]
combined_genes$Number_Cancer_Types_Both <- combined_gene_count[combined_genes$Gene]


dim(combined_genes %>% filter(Number_Cancer_Types_TARGET>=2))
# filter for genes common in both cancers
combined_genes <-  combined_genes %>% filter(Number_Cancer_Types_TCGA>0, Number_Cancer_Types_TARGET>=2)

# if(project=='TCGA'){
#   common_genes <- analysis_genes%>%filter(Number_Cancer_Types>=7)
# } else{
#   common_genes <- analysis_genes%>%filter(Number_Cancer_Types>=2)
# }

# compile matrix of gene frequency among cancers
gene_matrix <- matrix(0,nrow=length(unique(combined_genes$Gene)), ncol=length(cancers))
rownames(gene_matrix) <- unique(combined_genes$Gene)
colnames(gene_matrix) <- cancers
for(i in 1:nrow(combined_genes)){
  gene_matrix[as.character(combined_genes[i,"Gene"]),combined_genes[i,"Cancer_Type"]] <- 1
}
rownames(gene_matrix) <- paste0(rownames(gene_matrix),' (',TCGA_gene_count[rownames(gene_matrix)],':',TARGET_gene_count[rownames(gene_matrix)],')')

# color annotations
col_fun = pal_npg("nrc",alpha=0.7)(9)[c(6,8)] # heatmap
annot_col <- pal_simpsons("springfield",alpha=0.7)(2) # annotation
names(annot_col) <- c('TCGA','TARGET')
ha = HeatmapAnnotation(' ' = c(rep('TCGA',length(unique(TCGA_genes$Cancer_Type))),
                                   rep('TARGET',length(unique(TARGET_genes$Cancer_Type)))),
                       col = list(' '=annot_col), 
                       annotation_legend_param = list(labels_gp=gpar(fontsize=15,fontfamily='Helvetica'),
                                                      border='black',
                                                      title_gp=gpar(fontsize=15,fontfamily='Helvetica'),
                                                      legend_height=unit(2.5,'in')),
                       show_annotation_name = F) # column annotation to identify cancer project
if(metapcna){
  prefix <- "TCGA_TARGET_metapcna_single_gene_overlap_heatmap"
  title <- 'Metapcna Prognostic Genes in Cancer Types'
} else{
  prefix<- "TCGA_TARGET_single_gene_overlap_heatmap"
  title <- 'Prognostic Genes in Cancer Types'
}

png(paste0(prefix,'.png'), units="in", width=15, height=10, res=500)
h <- Heatmap(gene_matrix,
             name = ' ',
             heatmap_legend_param = list(labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                         labels=c('Included','Not Included'), border='black',
                                         title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                         legend_height = unit(2.5, "in")),
             column_title = 'Cancer Type',
             column_title_gp = gpar(fontsize=20,fontfamily='Helvetica'),
             column_title_side = 'bottom',
             row_title = 'Genes',
             row_title_gp = gpar(fontsize=20,fontfamily='Helvetica'),
             row_names_gp = gpar(fontsize=15,fontfamily='Helvetica'),
             column_names_gp = gpar(fontsize=15,fontfamily='Helvetica'),
             column_names_rot = 90,
             bottom_annotation = ha, # column annotation
             col=col_fun)
# add more titles to heatmap
ht <- draw(h,
           column_title=paste0(title), column_title_gp=gpar(fontsize=25, fontfamily='Helvetica')
)
ht
dev.off()
write.csv(gene_matrix[row_order(ht),column_order(ht)], paste0(prefix,".csv"))


