args <- commandArgs(trailingOnly=TRUE)
if (length(args)<=1){
  stop('Usage: Rscript 01_raw_counts_dataframe.R [cancer1] [cancer2] ...\nAt least 2 cancers must be supplied (i.e.)')
} else{
  print(args)
}


dir.create('../../TCGA_PAN', showWarnings=F)
dir.create('../../TCGA_PAN/gene_predictive_clinical_prediction_analysis',showWarnings=F)
setwd('../../TCGA_PAN/gene_predictive_clinical_prediction_analysis/')

library(dplyr)
library(stringr)
# read in pfi-significant genes into one data frame

pfigene_df <- data.frame(gene=c(), cancer=c())
for (cancer in args){
  genes <- rownames(read.csv(paste0('../../TCGA_',cancer,'/05_gene_clinical_prediction_analysis/',cancer,'_significant_pval_differential_genes_PFI.csv'),row.names=1))
  df <- data.frame(gene=genes, cancer=rep(cancer,length(genes)))
  pfigene_df <- rbind(pfigene_df,df)
}

# find cancer overlap for each gene
gene_overlap <- aggregate(cancer~gene, data=pfigene_df, FUN=paste, collapse='; ')
gene_overlap$num_cancer <- sapply(gene_overlap$cancer, FUN=str_count, pattern='; ')+1

gene_overlap <- arrange(gene_overlap, num_cancer, cancer)
gene_overlap2 <- gene_overlap %>% filter(num_cancer>1)
final_gene_overlap <- aggregate(gene~cancer+num_cancer, data=gene_overlap2, FUN=paste, collapse='; ')
write.csv(final_gene_overlap, file=paste0(paste(args,collapse='_'),'_overlapping_genes.csv'),quote=F,row.names=T)

# Venn Diagram visualization
library(VennDiagram)
library(RColorBrewer)
# create list with each entry as cancer's pfigenes
x <- aggregate(gene~cancer, data=pfigene_df, FUN=list)
l <- as.list(x[,2])
names(l) <- x[,1]

# plot venn diagram
f <- futile.logger::flog.threshold(futile.logger::ERROR, name="VennDiagramLogger")
colors <- brewer.pal(length(l),"Accent")
cat_labels <- paste0(names(l),': ',lapply(l,FUN=length))
v <- venn.diagram(x=l,
             category.names=cat_labels,
             filename=paste0(paste(args,collapse="_"),'_pfi_genes_venn_diagram.png'),
             ### Output features
             imagetype='png',
             resolution=300,
             margin=0.05,
             ### Title
             main = 'PFI Significant Genes Overlap',
             main.cex=3,
             main.fontfamily='sans',
             ### Circles
             lwd = 2,
             col = 'black',
             fill=colors,
             ### Numbers
             cex = 2,
             fontface = 'bold',
             fontfamily='sans',
             ### Set names
             cat.cex=2,
             cat.fontface = 'bold',
             cat.fontfamily='sans'
             )

# find functional and complexes groups

complexes <- read.delim('../../metadata/epigenes.csv')
complexes$Function <- gsub(',', ';', complexes$Function)
complexes$Protein.complex <- gsub(',', ';', complexes$Protein.complex)

# find main protein complexes and functional types
# subset epigene complex df with top_genes and order by protein complex
complex_subset <- complexes %>% filter(HGNC.approved.symbol%in%common_genes) %>% select(HGNC.approved.symbol,Function, Protein.complex)
complex_subset_prot_order <- arrange(complex_subset, desc(Protein.complex))
true_prot.complex <- complex_subset_prot_order %>% filter(Protein.complex!='#') 
if(nrow(true_prot.complex)!=0){
  complex_subset_prot_order[1:nrow(true_prot.complex),] <- arrange(true_prot.complex,Protein.complex)
}
# find enriched functinal groups
functional_groups <- aggregate(HGNC.approved.symbol~Function, data=complex_subset,FUN=paste,collapse='; ')

# write tables
write.csv(complex_subset_prot_order, file=paste0(cancer1, '_', cancer2, '_', cancer3, '_pfi_genes_protein_complexes.csv'),quote=F,row.names=F)
write.csv(functional_groups, file=paste0(cancer1, '_', cancer2, '_', cancer3, '_pfi_genes_functional_groups.csv'),quote=F,row.names=F)



#########
# find pfi significance of overlapping genes
#########

# common genes among cancers
common_genes <- str_split(final_gene_overlap[nrow(final_gene_overlap),'gene'],'; ')[[1]]

pfi_dir_df <- data.frame(row.names=common_genes)
read_pfi <- function(cancer){
  pfi_promote <- rownames(read.csv(paste0('../../TCGA_',cancer,'/05_gene_clinical_prediction_analysis/pfi_promoting_genes_protein_complexes.csv'),row.names=1))
  pfi_inhibit <- rownames(read.csv(paste0('../../TCGA_',cancer,'/05_gene_clinical_prediction_analysis/pfi_inhibiting_genes_protein_complexes.csv'),row.names=1))
  pfi_dir <- data.frame(row.names=common_genes, direction=rep('none', length(common_genes)))
  for (gene in common_genes){
    if (gene %in% pfi_promote){
      pfi_dir[gene,'direction'] <- 'longer pfi'
    } else if (gene %in% pfi_inhibit){
      pfi_dir[gene,'direction'] <- 'shorter pfi'
    }
  }
  pfi_dir_df$cancer <<- pfi_dir[,'direction']
  colnames(pfi_dir_df)[ncol(pfi_dir_df)] <<- cancer
}

# find direction of genes' association to pfi
lapply(args, read_pfi)

write.csv(pfi_dir_df, paste0(paste(args, collapse='_'),'_pfi_direction.csv'), quote=F)



