library(argparse)
library(dplyr)

# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
# parse the arguments
args <- parser$parse_args()
print(args)
cancer <- args$cancer

setwd(paste0("../..//TCGA_", cancer, "/06_diff_gene_info/"))


depmap_dependencies <- read.csv(paste0('../supplemental/depmap_dependencies_enriched_in_',cancer,'.csv'),header=T) %>%
  filter(Type=='gene')

pfi_epigene <- read.csv(paste0('../05_gene_clinical_prediction_analysis/',cancer,'_significant_pval_differential_genes_PFI.csv'))$X


depmap_epigene <- depmap_dependencies %>% filter(Gene.Compound %in% pfi_epigene)

depmap_epigene$Dataset <- gsub(',', ';', depmap_epigene$Dataset)

write.csv(depmap_epigene, 'depmap_dependencies_epigenes_cell_lines.csv', quote=F, row.names=F)

