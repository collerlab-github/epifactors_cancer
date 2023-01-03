library(argparse)

# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("project", nargs=1, help='project name (e.g. TCGA or TARGET)')
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
# parse the arguments
args <- parser$parse_args()
print(args)
project <- args$project
cancer <- args$cancer


setwd(paste0("../../",project,"_", cancer, "/06_diff_gene_info/"))


if(project=='TCGA'){
  event<-'PFI'
  # 
  # cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
  #                  'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
  #                  'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
}else if(project=='TARGET'){
  event <- 'OS'
  # cancer_list <- c('AML','NBL','OS','WT')
}

library(dplyr)

# read in mutated genes brca table from cbioportal
survival_patients <- read.csv('../05_gene_clinical_prediction_analysis/survival_patients.csv')
print('Read in mutation data')
mutgenes <- read.delim('Mutated_Genes.txt', header=T, row.names=NULL, sep='\t')
colnames(mutgenes) <- c('Gene','MutSig_Q_value', 'Total_Num_Mutations', 'Num_Mutant_Patients', 
                        'Total_Num_Patients', 'Freq_of_Mutant_Patients', 'Is_Cancer_Gene_(source:OncoKB)')
mutgenes$Total_Num_Patients <- nrow(survival_patients)
mutgenes$Freq_of_Mutant_Patients <- mutgenes$Num_Mutant_Patients / mutgenes$Total_Num_Patients
# read in cna (copy number aberation) brca table from cbioportal
cnagenes <- read.delim('CNA_Genes.txt', header=T, row.names=NULL, sep='\t')
colnames(cnagenes) <- c('Gene','Gistic_Q_value', 'Cytoband', 'CNA','Num_CNA_Patients', 'Total_Num_Patients', 'Freq_of_CNA_Patients', 'Is_Cancer_Gene_(source:OncoKB)')
cnagenes$Total_Num_Patients <-  nrow(survival_patients)
cnagenes$Freq_of_Mutant_Patients <- cnagenes$Num_CNA_Patients / cnagenes$Total_Num_Patients



# # add gene discription from bioMart to mutated and cna gene tables
biomart <- read.csv('../../metadata/epigenes_biomart.csv', header=T) %>% select(Gene.description,Gene.name)
biomart <- unique(biomart)
biomart$Gene.description <- gsub(' \\[.+\\]','',biomart$Gene.description)
biomart$Gene.description <- gsub(',', ';', biomart$Gene.description)

mutgenes <- merge(x=mutgenes,y=biomart, by.x='Gene', by.y='Gene.name')
cnagenes <- merge(x=cnagenes, y=biomart, by.x='Gene', by.y='Gene.name')


# read in event significant genes
pfi <- read.csv(paste0('../05_gene_clinical_prediction_analysis/',cancer,'_significant_pval_differential_genes_',event,'.csv'), header=T, row.names=1)


# extract event genes from mutated and cna genes
mutpfi <- mutgenes %>% filter(Gene%in%rownames(pfi)) %>% arrange(desc(Freq_of_Mutant_Patients), Gene)
write.csv(mutpfi, file=paste0(event,'_mutated_genes.csv'), quote=F, row.names=F)

cnapfi <- cnagenes %>% filter(Gene%in%rownames(pfi)) %>% arrange(desc(Freq_of_CNA_Patients), Gene)
write.csv(cnapfi, file=paste0(event,'_cna_genes.csv'), quote=F, row.names=F)

