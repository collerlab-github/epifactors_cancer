library(argparse)

# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("project", nargs=1, help="project name (TCGA or TARGET)")
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
# parse the arguments
args <- parser$parse_args()
print(args)
project <- args$project
cancer <- args$cancer
setwd(paste0("../../",project,"_", cancer, "/02_filtering/"))


library(readxl)
library(dplyr)
# read in log2 normalized and normalized counts
log2norm_df <- read.csv(paste0("../RNA-seq_datasets/", cancer, "_log2normalized_counts.csv"), header=T,row.names=1)
norm_df <- read.csv(paste0("../RNA-seq_datasets/",cancer,"_normalized_counts.csv"),row.names=1, header=T)

# read in gencode genename index
gencode <- read_excel("../../metadata/gencode.v22.annotation.bed.xlsx",sheet="gencode.v22.annotation.bed", range = cell_cols(4:5),col_names = F)
gencode <- as.data.frame(gencode)

# fix gene names
gencode[,1]<- gsub(';','',gencode[,1])
gencode[,2] <- gsub(';','',gencode[,2])

colnames(gencode) <- c('ENSG','HUGO')
rownames(gencode) <- gencode$ENSG

# read in epigenes list
epigenes <- read.csv("../../metadata/epigenes.csv",sep = '\t')
colnames(epigenes)[1] <- 'HUGO'

# function to filter out epigenes from normalized data frame
print('Filtering epigenes')
filter_func <- function(norm_df, outfilename){

  # add corresponding hugo name column to normalized counts data
  norm_df$HUGO <- gencode[match(rownames(norm_df), rownames(gencode)),2]
  
  # check 
  head(norm_df %>%
    select(HUGO))
  
  cancer_epigenes <- norm_df[which(norm_df$HUGO %in% epigenes$HUGO),]
  # identify ENSG's with the same HUGO ID, use the refseq ENSG
  dup_ids <- cancer_epigenes$HUGO[duplicated(cancer_epigenes$HUGO)]
  cancer_epigenes[cancer_epigenes$HUGO %in% dup_ids,] %>% select(HUGO)
  
  
  # look up genes in refseq website, keep the corresponding ENSG
  
  # TAF9 --- ENSG00000273841
  # KDM4C --- ENSG00000107077
  # SMYD3 --- ENSG00000185420
  # PAGR1 -- ENSG00000280789
  
  # keep the ids listed above
  extract_ENSG <- c('ENSG00000085231.12','ENSG00000274527.3','ENSG00000280657.1','ENSG00000280893.1')
  cancer_epigenes2 <- cancer_epigenes[!rownames(cancer_epigenes)%in%extract_ENSG,]
  
  # check duplicated ids, there should only be one occurance of each HUGO gene name now
  cancer_epigenes2[cancer_epigenes2$HUGO %in% dup_ids,] %>% select(HUGO)
  
  # set rownames to HUGO names
  rownames(cancer_epigenes2) <- cancer_epigenes2$HUGO
  # write table
  write.csv(cancer_epigenes2[,-ncol(cancer_epigenes2)],file=outfilename,quote=F)
  # return epigenes dataframe
  return(cancer_epigenes2[,-ncol(cancer_epigenes2)])
}

# filter epigenes for log2 and normalized counts
log2_outfile <- paste0('../RNA-seq_datasets/', cancer, '_epigenes_log2norm_counts.csv')
norm_outfile <- paste0('../RNA-seq_datasets/', cancer, '_epigenes_norm_counts.csv')
log2_epi <- filter_func(log2norm_df, log2_outfile)
norm_epi <- filter_func(norm_df, norm_outfile)



###################################
##### Median Counts Filtering #####
###################################
# find median of normalized var genes
print('Finding median and 75th percentile counts epigenes')
norm_meds <- sort(round(apply(norm_epi,MARGIN = 1, FUN=median),digits=3))
write.csv(norm_meds, file=paste0("../RNA-seq_datasets/", cancer, "_median_normalized_counts.csv"), quote=F)
# genes with at least 25% of tumors >= 5 counts (aka 75th percentile >=5)
counts_75th <- round(apply(norm_epi, MARGIN=1, quantile, .75),digits=3)
write.csv(counts_75th, file=paste0(cancer,"_75th_percentiles.csv"), quote=F)

# filter 75th genes >=5
print('Filtering epigenes for 75th percentile >=5 based on normalized counts')
filt_epi_75th <- norm_epi %>% filter(counts_75th>=5)
write.csv(filt_epi_75th, file=paste0('../RNA-seq_datasets/', cancer, "_75th_percentile>5_counts.csv"), quote=F)

# write log2 transformed 75th norm genes
print('Outputting log2 normalized counts of the filtered epigenes')
log2_filt_epi_75th <- log2_epi[rownames(filt_epi_75th),]
write.csv(log2_filt_epi_75th, file=paste0('../RNA-seq_datasets/', cancer, "_log2_75th_percentile>5_counts.csv"), quote=F)

# find number of genes before, after filter, avg standard deviation
library(genefilter)
df <- data.frame('epigenes'=nrow(norm_epi), 'epigenes_75th_counts>=5'=nrow(filt_epi_75th), 'filtered_avg_sd'=mean(rowSds(as.matrix(filt_epi_75th))), 'log2_filtered_avg_sd'=mean(rowSds(as.matrix(log2_filt_epi_75th))))
write.csv(df,file=paste0(cancer,'_number_genes_comparison.csv'), quote=F)

# find number of genes at each sd cutoff for log2 counts
vargenes_df <- log2_filt_epi_75th
# filter out y-chromosome genes for non sex-specific cancers
if (cancer %in% c('BRCA','PRAD','OV','UCEC','CESC','TGCT')==F){
  gencode <- read_excel("../../metadata/gencode.v22.annotation.bed.xlsx",sheet="gencode.v22.annotation.bed",col_names = F)
  gencode <- as.data.frame(gencode)[,c(1,5)]
  # fix gene names
  gencode[,2]<- gsub(';','',gencode[,2])
  # set column names
  colnames(gencode) <- c('Chromosome','HUGO')
  # filter out any sd cutoff gene on y chromosome
  non_ychrom_vargenes <- unique(gencode %>% filter(HUGO %in% rownames(vargenes_df), Chromosome != 'chrY'))$HUGO
  vargenes_df <- vargenes_df[non_ychrom_vargenes,]
}
# calculate row SD 
print('Generating SD cutoff table')
gene_sd <- rowSds(as.matrix(vargenes_df))
names(gene_sd) <- rownames(vargenes_df)
cutoffs <- data.frame(matrix(ncol=6,nrow=1))
for (i in seq(1,6)){
  cutoffs[i] <- length(gene_sd[which(gene_sd>=(0.1*i))])
  if (cutoffs[i]>=450){
    sd_cutoff <<- 0.1*i
  }
}
print(paste0('best sd cutoff: ',sd_cutoff))
colnames(cutoffs) <- c(0.1,0.2,0.3,0.4,0.5,0.6)
write.csv(cutoffs,file="epigene_sd_cutoffs.csv",row.names=F, quote=F)

# filter genes with sd cutoff (for log2 genes, and normalized genes)
vargenes <- gene_sd[which(gene_sd>=sd_cutoff)]

length(vargenes)

final_vargenes_df <- vargenes_df[which(rownames(vargenes_df)%in%names(vargenes)),]

print(paste0('number of sd>=', sd_cutoff,' genes:',dim(final_vargenes_df)[1]))


# save table
print('saving filtered dataset')
write.csv(final_vargenes_df,file=paste0("../RNA-seq_datasets/",cancer,"_epigenes_log2norm_sd",sd_cutoff,"_counts.csv"),quote=F)
norm_vargenes_df <- filt_epi_75th[rownames(final_vargenes_df),]
write.csv(norm_vargenes_df, file=paste0("../RNA-seq_datasets/",cancer,"_epigenes_norm_sd",sd_cutoff,"_counts.csv"),quote=F)


