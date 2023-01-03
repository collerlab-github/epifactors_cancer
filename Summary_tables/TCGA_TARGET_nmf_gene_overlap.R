library(dplyr)
library(reshape2)
library(ggplot2)
# cancer lists
TCGA_cancers <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                                 'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                                 'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
TARGET_cancers <- c('AML','NBL','OS','WT')
# nmf rank
rank <- 2

# initialize list of TCGA cancers and TARGET cancers
TCGA_list <- vector(mode = "list", length=length(TCGA_cancers))
names(TCGA_list) <- TCGA_cancers
TARGET_list <- vector(mode='list',length=length(TARGET_cancers))
names(TARGET_list) <- TARGET_cancers
# extract top nmf genes of each cancer
for (cancer in TCGA_cancers){
  top_list <- c(rownames(read.csv(paste0('../../TCGA_',cancer,'/03_nmf/Rank_',rank,'/nmf_lee_rank',rank,'_feature1_genes.csv'),row.names=1)),
                rownames(read.csv(paste0('../../TCGA_',cancer,'/03_nmf/Rank_',rank,'/nmf_lee_rank',rank,'_feature2_genes.csv'),row.names=1)))
  TCGA_list[[cancer]] <- top_list
  
}
for (cancer in TARGET_cancers) {
  top_list <- c(rownames(read.csv(paste0('../../TARGET_',cancer,'/03_nmf/Rank_',rank,'/nmf_lee_rank',rank,'_feature1_genes.csv'),row.names=1)),
                rownames(read.csv(paste0('../../TARGET_',cancer,'/03_nmf/Rank_',rank,'/nmf_lee_rank',rank,'_feature2_genes.csv'),row.names=1)))
  TARGET_list[[cancer]] <- top_list
  
}

# table of TARGET cancer and top nmf genes
top_genes_TARGET <- melt(TARGET_list)
colnames(top_genes_TARGET) <- c('Gene','Cancer_type')
write.csv(top_genes_TARGET, '../../Summary_tables/TARGET_nmf_top_genes_all.csv', row.names = F)
# matrix of genes x cancer indicating which gene is in top NMF genes of each cancer
top_genes_expanded <- data.frame(gene=unique(top_genes_TARGET$Gene))
rownames(top_genes_expanded) <- top_genes_expanded$gene
# place 1 where gene is in the top NMF genes of the cancer
for (cancer in TARGET_cancers){
  top_genes_expanded[,cancer] <- 0
  top_genes_expanded[TARGET_list[[cancer]],cancer] <- 1
}
# remove first col
top_genes_expanded <- top_genes_expanded[,-1]
# sort by highest frequency gene
top_genes_expanded$total <- apply(top_genes_expanded,1,sum)
top_genes_expanded <- rbind(apply(top_genes_expanded,2,sum), top_genes_expanded)
rownames(top_genes_expanded)[1] <- 'total'
top_genes_expanded <- top_genes_expanded %>% arrange(desc(total))
write.csv(top_genes_expanded, '../../Summary_tables/TARGET_nmf_top_genes_x_cancer.csv')


# overlap btwn all top NMF genes in TARGET vs TCGA
top_TARGET <- rownames(top_genes_expanded)
top_TCGA <- unique(melt(TCGA_list)[,1])
overlap_genes_total <- intersect(top_TARGET,top_TCGA)


# function to find overlapping gene count for in each TCGA cancers for one TARGET cancer
find_overlap <- function(cancer, TARGET_list, TCGA_list){
  c(length(TARGET_list[[cancer]]),sapply(TCGA_list,function(x) length(intersect(x,TARGET_list[[cancer]])), simplify=T))
  }
gene_overlap_df <- sapply(names(TARGET_list), find_overlap, TARGET_list, TCGA_list)
rownames(gene_overlap_df)[1] <- 'None' # overlap with no TCGA cancer

# use melt function from reshape2 package to convert dataframe into long form (compatible with ggplot2)
# we do this so that one measurement is in each row, not each cell
gene_overlap_df <- melt(gene_overlap_df,varnames = c('TCGA Cancer','TARGET'), 
                        value.name = 'Number of Overlapping Genes')
# plot overlaps
ggplot(gene_overlap_df, aes(x=`TCGA Cancer`,y=`Number of Overlapping Genes`, fill=`TARGET`))+
  geom_bar(stat = 'identity',position="dodge")+
  ggtitle("Overlapping Genes Between TARGET and TCGA Cancers")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  theme(legend.title = element_text(size = 8), legend.key.size = unit(0.1,'cm'))
ggsave('../../Summary_tables/TCGA_TARGET_nmf_overlapping_genes.png',width=10,height=7)

# overlaps within each project
# genes overlapping among all TCGA
Reduce(intersect, TCGA_list)
Reduce(intersect, TARGET_list)
