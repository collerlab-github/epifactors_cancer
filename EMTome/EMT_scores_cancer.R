cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                 'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                 'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')


emt_scores <- read.csv('../../metadata/EMTome_EMT_scores.csv', row.names = 1)

emt_scores$cancer <- 'None'
emt_scores$NMF <- 0
df <- data.frame()
for (cancer in cancer_list){
  nmf <- read.csv(paste0('../../TCGA_',cancer,'/03_nmf/Rank_2/nmf_lee_rank2_cluster_membership.csv'), row.names = 1)
  nmf$cancer <- cancer
  df <- rbind(df,nmf)
  # print(cancer)
  # print(table(rownames(nmf) %in% rownames(emt_scores)))
  # emt_scores[rownames(nmf),"cancer"] <- cancer
  # emt_scores[rownames(nmf),"NMF"] <- nmf$cluster
}

df$emt <- emt_scores[rownames(df),'EMT']
write.csv(df, '../../metadata/EMTome_EMT_scores.csv')
