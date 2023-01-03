library(circlize)
library(reshape2)
library(dplyr)

setwd('../../Summary_tables/')
project <- 'TCGA'
# circos plot for cancers and pfi-sig genes
cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                 'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                 'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
top_num <- 10
top <- rownames(read.csv('../Summary_tables/TCGA_single_gene_PFI_epigene_summary.csv',row.names=1)[1:top_num,])

# create adjacency list
adj_df <- data.frame(matrix(ncol=2,nrow=0))
colnames(adj_df) = c('cancer','gene')

for (cancer in cancer_list) {
  fname <- paste0('../',project,'_',cancer,'/05_gene_clinical_prediction_analysis/',cancer,'_significant_pval_differential_genes_PFI.csv')
  genes <- top[which(top %in% rownames(read.csv(fname,row.names=1)))]
  if (length(genes)>1){
    c_df <- data.frame(cancer=cancer,gene=top[which(top%in%genes)])
    adj_df <- rbind(adj_df,c_df)
  }
}
adj_df$value <- 1

# specify sector colors and ordering
grid.col = c(rep('blue',length(unique(adj_df$cancer))),rep('grey',length(unique(adj_df$gene))))
names(grid.col) <- c(unique(adj_df$cancer),unique(adj_df$gene))
grid.col[c('ACC','KIRC','LGG','LIHC','LUAD')] <- 'red'

color_order <- c("red","blue","grey")
order <- order(match(grid.col,color_order))
grid.col <- grid.col[order]
border_df <- adj_df %>% filter(cancer%in%c('ACC','KIRC','LGG','LIHC','LUAD')) %>% select(cancer,gene)
border_df$color <- 'black'
# plot chord diagram
png(paste0(project,'_single_gene_circos_plot.png'),width=1000,height=800)
circos.par(canvas.ylim=c(-1.15,1.25), # edit  canvas size 
           track.margin = c(0.01, 0))
chordDiagram(x = adj_df, grid.col=grid.col, order=names(grid.col),annotationTrack = c("grid"),link.border = border_df)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE,adj=c(-0.1,0.5),cex = 2)
}, bg.border = NA) # here set bg.border to NA is important
title("Top PFI Prognostic Genes in Cancer Types", cex.main=3, line=-2)
circos.clear()
dev.off()

