library(circlize)
library(reshape2)
library(dplyr)
library(RColorBrewer)

setwd('../../prognostic_ML_classifier/')

top_nmf <- read.csv('../Summary_tables/nmf_top_genes_all.csv')
model_features <- read.csv('top_feature_SV_RF_XGB.csv')


filtered_nmf <- top_nmf %>% filter(Gene %in% model_features$Feature) %>% select(Cancer_Type, Gene)

# specify sector colors and ordering
grid.col = c(rep('blue',length(unique(filtered_nmf$Cancer_Type))),rep('grey',length(unique(filtered_nmf$Gene))))
names(grid.col) <- c(unique(filtered_nmf$Cancer_Type),unique(filtered_nmf$Gene))
grid.col[c('ACC','KIRC','LGG','LIHC','LUAD')] <- 'red'

color_order <- c("red","blue","grey")
order <- order(match(grid.col,color_order))
grid.col <- grid.col[order]
border_df <- filtered_nmf %>% filter(Cancer_Type%in%c('ACC','KIRC','LGG','LIHC','LUAD')) %>% select(Cancer_Type,Gene)
border_df$color <- 'black'

png('top_ML_features_nmf_circos.png', width=1000, height=800)
circos.par(canvas.ylim=c(-1.25,1.5), # edit  canvas size 
           track.margin = c(0.01, 0))
chordDiagram(x = filtered_nmf, grid.col=grid.col, order=names(grid.col),annotationTrack = c("grid")) #,link.border = border_df)

# par(mar = c(1, 1, 1, 1))
# circos.par("canvas.xlim" = c(0,0.5),
#            "canvas.ylim" = c(0,0.5))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE,adj=c(-0.1,0.5),cex = 1.75)
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()
dev.off()



################
## SHAP circos
###############
model_features <- read.csv('batch_corrected/top_shap_feature_SV.csv')[c(1:20),]
filtered_nmf <- top_nmf %>% filter(Gene %in% model_features$col_name) %>% select(Cancer_Type, Gene)
# specify sector colors and ordering
grid.col = c(rep('purple',length(unique(filtered_nmf$Cancer_Type))),rep('grey',length(unique(filtered_nmf$Gene))))
names(grid.col) <- c(unique(filtered_nmf$Cancer_Type),model_features$col_name)
grid.col[c('ACC','KIRC','LGG','LIHC','LUAD')] <- 'orange'
color_order <- c("orange","purple","grey")
order <- order(match(grid.col,color_order))
grid.col <- grid.col[order]
grid.col[which(grid.col=='grey')] <- colorRampPalette(c('darkgreen','lightgrey'))(length(which(grid.col=='grey')))

border_df <- filtered_nmf %>% filter(Cancer_Type%in%c('ACC','KIRC','LGG','LIHC','LUAD')) %>% select(Cancer_Type,Gene)
border_df$color <- 'black'

png('batch_corrected/top_shap_SVM_ML_features_nmf_circos.png', width=1000, height=800)
circos.par(canvas.ylim=c(-1.25,1.5), # edit  canvas size 
           track.margin = c(0.01, 0))
chordDiagram(x = filtered_nmf, grid.col=grid.col, order=names(grid.col),annotationTrack = c("grid"),link.border = border_df)

# par(mar = c(1, 1, 1, 1))
# circos.par("canvas.xlim" = c(0,0.5),
#            "canvas.ylim" = c(0,0.5))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE,adj=c(-0.1,0.5),cex = 1.75)
}, bg.border = NA) # here set bg.border to NA is important
circos.clear()
dev.off()
