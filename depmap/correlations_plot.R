library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)
setwd('../../depmap/')
# 
# # read in top 20 pfi-significant epigenes in each cancer
# epigenes <- read_excel('../Summary_tables/single_gene_adjpvalue_and_pfi_effect_top20.xlsx')
# 
# top20 <- sapply(epigenes$Epigene, function(x) length(which(epigenes$Epigene==x)))
# top20 <- data.frame(sort(top20[unique(names(top20))], decreasing=T))
# colnames(top20) <- 'Frequency'
# 
# write.csv(top20, '../Summary_tables/epigenes_frequency.csv', quote=F)


# read in crispr dependency scores
crispr <- read.csv('depmap_TCGA_matched_epigenes_filtered_crispr.csv')
# primary disease type and subtype for cell lines in depmap
subtypes <- read.csv('TCGA_depmap_matchings.csv')

# # filter crispr scores for top20 epigenes
# crispr20 <- filter(crispr,gene_name%in%rownames(top20))
# crispr20 <- merge(crispr20,subtypes[,c('Depmap_id','Primary_Disease_Depmap','Subtype_Disease_Depmap')], by.x='depmap_id', by.y='Depmap_id')
# crispr20$disease <- gsub(' ','-',paste(crispr20$Primary_Disease_Depmap,crispr20$Subtype_Disease_Depmap))
# # append epigene frequency among top20
# crispr20$frequency <- top20[crispr20$gene_name,]
# # median dependency score for each epigene
# crispr20_med <- crispr20 %>% group_by(gene_name) %>% summarise(median(dependency,na.rm = T))
# colnames(crispr20_med)[2] <- 'dependency'
# crispr20_med$frequency <- top20[crispr20_med$gene_name,]

# Plotting Function
cor_plot <- function(df, title, filename){
  ggplot(df,aes(frequency,dependency))+
    geom_point()+
    stat_cor(method='pearson', label.x.npc = 'center', label.y.npc = 'top', size=5)+
    geom_smooth(method='lm',formula=y~x)+
    ggtitle(title)+
    xlab('Frequency of Epigene')+
    ylab('Crispr Dependency Score')+
    theme_bw()+
    theme(plot.title = element_text(size=12))
  ggsave(filename, width=7, height=7, units='in')
}

# cor_plot(crispr20_med, 'Median Dependency Score of Top 20 PFI Epigenes for Depmap Cell Lines', 'Crispr_Median_Epigene_Freq_correlation.png')
# 
# # pfi_promoting and inhibiting epigene correlation plots for each cancer type
# cancer_list <- unique(subtypes$Study_Abbreviation_TCGA)
# # function for extracting medians of subsetted epigenes
# get_medians <- function(crispr_df, depID, epigenes){
#   df <- filter(crispr_df, depmap_id%in%depID, gene_name%in%epigenes)
#   df_med <- df %>% group_by(gene_name) %>% summarise(median(dependency, na.rm=T))
#   colnames(df_med)[2] <- 'dependency'
#   df_med$frequency <- top20[df_med$gene_name,]
#   df_med$frequency[is.na(df_med$frequency)]<-0
#   return(df_med)
# }
# for (c in cancer_list) {
#   # extract depmap id's for each cancer type
#   depID <- filter(subtypes, Study_Abbreviation_TCGA==c)[,'Depmap_id']
#   # move to next iteration if no depIDs
#   if(all(is.na(depID))==T){next}
#   # files for each effect (promoting and inhibiting)
#   pro_file <- paste0('../TCGA_',c,'/05_gene_clinical_prediction_analysis/pfi_promoting_genes_protein_complexes.csv')
#   inh_file <- paste0('../TCGA_',c,'/05_gene_clinical_prediction_analysis/pfi_inhibiting_genes_protein_complexes.csv')
#   # gather genes from each effect
#   pro_genes <- rownames(read.csv(pro_file,header=T,row.names=1))
#   inh_genes <- rownames(read.csv(inh_file,header=T,row.names=1))
#   
#   # # gather genes of each effect among the top 20 pfi epigenes in the cancer
#   # pro_genes_20 <- (epigenes %>% filter(`Cancer Type`==c, Epigene %in% pro_genes))[,'Epigene',drop=T]
#   # inh_genes_20 <- (epigenes %>% filter(`Cancer Type`==c, Epigene %in% inh_genes))[,'Epigene',drop=T]
#   
#   # extract dependency scores associated with each group of epigenes and depmapID, and then return median dependency
#   crispr_pro_med <- get_medians(crispr,depID,pro_genes)
#   crispr_inh_med <- get_medians(crispr,depID,inh_genes)
#   crispr20_pro_med <- get_medians(crispr20,depID,pro_genes)
#   crispr20_inh_med <- get_medians(crispr20,depID,inh_genes)
#   # plot
#   cor_plot(crispr_pro_med, paste0('Median Dependency Score of ',c,' PFI-promoting Epigenes for Depmap Cell Lines'), 
#                                   paste0('pfi_promoting_epigene_figures/',c,' Crispr_Median_Epigene_Freq_correlation.png'))
#   cor_plot(crispr_inh_med, paste0('Median Dependency Score of ',c,' PFI-inhibiting Epigenes for Depmap Cell Lines'), 
#            paste0('pfi_inhibiting_epigene_figures/',c,' Crispr_Median_Epigene_Freq_correlation.png'))
#   cor_plot(crispr20_pro_med, paste0('Median Dependency Score of ',c,' PFI-promoting in Top 20 Epigenes for Depmap Cell Lines'), 
#            paste0('pfi_promoting_epigene_figures/',c,' Crispr_Median_top20_Epigene_Freq_correlation.png'))
#   cor_plot(crispr20_inh_med, paste0('Median Dependency Score of ',c,' PFI-inhibiting in Top 20 Epigenes for Depmap Cell Lines'), 
#            paste0('pfi_inhibiting_epigene_figures/',c,' Crispr_Median_top20_Epigene_Freq_correlation.png'))
# }

#####################################
### All Epigenes and Cancer Table ###
#####################################
library(xlsx) # use this library in order to extract color information of excel
library(stringr)
# load workbook and initialize rows and cells in order to access sheet contents
epigenes <- loadWorkbook('../Summary_tables/single_gene_adjpvalue_and_pfi_effect_transpose.xlsx')
sheet1 <- getSheets(epigenes)[[1]]
rows <- getRows(sheet1)
cells <- getCells(rows)
# get contents of excel sheet (the pvalues)
pvals <- sapply(cells, getCellValue)
pvals <- pvals[-c(1:7)]
styles <- sapply(cells, getCellStyle)
# function to extract color
cellColor <- function(style) {
  fg  <- style$getFillForegroundXSSFColor()
  rgb <- tryCatch(fg$getRgb(), error = function(e) NULL)
  rgb <- paste(rgb, collapse = "")
  return(rgb)
}
colors <- sapply(styles,cellColor)
colors <- colors[-c(1:7)]
# transform to data frame
pval_df <- data.frame(matrix(NA,nrow=722,ncol=26))
colors_df <- data.frame(matrix(NA,nrow=722,ncol=26))
# fill data frame function
fill_pval <- function(pval_vector_name, pval_vector, pval_df, colors_vector,colors_df){
  i <- str_split(pval_vector_name,'\\.')[[1]][1]
  j <- str_split(pval_vector_name,'\\.')[[1]][2]
  # fill in dataframe (adjust row index accordingly)
  pval_df[as.numeric(i)-4,as.numeric(j)] <<- pval_vector[paste0(i,'.',j)]
  # fill in color data frame
  colors_df[as.numeric(i)-4,as.numeric(j)] <<- colors_vector[paste0(i,'.',j)]
}
s <- sapply(names(pvals),fill_pval,pval_vector=pvals,pval_df=pval_df, colors_vector=colors, colors_df=colors_df)

# format df's
columns <- as.character(pval_df[1,-1])
rows <- as.character(pval_df[-1,1])
pval_df <- pval_df[-1,-1]
colors_df <- colors_df[-1,-1]
rownames(pval_df) <- rows
colnames(pval_df) <- columns
rownames(colors_df) <- rows
colnames(colors_df) <- columns
# 1 get epigenes
epigenes <- rownames(pval_df)[-1]
# 2 add TCGA cancer type to crispr df
crispr1 <- merge(crispr,subtypes[,c('Study_Abbreviation_TCGA','Cell_Line_Depmap')],by.x='cell_line',by.y='Cell_Line_Depmap')
# 3 add pvalue to crispr df
crispr1$adj_pval <- sapply(1:nrow(crispr1),function(x) pval_df[crispr1$gene_name[x],crispr1$Study_Abbreviation_TCGA[x]])
# filter out entries with no dependency score or adj_pval
crispr2 <- crispr1 %>% filter(adj_pval!='NA',is.na(dependency)==F)
crispr2 <- crispr2[,c('cell_line','gene_name','Study_Abbreviation_TCGA','adj_pval','dependency')]
# add color
crispr2$effect <- sapply(1:nrow(crispr2),function(x) ifelse(colors_df[crispr2$gene_name[x],crispr2$Study_Abbreviation_TCGA[x]]=='90ee90','pfi-promoting','pfi-inhibiting'))
# obtain frequencies
epi_frequencies <- sapply(unique(crispr2[,c('gene_name','Study_Abbreviation_TCGA')]$gene_name), 
                          function(x) length(which(unique(crispr2[,c('gene_name','Study_Abbreviation_TCGA')])$gene_name==x)))
crispr2$frequency <- epi_frequencies[crispr2$gene_name]
head(crispr2)
# write to csv and xlsx
write.csv(crispr2,'depmap_crispr_dependency_epigene_pfi_adj_pval_summary.csv',quote=F,row.names=F)

# corr plot for each pfi-effect
crispr2$adj_pval <- as.numeric(crispr2$adj_pval)
crispr2_promote <- crispr2%>%filter(effect=='pfi-promoting')
crispr2_inhibit <- crispr2%>%filter(effect=='pfi-inhibiting')
cor_plot(crispr2_promote,title = 'Dependency Scores of PFI-promoting Epigenes in Depmap Cell Lines',
         filename = 'crispr_pfi_promote_epigenes_dependency_scores_correlation.png')

cor_plot(crispr2_inhibit,title = 'Dependency Scores of PFI-inhibiting Epigenes in Depmap Cell Lines',
         filename = 'crispr_pfi_inhibit_epigenes_dependency_scores_correlation.png')

########################################
### plot dependency score vs -log10 padj value
########################################
ggplot(crispr2_promote,aes(dependency,-log10(adj_pval)))+
  geom_point()+
  ggtitle('Dependency Scores and Padj Values of PFI-promoting Epigenes in Depmap Cell Lines')+
  xlab('Dependency Score')+
  ylab('-log10(padj)')+
  theme_bw()+
  theme(plot.title = element_text(size=12))
ggsave('crispr_pfi_promote_epigenes_dependency_scores_padj_correlation.png', width=7, height=7, units='in')

ggplot(crispr2_inhibit,aes(dependency,-log10(adj_pval)))+
  geom_point()+
  ggtitle('Dependency Scores and Padj Values of PFI-inhibiting Epigenes in Depmap Cell Lines')+
  xlab('Dependency Score')+
  ylab('-log10(padj)')+
  theme_bw()+
  theme(plot.title = element_text(size=12))
ggsave('crispr_pfi_inhibit_epigenes_dependency_scores_padj_correlation.png', width=7, height=7, units='in')

########################################
### plot dependency score vs -log10 padj value with dependency cutoff
########################################

# dependency < -0.5
crispr0.5_promote <- crispr2%>%filter(dependency<(-0.5),effect=='pfi-promoting')
crispr0.5_inhibit <- crispr2%>%filter(dependency<(-0.5),effect=='pfi-inhibiting')
ggplot(crispr0.5_promote,aes(-(dependency),-log10(adj_pval)))+
  geom_point()+
  ggtitle('Dependency Scores and Padj Values of PFI-promoting Epigenes in Depmap Cell Lines')+
  xlab('-(Dependency Score)')+
  ylab('-log10(padj)')+
  # stat_cor(method='pearson', label.x.npc = 'center', label.y.npc = 'top', size=5)+
  # geom_smooth(method='lm',formula=y~x)+
  theme_bw()+
  theme(plot.title = element_text(size=12))
ggsave('crispr_pfi_promote_epigenes_dependency_scores_-0.5cutoff_padj_correlation.png', width=7, height=7, units='in')

ggplot(crispr0.5_inhibit,aes(-(dependency),-log10(adj_pval)))+
  geom_point()+
  ggtitle('Dependency Scores and Padj Values of PFI-inhibiting Epigenes in Depmap Cell Lines')+
  xlab('-(Dependency Score)')+
  ylab('-log10(padj)')+
  stat_cor(method='pearson', label.x.npc = 'center', label.y.npc = 'top', size=5)+
  geom_smooth(method='lm',formula=y~x)+
  theme_bw()+
  theme(plot.title = element_text(size=12))
ggsave('crispr_pfi_inhibit_epigenes_dependency_scores_-0.5cutoff_padj_correlation.png', width=7, height=7, units='in')

# dependency < -1
crispr1_promote <- crispr2%>%filter(dependency<(-1),effect=='pfi-promoting')
crispr1_inhibit <- crispr2%>%filter(dependency<(-1),effect=='pfi-inhibiting')
ggplot(crispr1_promote,aes(-(dependency),-log10(adj_pval)))+
  geom_point()+
  ggtitle('Dependency Scores and Padj Values of PFI-promoting Epigenes in Depmap Cell Lines')+
  xlab('-(Dependency Score)')+
  ylab('-log10(padj)')+
  # stat_cor(method='pearson', label.x.npc = 'center', label.y.npc = 'top', size=5)+
  # geom_smooth(method='lm',formula=y~x)+
  theme_bw()+
  theme(plot.title = element_text(size=12))
ggsave('crispr_pfi_promote_epigenes_dependency_scores_-1cutoff_padj_correlation.png', width=7, height=7, units='in')

ggplot(crispr1_inhibit,aes(-(dependency),-log10(adj_pval)))+
  geom_point()+
  ggtitle('Dependency Scores and Padj Values of PFI-inhibiting Epigenes in Depmap Cell Lines')+
  xlab('-(Dependency Score)')+
  ylab('-log10(padj)')+
  stat_cor(method='pearson', label.x.npc = 'center', label.y.npc = 'top', size=5)+
  geom_smooth(method='lm',formula=y~x)+
  theme_bw()+
  theme(plot.title = element_text(size=12))
ggsave('crispr_pfi_inhibit_epigenes_dependency_scores_-1cutoff_padj_correlation.png', width=7, height=7, units='in')

########################################
### plot median dependency score across cell type and epigene vs -log10 padj value with dependency cutoff 
########################################

# -0.5  cutoff
crispr0.5_medians_promote <- crispr0.5_promote %>% 
  group_by(gene_name,Study_Abbreviation_TCGA) %>%
  summarise(med_dep=median(dependency),adj_pval=unique(adj_pval),.groups = "keep")
crispr0.5_medians_inhibit <- crispr0.5_inhibit %>% 
  group_by(gene_name,Study_Abbreviation_TCGA) %>%
  summarise(med_dep=median(dependency),adj_pval=unique(adj_pval),.groups= "keep")

ggplot(crispr0.5_medians_promote,aes(-(med_dep),-log10(adj_pval)))+
  geom_point()+
  ggtitle('Dependency Scores and Padj Values of PFI-promoting Epigenes in Depmap Cell Lines')+
  xlab('-(Median Dependency Score)')+
  ylab('-log10(padj)')+
  stat_cor(method='pearson', label.x.npc = 'center', label.y.npc = 'top', size=5)+
  geom_smooth(method='lm',formula=y~x)+
  theme_bw()+
  theme(plot.title = element_text(size=12))
ggsave('crispr_pfi_promote_epigenes_median_dependency_scores_-0.5cutoff_padj_correlation.png', width=7, height=7, units='in')

ggplot(crispr0.5_medians_inhibit,aes(-(med_dep),-log10(adj_pval)))+
  geom_point()+
  ggtitle('Dependency Scores and Padj Values of PFI-inhibiting Epigenes in Depmap Cell Lines')+
  xlab('-(Median Dependency Score)')+
  ylab('-log10(padj)')+
  stat_cor(method='pearson', label.x.npc = 'center', label.y.npc = 'top', size=5)+
  geom_smooth(method='lm',formula=y~x)+
  theme_bw()+
  theme(plot.title = element_text(size=12))
ggsave('crispr_pfi_inhibit_epigenes_median_dependency_scores_-0.5cutoff_padj_correlation.png', width=7, height=7, units='in')

# -1  cutoff
crispr1_medians_promote <- crispr1_promote %>% 
  group_by(gene_name,Study_Abbreviation_TCGA) %>%
  summarise(med_dep=median(dependency),adj_pval=unique(adj_pval),.groups = "keep")
crispr1_medians_inhibit <- crispr1_inhibit %>% 
  group_by(gene_name,Study_Abbreviation_TCGA) %>%
  summarise(med_dep=median(dependency),adj_pval=unique(adj_pval),.groups= "keep")

ggplot(crispr1_medians_promote,aes(-(med_dep),-log10(adj_pval)))+
  geom_point()+
  ggtitle('Dependency Scores and Padj Values of PFI-promoting Epigenes in Depmap Cell Lines')+
  xlab('-(Median Dependency Score)')+
  ylab('-log10(padj)')+
  stat_cor(method='pearson', label.x.npc = 'center', label.y.npc = 'top', size=5)+
  geom_smooth(method='lm',formula=y~x)+
  theme_bw()+
  theme(plot.title = element_text(size=12))
ggsave('crispr_pfi_promote_epigenes_median_dependency_scores_-1cutoff_padj_correlation.png', width=7, height=7, units='in')

ggplot(crispr1_medians_inhibit,aes(-(med_dep),-log10(adj_pval)))+
  geom_point()+
  ggtitle('Dependency Scores and Padj Values of PFI-inhibiting Epigenes in Depmap Cell Lines')+
  xlab('-(Median Dependency Score)')+
  ylab('-log10(padj)')+
  stat_cor(method='pearson', label.x.npc = 'center', label.y.npc = 'top', size=5)+
  geom_smooth(method='lm',formula=y~x)+
  theme_bw()+
  theme(plot.title = element_text(size=12))
ggsave('crispr_pfi_inhibit_epigenes_median_dependency_scores_-1cutoff_padj_correlation.png', width=7, height=7, units='in')
