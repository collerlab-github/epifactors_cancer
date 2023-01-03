setwd('../../Summary_tables/')

library(xlsx) # use this library in order to extract color information of excel
library(stringr)
library(dplyr)
library(readxl)

# read in crispr dependency scores
crispr <- read.csv('../depmap/depmap_TCGA_matched_epigenes_filtered_crispr.csv')
# primary disease type and subtype for cell lines in depmap
subtypes <- read.csv('../depmap/TCGA_depmap_matchings.csv')

# load workbook and initialize rows and cells in order to access sheet contents
epigenes <- loadWorkbook('single_gene_adjpvalue_and_pfi_effect_transpose.xlsx')
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
pval_df <- data.frame(matrix(NA,nrow=720*24,ncol=4))
# fill data frame function
fill_pval <- function(pval_vector_name, pval_vector, pval_df, colors_vector){
  # print(pval_vector_name)
  i <- str_split(pval_vector_name,'\\.')[[1]][1]
  j <- str_split(pval_vector_name,'\\.')[[1]][2]
  # fill in dataframe (adjust row index accordingly)
  if (i %in% c('5')){
    if (!(j %in% c('1','26'))){
      pval_df[(0:719*24)+as.numeric(j)-1,as.numeric(i)-4] <<- pval_vector[paste0(i,'.',j)]
    }
  }
  else if (i!='6'){
    if (j =='1'){
      index <- (24*(as.numeric(i)-7))+1
      pval_df[index:(index+23),2] <<- pval_vector[paste0(i,'.',j)]
    }
    else if (j!='26'){
      index <- (24*(as.numeric(i)-7))+as.numeric(j)-1
      pval_df[index,3] <<- pval_vector[paste0(i,'.',j)]
      pval_df[index,4] <<- colors_vector[paste0(i,'.',j)]
    }
  }

}
s <- sapply(names(pvals),fill_pval,pval_vector=pvals,pval_df=pval_df, colors_vector=colors)

# format df's
colnames(pval_df) <- c("Cancer_Type","Epigene","Adj_pval","PFI_effect")
pval_df[which(pval_df$PFI_effect=='90ee90'),'PFI_effect'] <- 'promoting'
pval_df[which(pval_df$PFI_effect=='ff6a6a'),'PFI_effect'] <- 'inhibiting'
pval_df[which(pval_df$PFI_effect==''),'PFI_effect'] <- 'NA'


# read in oncoKB data
oncokb <- read.delim('../metadata/oncokb_cancerGeneList.tsv')
pval_df1 <- merge(pval_df,oncokb[,c('Hugo.Symbol', 'Is.Oncogene', 'Is.Tumor.Suppressor.Gene')],by.x='Epigene',by.y='Hugo.Symbol',all.x = T)
oncokb_combos <- paste(pval_df1$Is.Oncogene,pval_df1$Is.Tumor.Suppressor.Gene,sep=';')
oncokb_combos_table <- c('NA;NA'=NA,'No;No'='None','No;Yes'='TSG','Yes;No'='oncogene','Yes;Yes'='both')
pval_df1$Oncokb_summary <- oncokb_combos_table[oncokb_combos]
# read in cosmic data
cosmic <- read.csv('../metadata/COSMIC_cancer_gene_census.csv')
pval_df1 <- merge(pval_df1, cosmic[,c('Gene.Symbol','Role.in.Cancer')],by.x='Epigene',by.y='Gene.Symbol',all.x=T)
colnames(pval_df1)[ncol(pval_df1)] <- 'COSMIC_role_in_cancer'
pval_df1$COSMIC_role_in_cancer <- gsub(', ', '_',pval_df1$COSMIC_role_in_cancer)
# oncogene summary column
pval_df1$'Oncokb;COSMIC' <- paste(pval_df1$Oncokb_summary,pval_df1$COSMIC_role_in_cancer,sep=';')
# read in immunotherapy essential genes
immuno <- read_excel('../metadata/nature23477-s3/SI Tab 1.xlsx')
immuno1 <- na.omit(immuno[-c(1:4),2])
immuno2 <- na.omit(immuno[-c(1:4),11])
immuno3 <- na.omit(immuno[-c(1:4),14])
full_immuno <- merge(immuno1,immuno2,by=1,all=T)
full_immuno <- merge(full_immuno,immuno3,by=1,all=T)

pval_df1$immunotherapy_essential <- ifelse(pval_df1$Epigene%in%full_immuno[,1],'yes','no')

write.csv(pval_df1,'single_gene_adjpvalue_and_pfi_effect_oncokb_COSMIC_immunotherapy.csv',quote=F,row.names=F)
