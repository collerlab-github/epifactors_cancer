setwd('../../Summary_tables/')
library(xlsx)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

project <- 'TCGA'
end <- 'PFI'

epigenes_meta <- read.csv('../metadata/epigenes.csv', sep='\t', row.names = 1)
marks <- names(table(epigenes_meta$Function))
# classify epigenetic marks
dna_mod <- marks[grepl('DNA modification', marks)]
hist_mod <- marks[grepl('Histone modification', marks)]
chrom_mod <- marks[grepl('Chromatin', marks)]
chrom_mod <- chrom_mod[!chrom_mod%in%c(dna_mod, hist_mod)]
chrom_mod <- chrom_mod[!grepl('RNA', chrom_mod)]
# Note: one entry is both Histone Modification erase and DNA modification cofactor. Which mark should it be?
epigenes_meta[epigenes_meta$Function=="Histone modification erase, DNA modification cofactor",]

epigenes_meta['Epigenetic.Mark'] <- 'None'
epigenes_meta[epigenes_meta$Function %in% dna_mod, 'Epigenetic.Mark'] <- 'DNA modification'
epigenes_meta[epigenes_meta$Function %in% hist_mod, 'Epigenetic.Mark'] <- 'Histone modification'
epigenes_meta[epigenes_meta$Function %in% chrom_mod, 'Epigenetic.Mark'] <- 'Chromatin Remodelling'

dna_mod_genes <- rownames(epigenes_meta[epigenes_meta$Epigenetic.Mark=='DNA modification',])
hist_mod_genes <- rownames(epigenes_meta[epigenes_meta$Epigenetic.Mark=='Histone modification',])
chrom_mod_genes <- rownames(epigenes_meta[epigenes_meta$Epigenetic.Mark=='Chromatin Remodelling',])

# Subset Prognostic Genes by epigenetic mark
effect_df <- read.csv(paste0(project,'_single_gene_effect_',end,'.csv'), row.names=1)
clust_df <- read.csv(paste0(project,'_single_gene_',end,'_adjpvalues_cox_cluster.csv'), row.names=1)

dna_clustdf <- clust_df[,colnames(clust_df)%in%dna_mod_genes]
dna_effectdf <- effect_df[,colnames(effect_df)%in%dna_mod_genes]
hist_clustdf <- clust_df[,colnames(clust_df)%in%hist_mod_genes]
hist_effectdf <- effect_df[,colnames(effect_df)%in%hist_mod_genes]
chrom_clustdf <- clust_df[,colnames(clust_df)%in%chrom_mod_genes]
chrom_effectdf <- effect_df[,colnames(effect_df)%in%chrom_mod_genes]



#######
# Write Excel Sheets
#######
# excel format
# create workbook
wb <- createWorkbook(type='xlsx')

# initialize a cell style. DEFAULT SETTING: CellStyle(wb, dataFormat = NULL, alignment=NULL, border=NULL, fill=NULL, font=NULL)
# Row and Columns style
ROWCOL_STYLE <- CellStyle(wb) + Font(wb, isBold = F)

# set style color for each cell according to the pfi effect
cs_promote <- CellStyle(wb)+Fill(foregroundColor='lightgreen',backgroundColor='lightgreen')
cs_inhibit <- CellStyle(wb)+Fill(foregroundColor='indianred1',backgroundColor='indianred1')
library(stringr)
add_excel_sheet <- function(wb, pval_df, effect_df, val_type){
  if(sum(is.na(pval_df))==(nrow(pval_df)*ncol(pval_df))){
    return('df contains all NA values')
  }
  # create new sheet to write
  sheet <- createSheet(wb, sheetName = paste0(end,"-Sig Epigene ",val_type))
  
  # add data frame
  # pval_df[is.na(pval_df)] <- 0
  addDataFrame(pval_df, sheet=sheet, colnamesStyle = ROWCOL_STYLE, rownamesStyle = ROWCOL_STYLE)
  rows <- getRows(sheet=sheet) # get all rows
  cells <- getCells(rows) # get all non empty cells
  values <- lapply(cells, getCellValue) # extract values
  
  # function to color excel sheet according to the effect of the gene in the cancer
  endpoint_effect_style <- function(value_name, effect_df){
    # print(values[[value_name]])
    # get index in comparison to effect_df (1 lower than the excel index)
    row_ind <- as.numeric(str_split(value_name,'\\.')[[1]][1])-1
    col_ind <- as.numeric(str_split(value_name,'\\.')[[1]][2])-1
    if(row_ind==0 || col_ind==0){
      return()
    }
    eff <- effect_df[row_ind,col_ind]
    na_val <- is.na(values[[value_name]]) | (values[[value_name]]=='')
    if (eff==1 && !na_val){
      setCellStyle(cell=cells[[value_name]],cellStyle = cs_promote)
    } else if (eff==-1 && !na_val) {
      setCellStyle(cell=cells[[value_name]],cellStyle = cs_inhibit)
    } else{
      setCellValue(cell=cells[[value_name]],value = 'NA')
    }
  }
  color_coding <- sapply(names(values),endpoint_effect_style, effect_df)
  
  # # save workbook to file
  # saveWorkbook(wb, filename)
  
}

add_excel_sheet(wb, dna_clustdf, dna_effectdf, 'DNA Modification')
add_excel_sheet(wb, hist_clustdf, hist_effectdf, 'Histone Modification')
add_excel_sheet(wb, chrom_clustdf, chrom_effectdf, 'Chromatin Remodelling')
# save workbook to file
saveWorkbook(wb, paste0(project,'_single_gene_',end,'_epigenetic_marks_cox_cluster.xlsx'))

# add single gene prognostic effect summary to epigenes metadata
count_effects <- function(gene,clust_df, effect_df){
  return(c('num_pfi_promote'=sum(effect_df[!is.na(clust_df[,gene]), gene]==1),
           'num_pfi_inhibit'=sum(effect_df[!is.na(clust_df[,gene]), gene]==-1)))
}
epigenes_meta[,c('num_cancers_pfi_promote', 'num_cancers_pfi_inhibit')]<- NA
effect_genes_count <- t(cbind(sapply(colnames(dna_clustdf), count_effects, dna_clustdf, dna_effectdf),
                        sapply(colnames(hist_clustdf), count_effects, hist_clustdf, hist_effectdf),
                        sapply(colnames(chrom_clustdf), count_effects, chrom_clustdf, chrom_effectdf)))
epigenes_meta[rownames(effect_genes_count), c('num_cancers_pfi_promote',"num_cancers_pfi_inhibit")] <- effect_genes_count

write.csv(epigenes_meta, '../metadata/epigenes_marks_prognostic_info.csv')
