
library(xlsx)
library(ComplexHeatmap)
library(circlize)

################
#####INPUTS#####
################
project <- 'TARGET' # TCGA or TARGET
if (project=='TCGA'){
  cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                   'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                   'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
  end <- 'PFI'
} else if (project=='TARGET'){
  cancer_list <- c('AML','NBL','WT')
  end <- 'OS'
} else{
  quit(save='no')
}################

setwd('../../Summary_tables/')


# epigenes
epigene_list <- read.csv('../metadata/epigenes.csv', sep='\t')[,1]

# create epigene x cancer mutation count data frame
mut_df <- as.data.frame(matrix(NA, ncol=length(epigene_list), nrow=length(cancer_list)))
rownames(mut_df) <- cancer_list
colnames(mut_df) <- epigene_list

# create same dimension data frame for CNA amplification
amp_df <- mut_df

# create same dimension data frame for CNA deletion
del_df <- mut_df

# read in each cancer's mutation and cnv data
prog_genes <- c()
for (c in cancer_list) {
  pval_file <- paste0('../',project,'_',c,'/05_gene_clinical_prediction_analysis/',c,'_significant_cox_cluster_differential_genes_',end,'.csv')
  prog_genes <- c(prog_genes, rownames(read.csv(pval_file, row.names=1)))
  # fill in mutation data frame
  # file
  cancer_mut_file <- paste0('../',project,'_',c,'/06_diff_gene_info/',end,'_mutated_genes.csv')
  # read and fill in cancer's gene mutation frequencies
  muts <- read.csv(cancer_mut_file, header=T, row.names=1)
  mut_df[c,rownames(muts)] <- muts$Num_Mutant_Patients/muts$Total_Num_Patients
  
  # fill in cna data frames
  # files for each effect (promoting and inhibiting)
  cancer_cna_file <- paste0('../',project,'_',c,'/06_diff_gene_info/',end,'_cna_genes.csv')
  # gather genes from each effect
  cna <- read.csv(cancer_cna_file, header=T)
  amp <- cna[which(cna$CNA=='AMP'),]
  rownames(amp) <- amp$Gene
  amp <- amp[,-1]
  del <- cna[which(cna$CNA=='HOMDEL'),]
  rownames(del) <- del$Gene
  del <- del[,-1]
  # update amp and del data frames
  amp_df[c,rownames(amp)] <- amp$Num_CNA_Patients/amp$Total_Num_Patients
  del_df[c,rownames(del)] <- del$Num_CNA_Patients/del$Total_Num_Patients
}

# write tables
prog_genes <- unique(prog_genes)
mut_df <- mut_df[,prog_genes]
amp_df <- amp_df[,prog_genes]
del_df <- del_df[,prog_genes]
write.csv(mut_df, paste0(project,'_single_gene_',end,'_mutation_summary.csv'), row.names=T, quote=F)
write.csv(amp_df, paste0(project,'_single_gene_',end,'_cna_amp_summary.csv'), row.names=T, quote=F)
write.csv(del_df, paste0(project,'_single_gene_',end,'_cna_del_summary.csv'), row.names=T, quote=F)

###########
# heatmaps
# NOTE: too many NA values prohibit clustering of cancers and genes
###########

# # function to label heatmap cells
# cell_function <- function(j,i,x,y,width,height,fill){
#   grid.text(label = round(pval_df[i,j],digits=7),x = x,y = y)
# }
# heatmap
png(paste0(project,'_single_gene_',end,'_mutation_no_cluster.png'),units='in',width=7,height=7,res=500)
h <- draw(Heatmap(as.matrix(mut_df), cluster_rows = F, cluster_columns = F,
        name = 'frequency',
        column_title = paste0('Single Gene ',end,'Mutation Frequency'),
        na_col='black',
        rect_gp = gpar(col='black',lwd=0.05),
        show_column_names = F
        )
)
dev.off()

###########
# heatmap of cna
###########

# amplification
png(paste0(project,'_single_gene_effect_',end,'_cna_amp_no_cluster.png'),units='in',width=7,height=7,res=500)
Heatmap(as.matrix(amp_df),
        name = 'CNA Amplification',
        column_title = paste0(end,'-significant CNA Amplification Epigenes in Cancers'),
        rect_gp = gpar(col='black',lwd=0.05),
        cluster_rows=F, cluster_columns=F,
        show_column_names = F
        )
dev.off()

# deletion
png(paste0(project,'_single_gene_effect_',end,'_cna_del_no_cluster.png'),units='in',width=7,height=7,res=500)
Heatmap(as.matrix(del_df),
        name = 'CNA Deletion',
        column_title = paste0(end,'-significant CNA Deletion Epigenes in Cancers'),
        rect_gp = gpar(col='black',lwd=0.05),
        cluster_rows=F, cluster_columns=F,
        show_column_names = F
)
dev.off()



#######
# Write Excel Sheets
#######

# read in effect xlsx table
effect_df <- read.csv(paste0(project,'_single_gene_effect_',end,'.csv'), row.names=1, header=T)
effect_df <- effect_df[rownames(mut_df), colnames(mut_df)]
effect_df[is.na(mut_df)] <- 0
# excel format
# create workbook
wb <- createWorkbook(type='xlsx')

# initialize a cell style. DEFAULT SETTING: CellStyle(wb, dataFormat = NULL, alignment=NULL, border=NULL, fill=NULL, font=NULL)
# Row and Columns style
ROWCOL_STYLE <- CellStyle(wb) + Font(wb, isBold = F)

# create new sheet to write
sheet <- createSheet(wb, sheetName = paste0(end,"-Sig Epigene Mutation Frequency"))

# add data frame
mut_df[is.na(mut_df)] <- 0
addDataFrame(mut_df, sheet=sheet, colnamesStyle = ROWCOL_STYLE, rownamesStyle = ROWCOL_STYLE)
rows <- getRows(sheet=sheet) # get all rows
cells <- getCells(rows) # get all non empty cells
values <- lapply(cells, getCellValue) # extract values

# set style color for each cell according to the pfi effect
cs_promote <- CellStyle(wb)+Fill(foregroundColor='lightgreen',backgroundColor='lightgreen')
cs_inhibit <- CellStyle(wb)+Fill(foregroundColor='indianred1',backgroundColor='indianred1')
library(stringr)
# function to color excel sheet according to the effect of the gene in the cancer
pfi_effect_style <- function(value_name){
  print(value_name)
  # get index in comparison to effect_df (1 lower than the excel index)
  row_ind <- as.numeric(str_split(value_name,'\\.')[[1]][1])-1
  col_ind <- as.numeric(str_split(value_name,'\\.')[[1]][2])-1
  if(row_ind==0 || col_ind==0){
    return()
  }
  eff <- effect_df[row_ind,col_ind]
  if (eff==1){
    setCellStyle(cell=cells[[value_name]],cellStyle = cs_promote)
  } else if (eff==-1) {
    setCellStyle(cell=cells[[value_name]],cellStyle = cs_inhibit)
  }
}
color_coding <- sapply(names(values),pfi_effect_style)

# save workbook to file
saveWorkbook(wb, paste0(project,'_single_gene_',end,'_mutation_summary.xlsx'))

###################################
# excel sheet for cna amplification
effect_df <- read.csv(paste0(project,'_single_gene_effect_',end,'.csv'), row.names=1, header=T)
effect_df <- effect_df[rownames(amp_df), colnames(amp_df)]
effect_df[is.na(amp_df)] <- 0

# create workbook
wb2 <- createWorkbook(type='xlsx')

# initialize a cell style. DEFAULT SETTING: CellStyle(wb, dataFormat = NULL, alignment=NULL, border=NULL, fill=NULL, font=NULL)
# Row and Columns style
ROWCOL_STYLE <- CellStyle(wb2) + Font(wb2, isBold = F)

# create new sheet to write
sheet <- createSheet(wb2, sheetName = paste0(end,"-Sig Epigene CNA Amplfication Frequency"))

# add data frame
amp_df[is.na(amp_df)] <- 0
addDataFrame(amp_df, sheet=sheet, colnamesStyle = ROWCOL_STYLE, rownamesStyle = ROWCOL_STYLE)
rows <- getRows(sheet=sheet) # get all rows
cells <- getCells(rows) # get all non empty cells
values <- lapply(cells, getCellValue) # extract values

# set style color for each cell according to the pfi effect
cs_promote <- CellStyle(wb2)+Fill(foregroundColor='lightgreen',backgroundColor='lightgreen')
cs_inhibit <- CellStyle(wb2)+Fill(foregroundColor='indianred1',backgroundColor='indianred1')
library(stringr)
# color excel sheet according to the effect of the gene in the cancer
color_coding <- sapply(names(values),pfi_effect_style)

# save workbook to file
saveWorkbook(wb2, paste0(project,'_single_gene_',end,'_cna_amp_summary.xlsx'))

###################################
# excel sheet for cna deletion
effect_df <- read.csv(paste0(project,'_single_gene_effect_',end,'.csv'), row.names=1, header=T)
effect_df <- effect_df[rownames(del_df), colnames(del_df)]
effect_df[is.na(del_df)] <- 0

# create workbook
wb3 <- createWorkbook(type='xlsx')

# initialize a cell style. DEFAULT SETTING: CellStyle(wb, dataFormat = NULL, alignment=NULL, border=NULL, fill=NULL, font=NULL)
# Row and Columns style
ROWCOL_STYLE <- CellStyle(wb3) + Font(wb3, isBold = F)

# create new sheet to write
sheet <- createSheet(wb3, sheetName = paste0(end,"-Sig Epigene CNA Deletion Frequency"))

# add data frame
del_df[is.na(del_df)] <- 0
addDataFrame(del_df, sheet=sheet, colnamesStyle = ROWCOL_STYLE, rownamesStyle = ROWCOL_STYLE)
rows <- getRows(sheet=sheet) # get all rows
cells <- getCells(rows) # get all non empty cells
values <- lapply(cells, getCellValue) # extract values

# set style color for each cell according to the pfi effect
cs_promote <- CellStyle(wb3)+Fill(foregroundColor='lightgreen',backgroundColor='lightgreen')
cs_inhibit <- CellStyle(wb3)+Fill(foregroundColor='indianred1',backgroundColor='indianred1')
library(stringr)
# color excel sheet according to the effect of the gene in the cancer
color_coding <- sapply(names(values),pfi_effect_style)

# save workbook to file
saveWorkbook(wb3, paste0(project,'_single_gene_',end,'_cna_del_summary.xlsx'))

