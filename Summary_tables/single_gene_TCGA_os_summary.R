setwd('../../Summary_tables/')

dir.create('TCGA_os',showWarnings=F)
setwd('TCGA_os/')

library(xlsx)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

################
#####INPUTS#####
################
project <- 'TCGA' # TCGA or TARGET
if (project=='TCGA'){
  cancer_list <- c('BRCA','THCA','OV','LGG','PRAD','SKCM','UCEC','KIRC','CRC',
                   'CESC','LIHC','SARC','HNSC','KIRP','GBM','PCPG','LUAD',
                   'STAD','LUSC','ESCA','TGCT','ACC','BLCA','PAAD')
  end <- 'OS'
} else if (project=='TARGET'){
  cancer_list <- c('AML','NBL','OS','WT')
  end <- 'OS'
} else{
  quit(save='no')
}

################

# epigenes
epigene_list <- read.csv('../../metadata/epigenes.csv', sep='\t')[,1]

# create cancer x epigene endpoint pvalue data frame
pval_df <- as.data.frame(matrix(NA, ncol=length(epigene_list), nrow=length(cancer_list)))
rownames(pval_df) <- cancer_list
colnames(pval_df) <- epigene_list

# create cancer x epigene endpoint cox pvalue data frames for cluster, age, and gender 
cluster_df <- pval_df
age_df <- pval_df
gender_df <- pval_df

# create same dimension data frame for effect of each gene (-1=endpoint-inhibiting, 1=endpoint-promoting, 0=not significant)
effect_df <- as.data.frame(matrix(0, ncol=length(epigene_list), nrow=length(cancer_list)))
rownames(effect_df) <- cancer_list
colnames(effect_df) <- epigene_list

# create a data frame with rows as cancer types, and columns 1: # sig epigenes, 2: # sig epigene with sig cox cluster, 3: # sig epigene with sig cox age, 
#                                                           4: # sig epigene with sig cox gender epigenes 5: # endpoint-promoting epigenes, 6: # endpoint-inhibiting epigenes
cancer_summary_df <- as.data.frame(matrix(0, ncol=6, nrow=length(cancer_list)))
rownames(cancer_summary_df) <- cancer_list
colnames(cancer_summary_df) <- c('Num Sig Epigenes', 'Num Sig Epigene and Cox Cluster', 
                                 'Num Sig Epigene and Cox Age', 'Num Sig Epigene and Cox Gender',
                                 paste0('Num ',end,'-promoting Sig Epigenes'), paste0('Num ',end,'-inhibiting Sig Epigenes'))

# create a data frame with rows as epigenes, and columns 1: # cancer types that have significance, 2: # cancer types that have epigene and cox cluster sig, 3: # cancer types that have epigene and cox age sig
                                                      # 4: # cancer types that have epigene and cox gender sig, 5: # cancer types with epigene as endpoint-promoting, 6: # cancer types with epigene as endpoint-inhibiting
epigene_summary_df <- as.data.frame(matrix(0, ncol=6, nrow=length(epigene_list)))
rownames(epigene_summary_df) <- epigene_list
colnames(epigene_summary_df) <- c('Num Sig Cancer Types', 'Num Sig Cancer Types with Sig Cox Cluster', 
                                  'Num Sig Cancer Types with Sig Cox Age', 'Num Sig Cancer Types with Sig Cox Gender',
                                  paste0('Num ',end,'-promoting Sig Cancer Types'), paste0('Num ',end,'-inhibiting Sig Cancer Types'))

# read in each cancer's survival curve p-values
get_sig_genes <- function(filename){
  endpoint_pvals <- read.csv(filename, header=T, row.names=1)
  endpoint_pvals <- endpoint_pvals[which(endpoint_pvals[,2]<0.05),2,drop=F]
  return(endpoint_pvals)
}

for (c in cancer_list) {
  # fill in pvalue data frame
  # files for pvalue and cox pvalues
  cancer_file <- paste0('../../',project,'_',c,'/05a_TCGA_clinical_prediction_analysis/',c,'_significant_pval_differential_genes_',end,'.csv')
  cluster_file <- paste0('../../',project,'_',c,'/05a_TCGA_clinical_prediction_analysis/',c,'_significant_cox_cluster_differential_genes_',end,'.csv')
  age_file <- paste0('../../',project,'_',c,'/05a_TCGA_clinical_prediction_analysis/',c,'_significant_cox_age_differential_genes_',end,'.csv')
  gender_file <- paste0('../../',project,'_',c,'/05a_TCGA_clinical_prediction_analysis/',c,'_significant_cox_gender_differential_genes_',end,'.csv')
  
  # read and fill in cancer significant endpoint genes; use adjusted p-value
  sig_pval_genes <- get_sig_genes(cancer_file)
  sig_pval_num <- nrow(sig_pval_genes)
  pval_df[c,rownames(sig_pval_genes)] <- sig_pval_genes$adjpval
  
  sig_cluster_genes <- get_sig_genes(cluster_file) # cox cluster pvalues
  sig_cluster_num <- nrow(sig_cluster_genes)
  cluster_df[c,rownames(sig_cluster_genes)] <- sig_cluster_genes$adjpval_cluster
  
  sig_age_genes <- get_sig_genes(age_file) # cox age pvalues
  sig_age_num <- nrow(sig_age_genes)
  age_df[c,rownames(sig_age_genes)] <- sig_age_genes$adjpval_age
  
  sig_gender_genes <- get_sig_genes(gender_file) # cox gender pvalues
  sig_gender_num <- nrow(sig_gender_genes)
  gender_df[c,rownames(sig_gender_genes)] <- sig_gender_genes$adjpval_gender
  
  
  # fill in effect data frame
  # files for each effect (promoting and inhibiting)
  pro_file <- paste0('../../',project,'_',c,'/05a_TCGA_clinical_prediction_analysis/',end,'_promoting_genes_protein_complexes.csv')
  inh_file <- paste0('../../',project,'_',c,'/05a_TCGA_clinical_prediction_analysis/',end,'_inhibiting_genes_protein_complexes.csv')
  # gather genes from each effect
  if(file.exists(pro_file)) {
    pro_genes <- rownames(read.csv(pro_file,header=T,row.names=1))
  } else pro_genes <- c()
  if(file.exists(inh_file)){
    inh_genes <- rownames(read.csv(inh_file,header=T,row.names=1))
  } else inh_genes <- c()
  # update effects_df
  effect_df[c,pro_genes] <- 1 # endpoint-promoting = 1
  effect_df[c,inh_genes] <- -1 # endpoint-inhibiting = -1
  # update cancer_summary_df
  cancer_summary_df[c,] <- c(sig_pval_num, sig_cluster_num, sig_age_num, sig_gender_num, length(pro_genes), length(inh_genes))
  # update epigene_summary_df
  epigene_summary_df[rownames(sig_pval_genes),1] <- epigene_summary_df[rownames(sig_pval_genes),1]+1
  epigene_summary_df[rownames(sig_cluster_genes),2] <- epigene_summary_df[rownames(sig_cluster_genes),2]+1
  epigene_summary_df[rownames(sig_age_genes),3] <- epigene_summary_df[rownames(sig_age_genes),3]+1
  epigene_summary_df[rownames(sig_gender_genes),4] <- epigene_summary_df[rownames(sig_gender_genes),4]+1
  epigene_summary_df[pro_genes,paste0('Num ',end,'-promoting Sig Cancer Types')] <- epigene_summary_df[pro_genes,paste0('Num ',end,'-promoting Sig Cancer Types')]+1
  epigene_summary_df[inh_genes,paste0('Num ',end,'-inhibiting Sig Cancer Types')] <- epigene_summary_df[inh_genes,paste0('Num ',end,'-inhibiting Sig Cancer Types')]+1
}

# write cancer_summary table
write.csv(cancer_summary_df[order(cancer_summary_df$`Num Sig Epigenes`, decreasing = T),], paste0(project,'_single_gene_',end,'_cancer_summary.csv'), row.names=T, quote=F)
write.csv(epigene_summary_df[order(epigene_summary_df$`Num Sig Cancer Types`, decreasing = T),], paste0(project,'_single_gene_',end,'_epigene_summary.csv'), row.names=T, quote=F)

###########
# heatmaps of p-values
# NOTE: too many NA values prohibit clustering of cancers and genes
###########

# # function to label heatmap cells
# cell_function <- function(j,i,x,y,width,height,fill){
#   grid.text(label = round(pval_df[i,j],digits=7),x = x,y = y)
# }
color_function <- colorRamp2(c(0,0.05,0.06,1),c('blue','red','black','black'))
# heatmap with clustering, NA values are imputed as 1 to allow clustering to work
pval_heatmap <- function(pval_df, filename, val_type){
  if(sum(is.na(pval_df))==(nrow(pval_df)*ncol(pval_df))){
    return('df contains all NA values')
  }
  pval_df2 <- pval_df
  pval_df2[is.na(pval_df2)] <- 1
  png(filename,units='in',width=7,height=7,res=500)
  h <- draw(Heatmap(as.matrix(pval_df2), cluster_rows = T, cluster_columns = T,
          name = 'adj p-value',
          heatmap_legend_param = list(title='p-value',at=c(0,0.05,1),labels=c('0','0.05','NA')),
          column_title = paste0('Single Gene ',end,' ', val_type),
          rect_gp = gpar(col='black',lwd=0.05),
          col = color_function,
          show_column_names = F
          )
  )
  dev.off()
}
pval_heatmap(pval_df, paste0(project,'_single_gene_',end,'_adjpvalues_clust.png'), 'P-values')
pval_heatmap(cluster_df, paste0(project,'_single_gene_',end,'_adjpvalues_cox_cluster_clust.png'),'Cox Cluster P-values')
pval_heatmap(age_df, paste0(project,'_single_gene_',end,'_adjpvalues_cox_age_clust.png'),'Cox Age P-values')
pval_heatmap(gender_df, paste0(project,'_single_gene_',end,'_adjpvalues_cox_gender_clust.png'),'Cox Gender P-values')

# heatmap without clutering
pval_heatmap_no_cluster <- function(pval_df, filename, val_type){
  if(sum(is.na(pval_df))==(nrow(pval_df)*ncol(pval_df))){
    return('df contains all NA values')
  }
  png(filename,units='in',width=7,height=7,res=500)
  h <- draw(Heatmap(as.matrix(pval_df), cluster_rows = F, cluster_columns = F,
                    name = 'adj p-value',
                    heatmap_legend_param = list(title='p-value',at=c(0,0.05),labels=c('0','0.05')),
                    column_title = paste0('Single Gene ',end,' ',val_type),
                    rect_gp = gpar(col='black',lwd=0.05),
                    na_col='white',
                    col = color_function,
                    show_column_names = F
        )
  )
  dev.off()
}
pval_heatmap_no_cluster(pval_df, paste0(project,'_single_gene_',end,'_adjpvalues_no_clust.png'), 'P-values')
pval_heatmap_no_cluster(cluster_df, paste0(project,'_single_gene_',end,'_adjpvalues_cox_cluster_no_clust.png'),'Cox Cluster P-values')
pval_heatmap_no_cluster(age_df, paste0(project,'_single_gene_',end,'_adj_pvalues_cox_age_no_clust.png'),'Cox Age P-values')
pval_heatmap_no_cluster(gender_df, paste0(project,'_single_gene_',end,'_adj_pvalues_cox_gender_no_clust.png'),'Cox Gender P-values')

# heatmap with row and column ordered by most number of significant epigenes/cancer types, respectively
pval_heatmap_desc <- function(pval_df, filename, val_type){
  if(sum(is.na(pval_df))==(nrow(pval_df)*ncol(pval_df))){
    return('df contains all NA values')
  }
  number_na <- function(l) {
    return(sum(is.na(l)))
  }
  row_order <- order(apply(pval_df,MARGIN=1, number_na))
  col_order <- order(apply(pval_df, MARGIN=2, number_na))
  order_pval_df <- pval_df[row_order,col_order]
  png(filename,units='in',width=7,height=7,res=500)
  h <- draw(Heatmap(as.matrix(order_pval_df), cluster_rows = F, cluster_columns = F,
                    name = 'adj p-value',
                    heatmap_legend_param = list(title='p-value',at=c(0,0.05),labels=c('0','0.05')),
                    column_title = paste0('Single Gene ',end,' ',val_type),
                    rect_gp = gpar(col='black',lwd=0.05),
                    na_col = 'white',
                    col = color_function,
                    show_column_names = F
  )
  )
  dev.off()
  
}

pval_heatmap_desc(pval_df, paste0(project,'_single_gene_',end,'_adjpvalues_descending_orders.png'), 'P-values')
pval_heatmap_desc(cluster_df, paste0(project,'_single_gene_',end,'_adjpvalues_cox_cluster_descending_orders.png'),'Cox Cluster P-values')
pval_heatmap_desc(age_df, paste0(project,'_single_gene_',end,'_adj_pvalues_cox_age_descending_orders.png'),'Cox Age P-values')
pval_heatmap_desc(gender_df, paste0(project,'_single_gene_',end,'_adj_pvalues_cox_gender_descending_orders.png'),'Cox Gender P-values')


# write p-value table
write.csv(pval_df, file = paste0(project,'_single_gene_',end,'_adjpvalues.csv'),row.names=T, quote=F)
write.csv(pval_df[,which(apply(is.na(pval_df),2,sum)<4)], file = paste0(project,'_single_gene_',end,'_adjpvalues_filtered.csv'),row.names=T, quote=F)
write.csv(cluster_df, file = paste0(project,'_single_gene_',end,'_adjpvalues_cox_cluster.csv'),row.names=T, quote=F)
write.csv(cluster_df[,which(apply(is.na(cluster_df),2,sum)<4)], file = paste0(project,'_single_gene_',end,'_adjpvalues_cox_cluster_filtered.csv'),row.names=T, quote=F)

write.csv(age_df, file = paste0(project,'_single_gene_',end,'_adjpvalues_cox_age.csv'),row.names=T, quote=F)
write.csv(gender_df, file = paste0(project,'_single_gene_',end,'_adjpvalues_cox_gender.csv'),row.names=T, quote=F)


###########
# heatmap of gene effects on endpoint
###########

# heatmap
png(paste0(project,'_single_gene_effect_',end,'.png'),units='in',width=7,height=7,res=500)
Heatmap(as.matrix(effect_df),
        name = paste0('effect on ',end),
        column_title = paste0('Gene Effect on ',end,' in Cancers'),
        col = structure(c('blue','white','red'), names = c(-1,0,1)),
        rect_gp = gpar(col='black',lwd=0.05),
        heatmap_legend_param = list(at=c(-1,0,1),labels=c(paste0(end,'-inhibiting'),'none',paste0(end,'-promoting'))),
        show_column_names = F
        )
dev.off()
# write effect table
write.csv(effect_df,paste0(project,'_single_gene_effect_',end,'.csv'),row.names=T,quote=F)

#######
# Write Excel Sheets
#######
library(stringr)
write_excel_sheet <- function(pval_df, filename, val_type){
  if(sum(is.na(pval_df))==(nrow(pval_df)*ncol(pval_df))){
    return('df contains all NA values')
  }
  # pval_df
  # create workbook
  wb <- createWorkbook(type='xlsx')
  
  # initialize a cell style. DEFAULT SETTING: CellStyle(wb, dataFormat = NULL, alignment=NULL, border=NULL, fill=NULL, font=NULL)
  # Row and Columns style
  ROWCOL_STYLE <- CellStyle(wb) + Font(wb, isBold = F)
  
  # create new sheet to write
  sheet <- createSheet(wb, sheetName = paste0(end,"-Sig Epigene ",val_type))
  
  # add data frame
  # pval_df[is.na(pval_df)] <- 0
  addDataFrame(pval_df, sheet=sheet, colnamesStyle = ROWCOL_STYLE, rownamesStyle = ROWCOL_STYLE)
  rows <- getRows(sheet=sheet) # get all rows
  cells <- getCells(rows) # get all non empty cells
  values <- lapply(cells, getCellValue) # extract values
  
  # set style color for each cell according to the endpoint effect
  cs_promote <- CellStyle(wb)+Fill(foregroundColor='lightgreen',backgroundColor='lightgreen')
  cs_inhibit <- CellStyle(wb)+Fill(foregroundColor='indianred1',backgroundColor='indianred1')
  
  # function to color excel sheet according to the effect of the gene in the cancer
  endpoint_effect_style <- function(value_name){
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
  color_coding <- sapply(names(values),endpoint_effect_style)
  
  # save workbook to file
  saveWorkbook(wb, filename)
  
}
write_excel_sheet(pval_df, paste0(project,'_single_gene_adjpvalue_and_',end,'_effect.xlsx'), 'P-values')
write_excel_sheet(cluster_df, paste0(project,'_single_gene_adjpvalue_cox_cluster_and_',end,'_effect.xlsx'), 'Cox Cluster P-values')
write_excel_sheet(age_df, paste0(project,'_single_gene_adjpvalue_cox_age_and_',end,'_effect.xlsx'), 'Cox Age P-values')
write_excel_sheet(gender_df, paste0(project,'_single_gene_adjpvalue_cox_gender_and_',end,'_effect.xlsx'), 'Cox Gender P-values')


# cox corrected cluster pvalue dataframe, setting negative pvalues for os-inhibiting effect (-1 in effect_df) (high expression/poor outcome)
cluster_df_sign <- cluster_df * effect_df
write.csv(cluster_df_sign, file = paste0(project,'_single_gene_',end,'_adjpvalues_cox_cluster_with_sign.csv'),row.names=T, quote=F)

