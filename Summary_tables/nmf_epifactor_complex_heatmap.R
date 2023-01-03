library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(reshape2)
library(RColorBrewer)
library(circlize)

project <- 'TARGET' # TCGA or TARGET
analysis <- 'single_gene' # NMF or single gene

if(project == 'TARGET'){
  end <- 'OS'
  if(file.exists(paste0('../../Summary_tables/',project,'_nmf_top_genes_all.csv'))==F){
    # top genes df
    analysis_genes <- data.frame()
    cancers <- c('AML','NBL','OS','WT')
    for(c in cancers){
      top_df <- data.frame('Gene'=c(rownames(read.csv(paste0('../../',project,'_',c,'/03_nmf/Rank_2/nmf_lee_rank2_feature1_genes.csv'),row.names=1)),
                                    rownames(read.csv(paste0('../../',project,'_',c,'/03_nmf/Rank_2/nmf_lee_rank2_feature2_genes.csv'),row.names=1))),
                           'Cancer_Type'=c)
      analysis_genes <- rbind(analysis_genes,top_df)
    }
    write.csv(analysis_genes,paste0('../../Summary_tables/',project,'_nmf_top_genes_all.csv'),row.names=F)
  }
} else{
  end <- 'PFI'
}

if (analysis == 'NMF') {
  analysis_genes <- read.csv(paste0('../../Summary_tables/',project,'_nmf_top_genes_all.csv'), header=T)
  a_title <- analysis
} else if (analysis == 'single_gene'){
  analysis_genes <- read.csv(paste0('../../Summary_tables/',project,'_single_gene_',end,'_adjpvalues_cox_cluster.csv'))
  analysis_genes <- melt(analysis_genes, id.vars = 'X', variable.name = 'Gene', value.name = 'pval')
  analysis_genes <- analysis_genes[is.na(analysis_genes$pval)==F, c(2,1,3)]
  colnames(analysis_genes)[2] <- 'Cancer_Type'
  rownames(analysis_genes) <- NULL
  a_title <- 'Prognostic'
}


# find enriched complexes and groups for significant genes
#################### Formatting epigenetic complex dataframe with columns ##############################
#################### explaining the gene and the epigenetic complexes it is involved in ####################
# map uniprot to hugo genes and complexes
complexes <- read.delim('../../metadata/epigenes.csv')
complexes$Function <- gsub(',', ';', complexes$Function)
# hugo to uniprot conversion vector
hugo_to_uniprot <- complexes[,c("HGNC.approved.symbol", "UniProt.ID..human.")]
rownames(hugo_to_uniprot) <- complexes$HGNC.approved.symbol
# give each protein complex its own row
complexes$Protein.complex <- gsub(',', ';', complexes$Protein.complex)
protein_complex_by_gene <- str_split(complexes$Protein.complex,", ")
names(protein_complex_by_gene) <- complexes$HGNC.approved.symbol

complex_df <- melt(protein_complex_by_gene)[,c(2,1)]
colnames(complex_df) <- c("HUGO_gene","Protein_complex")
complex_df$Uniprot <- hugo_to_uniprot[complex_df$HUGO_gene,2]

genes <- unique(analysis_genes$Gene)
# find main protein complexes and functional types
# subset epigene complex df with top_genes and order by protein complex
complex_subset <- complex_df %>% filter(HUGO_gene%in%genes) %>% select(HUGO_gene, Protein_complex)
# # eliminate NA complexes
complex_subset <- complex_subset %>% filter(Protein_complex!="#") %>% arrange(Protein_complex)
#########################################################################################

############### Make dataframe with columns explaining the genes and ######################
############### the epigenetic complex family/group it is involved in #####################

# map uniprot to protein groups
groups <- read.csv('../../metadata/EpiGenes_complexes.csv', header=T, row.names=1)
# id to group
id_to_group <- groups[,"Group_name"]
names(id_to_group) <- rownames(groups)
# separate uniprot id's 
group_by_uniprot <- str_split(groups$UniProt_ID, ', ')
names(group_by_uniprot) <- rownames(groups)
group_df <- melt(group_by_uniprot)[,c(2,1)]
colnames(group_df) <- c("ID","Uniprot")
group_df$Group <- id_to_group[group_df$ID]
# shorten group names
group_df['Group'][group_df['Group']=='APOB mRNA editosome'] <- 'APOB'
group_df['Group'][group_df['Group']=='PcG and PcG-like'] <- 'PcG'

uniprot_to_hugo <- hugo_to_uniprot%>%filter(!duplicated(UniProt.ID..human.))
rownames(uniprot_to_hugo) <- uniprot_to_hugo$UniProt.ID..human.
group_df$HUGO_gene <- uniprot_to_hugo[group_df$Uniprot,1]
group_subset <- group_df %>% filter(HUGO_gene%in%genes) %>% select(HUGO_gene, Group)
#########################################################################################

cancers <- unique(analysis_genes$Cancer_Type)
#### Visualizations
# Protein Complexes

# initialize cancer_by_complex data frame denoting number of nmf genes in that complex
cancer_by_complex <- as.data.frame(matrix(0,nrow=length(cancers), ncol=length(unique(complex_subset$Protein_complex))))
rownames(cancer_by_complex) <- cancers
colnames(cancer_by_complex) <- unique(complex_subset$Protein_complex)

# for each cancer, find number of top nmf genes in each complex
for (c in cancers) {
  # top nmf genes
  top_genes <- (analysis_genes%>%filter(Cancer_Type==c))[,"Gene"]
  top_complexes_freq <- table((complex_subset %>% filter(HUGO_gene%in%top_genes))[,"Protein_complex"])
  cancer_by_complex[c,names(top_complexes_freq)] <- top_complexes_freq
}

# normalize cancer_by_complex values by the number of genes in each epifactor complex
total_complex_genes <- apply(cancer_by_complex,2,sum)

normalize <- function(complex){
  cancer_by_complex[,complex] <<- cancer_by_complex[,complex]/total_complex_genes[complex]
}

invisible(sapply(names(total_complex_genes),normalize))

# heatmap
library(RColorBrewer)
col_fun = brewer.pal(n=9,name="Reds")
# col_fun(seq(-3, 3))
png(paste0("../../Summary_tables/",project,"_",analysis,"_top_genes_epigenetic_complexes.png"), units="in", width=7, height=7, res=500)
Heatmap(as.matrix(cancer_by_complex),
        name = 'count',
        column_title = paste0('Top ',a_title,' Genes in Each Complex'),
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=5),
        column_names_rot = 45,
        col=col_fun)
dev.off()


#################### Complex Groups ######################
# size of each complex group
group_size <- table(group_df$Group)
# initialize vector for total number of genes in analysis for each cancer type
total_epigenes <- rep(0,length(cancers))
names(total_epigenes) <- cancers
# initialize cancer_by_group data frame denoting number of top genes in that group
cancer_by_group <- as.data.frame(matrix(0,nrow=length(cancers), ncol=length(unique(group_df$Group))))
rownames(cancer_by_group) <- cancers
colnames(cancer_by_group) <- unique(group_df$Group)
# initialize df denoting total number of epigenes in the group for each cancer
total_by_group <- cancer_by_group

# for each cancer, find number of top genes in each group and number of epigenes used in analysis
for (c in cancers) {
  # top genes
  top_genes <- (analysis_genes%>%filter(Cancer_Type==c))[,"Gene"]
  top_groups_freq <- table((group_subset %>% filter(HUGO_gene%in%top_genes))[,"Group"])
  cancer_by_group[c,names(top_groups_freq)] <- top_groups_freq
  # total epigenes
  epigenes_directory <- paste0('../../',project,'_',c,'/RNA-seq_datasets/')
  epigenes_file <- list.files(epigenes_directory,pattern = paste0(c,'_epigenes_log2norm_sd'))
  epigenes_df <- read.csv(paste0(epigenes_directory,epigenes_file),row.names=1)
  epigenes <- rownames(epigenes_df)
  total_epigenes[c] <- length(epigenes)
  
  epi_group_freq <- table((group_df %>% filter(HUGO_gene%in%epigenes))[,"Group"])
  total_by_group[c,names(epi_group_freq)] <- epi_group_freq
  }

# normalize by total number of top epigenes in cancer (i.e. row-sums = 1)
total_cancer_genes <- table(analysis_genes$Cancer_Type)
normalized_cancer_by_group <- cancer_by_group
normalize <- function(cancer){
  normalized_cancer_by_group[cancer,] <<- normalized_cancer_by_group[cancer,]/total_cancer_genes[cancer]
}
invisible(sapply(names(total_cancer_genes),normalize))
rownames(normalized_cancer_by_group) <- mapply(paste0, rownames(normalized_cancer_by_group),'(', total_cancer_genes[rownames(normalized_cancer_by_group)],')')
colnames(normalized_cancer_by_group) <- mapply(paste0, colnames(normalized_cancer_by_group),'(', group_size[colnames(normalized_cancer_by_group)],')')

# heatmap (each cell is colored by the percentage of top genes of the cancer that is in the complex group)
col_fun = brewer.pal(n=9,name="Reds")
png(paste0("../../Summary_tables/",project,"_",analysis,"_top_genes_epigenetic_groups_rownorm.png"), units="in", width=7, height=7, res=500)
h <- Heatmap(normalized_cancer_by_group,
        name = 'count',
        column_title = paste0('Top ',a_title,' Genes in Each Complex Group'),
        row_names_gp = gpar(fontsize=7),
        column_names_gp = gpar(fontsize=7),
        column_names_rot = 90,
        col=col_fun)
ht <- draw(h)
dev.off()
write.csv(normalized_cancer_by_group[row_order(ht),column_order(ht)],paste0("../../Summary_tables/",project,"_",analysis,"_top_genes_epigenetic_groups_rownorm.csv") )


##### Fisher's exact test for each cancer and complex
# function to create contingency table and fisher test
cont_table <- function(cancer,cgroup){
                       # ,total_by_group=total_by_group, 
                       # cancer_by_group=cancer_by_group, group_size=group_size, total_cancer_genes=total_cancer_genes){
  # make contingency table
  ct <- matrix(0,nrow=2,ncol=2)
  rownames(ct) <- c('In Group','Not In Group')
  colnames(ct) <- c('Top NMF', 'Not Top NMF')
  ct['In Group','Top NMF'] <- cancer_by_group[cancer,cgroup]
  ct['In Group','Not Top NMF'] <- total_by_group[cancer,cgroup]-ct['In Group','Top NMF']
  ct["Not In Group","Top NMF"] <- total_cancer_genes[cancer]-ct['In Group','Top NMF']
  ct["Not In Group","Not Top NMF"] <- total_epigenes[cancer]-total_cancer_genes[cancer]-ct["In Group","Not Top NMF"]
  print(ct)
    # fishers
  fish <- fisher.test(ct)
  # pval <- format(signif(fish$p.value,5),scientific = T)
  return(fish$p.value)
  }
v_cont_table <- Vectorize(cont_table)
fisher_heatmap <- outer(rownames(cancer_by_group),colnames(cancer_by_group), v_cont_table)
rownames(fisher_heatmap) <- mapply(paste0,rownames(cancer_by_group),'(',total_cancer_genes[rownames(cancer_by_group)],')')
colnames(fisher_heatmap) <- mapply(paste0,colnames(cancer_by_group),'(',group_size[colnames(cancer_by_group)],')')
fisher_heatmap <- t(apply(fisher_heatmap,1,p.adjust,method="BH")) # BH correction within each cancer type
# heatmap
color_border <- brewer.pal(n=9,"Reds")[c(8,1)]
col_fun <- colorRamp2(c(0,0.049999999,0.05),c(color_border,'grey'))
png(paste0("../../Summary_tables/",project,"_",analysis,"_top_genes_epigenetic_groups_fishers.png"), units="in", width=15, height=10, res=500)
h <- Heatmap(fisher_heatmap,
        name = 'pvalue',
        heatmap_legend_param = list(title='p-value', at=c(0:5)/100, labels=c(c(0:4)/100,'>=0.05'),
                                    labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                    title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                    legend_height = unit(2.5, "in")),
        column_title = 'Epifactor Group (# of genes)',
        column_title_gp = gpar(fontsize=20,fontfamily='Helvetica'),
        column_title_side = "bottom",
        row_title = 'Cancer Type (# of genes)',
        row_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
        row_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
        column_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
        column_names_rot = 60,
        col=col_fun)
# add more titles to heatmap
ht <- draw(h,
     column_title=paste0("Top ",a_title," Genes in Each Complex Group Fisher's P-value"), column_title_gp=gpar(fontsize=25, fontfamily='Helvetica')
)
ht
dev.off()
write.csv(fisher_heatmap[row_order(ht),column_order(ht)],paste0('../../Summary_tables/',project,'_',analysis,'_top_genes_epigenetic_groups_fishers.csv'))

#### find most significant genes in SWI/SNF for TCGA LIHC and ACC for single gene analysis
if (analysis=='single_gene' & project=='TCGA'){
  group_list <- c('SWI/SNF','HAT','SWI/SNF','MLL')
  cancer_list <- c('LIHC','LIHC','ACC','ACC')
  for (i in 1:length(cancer_list)){
    group_genes <- group_subset %>% filter(Group==group_list[i])
    c_df <- analysis_genes %>% filter(Cancer_Type==cancer_list[i], Gene %in% group_genes$HUGO_gene) %>% arrange(pval)
    write.csv(c_df,paste0('../../',project,'_',cancer_list[i],'/05_gene_clinical_prediction_analysis/',gsub('/','',group_list[i]),'_sig_genes.csv'))
  }
}


# raw counts heatmap
rownames(cancer_by_group) <- rownames(fisher_heatmap)
colnames(cancer_by_group) <- colnames(fisher_heatmap)
# add column for nmf genes that are not in any complex group
cancer_by_group$`nogroup(NA)` <- total_cancer_genes[cancers] - apply(cancer_by_group,1,sum)

col_fun <- brewer.pal(n=9,"Reds")
png(paste0("../../Summary_tables/",project,"_",analysis,"_top_genes_epigenetic_groups_raw.png"), units="in", width=15, height=10, res=500)
h <- Heatmap(cancer_by_group,
        name = 'Count',
        heatmap_legend_param = list(labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                    title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                    legend_height = unit(2.5, "in")),
        column_title = 'Epifactor Group (# of genes)',
        column_title_gp = gpar(fontsize=20,fontfamily='Helvetica'),
        column_title_side = "bottom",
        row_title = 'Cancer Type (# of genes)',
        row_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
        row_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
        column_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
        column_names_rot = 60,
        col=col_fun)
# add more titles to heatmap
ht <- draw(h,
     column_title=paste0('Top ',a_title,' Genes in Each Complex Group'), column_title_gp=gpar(fontsize=25, fontfamily='Helvetica')
     )
ht
dev.off()
write.csv(cancer_by_group[row_order(ht),column_order(ht)],paste0("../../Summary_tables/",project,"_",analysis,"_top_genes_epigenetic_groups_raw.csv") )


#################################
## top genes by cancer type

analysis_genes$Number_Cancer_Types <- table(analysis_genes$Gene)[analysis_genes$Gene]
if(project=='TCGA'){
  common_genes <- analysis_genes%>%filter(Number_Cancer_Types>=7)
} else{
  common_genes <- analysis_genes%>%filter(Number_Cancer_Types>=2)
}

gene_nmf_matrix <- matrix(0,nrow=length(unique(common_genes$Gene)), ncol=length(cancers))
rownames(gene_nmf_matrix) <- unique(common_genes$Gene)
colnames(gene_nmf_matrix) <- cancers

for(i in 1:nrow(common_genes)){
  gene_nmf_matrix[as.character(common_genes[i,"Gene"]),common_genes[i,"Cancer_Type"]] <- 1
}

# col_fun = structure(c(1,2),names=c("Not Present", "Present"))
col_fun = brewer.pal(n=9,name="Reds")[c(1,9)]
png(paste0("../../Summary_tables/",project,"_",analysis,"_common_genes_heatmap.png"), units="in", width=15, height=10, res=500)
h <- Heatmap(gene_nmf_matrix,
        name = ' ',
        heatmap_legend_param = list(labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                    labels=c('Included','Not Included'), border='black',
                                    title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                    legend_height = unit(2.5, "in")),
        column_title = 'Cancer Type',
        column_title_gp = gpar(fontsize=20,fontfamily='Helvetica'),
        column_title_side = 'bottom',
        row_title = 'Genes',
        row_title_gp = gpar(fontsize=20,fontfamily='Helvetica'),
        row_names_gp = gpar(fontsize=15,fontfamily='Helvetica'),
        column_names_gp = gpar(fontsize=15,fontfamily='Helvetica'),
        column_names_rot = 90,
        col=col_fun)
# add more titles to heatmap
ht <- draw(h,
     column_title=paste0('Common ',a_title,' Genes in Cancer Types'), column_title_gp=gpar(fontsize=25, fontfamily='Helvetica')
)
ht
dev.off()
write.csv(gene_nmf_matrix[row_order(ht),column_order(ht)], paste0("../../Summary_tables/",project,"_",analysis,"_common_genes_heatmap.csv"))


