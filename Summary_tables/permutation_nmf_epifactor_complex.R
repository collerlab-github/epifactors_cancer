library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(reshape2)
library(RColorBrewer)
library(circlize)

project <- 'TCGA' # TCGA or TARGET
analysis <- 'NMF' # NMF or single gene

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

# read in significant genes
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

genes <- unique(analysis_genes$Gene)
cancers <- unique(analysis_genes$Cancer_Type)

# map uniprot to hugo genes and complexes
complexes <- read.delim('../../metadata/epigenes.csv')
# hugo to uniprot conversion vector
hugo_to_uniprot <- complexes[,c("HGNC.approved.symbol", "UniProt.ID..human.")]
rownames(hugo_to_uniprot) <- complexes$HGNC.approved.symbol

# find enriched epigenetic protein groups for significant genes

#########################################################################################

############### Make dataframe with columns explaining the genes and ######################
############### the epigenetic complex family/group it is involved in #####################

# map uniprot to groups (family of epifactor proteins e.g. SWI/SNF)
groups <- read.csv('../../metadata/EpiGenes_protein_complexes.csv', sep='\t', header=T)[,c('Complex.group.name','Genes.in.complex')]
colnames(groups) <- c('group','genes')
# separate genes
group_by_genes <- str_split(groups$genes, ', ')
names(group_by_genes) <- groups$group
# remove genes with '?' or '()'
group_by_genes <- lapply(group_by_genes, function(x){x[!grepl('[\\(,\\),\\?]',x)]})

group_df <- melt(group_by_genes)[,c(2,1)]
colnames(group_df) <- c("group","gene")
group_df <- unique(group_df)
# replace ',' and ' ' in group names
group_df$group <- gsub(' ','_',group_df$group)
group_df$group <- gsub(',','',group_df$group)


# Include all epifactor genes in df
# Genes with no assigned complex are in the "None" group.
epigenes <- unique(read.csv('../../metadata/epigenes.csv', sep='\t')$HGNC.approved.symbol)
group_df <- group_df %>% filter(gene %in% epigenes)
group_df <- rbind(group_df, data.frame(group='None',gene=epigenes[epigenes%in%group_df$gene==F]))
group_df$group_size <- table(group_df$group)[group_df$group]
group_df$gene_freq <- table(group_df$gene)[group_df$gene]

group_df <- group_df %>% arrange(desc(gene_freq),gene, desc(group_size))
# Genes assigned to multiple complexes are also assigned "None" group.
group_df$revised_group <- group_df$group
group_df$revised_group[group_df$gene_freq>1] <- 'None'
group_df$revised_group_size <- table(group_df$revised_group)[group_df$revised_group]
write.csv(group_df,'../../metadata/epigenes_protein_complex_genes.csv', row.names=F, quote=F)


# summary numbers for gene frequencies
length(table(group_df$gene)) # 240 genes map to at least 1 complex
sum(genes %in% group_df$gene) # 105/345 NMF, 224/687 prognostic genes map to at least 1 complex
sum(group_df$gene_freq==1) # 202 genes map to 1 complex
sum(genes %in% (group_df%>%filter(gene_freq==1))[,'gene']) # 85/345NMF, 186/687 prognostic genes map to only 1 complex
sum(group_df$gene_freq==2)/2 # 30 genes map to 2 complexes
sum(group_df$gene_freq==3)/3 # 2 genes map to 3 complexes
sum(group_df$gene_freq==4)/4 # 6 genes map to 4 complexes

########################################################################################
# shorten group names
group_df['revised_group'][group_df['revised_group']=='APOB_mRNA_editosome'] <- 'APOB'
group_df['revised_group'][group_df['revised_group']=='PcG_and_PcG-like'] <- 'PcG'
group_df['revised_group'][group_df['revised_group']=='CRL4-DDB2_E3_ubiquitin_ligase_complex_CUL4A_CUL4B_variants'] <- 'CRL4'


# unique entries only
group_df <- unique(group_df[,c("revised_group","gene")])
colnames(group_df)[1] <- 'group'
# group_subset <- group_df %>% filter(gene%in%genes) %>% select(gene, group)
#########################################################################################

# size of each group
group_size <- table(group_df$group)
# initialize vector for total number of genes in analysis for each cancer type
total_epigenes <- rep(0,length(cancers))
names(total_epigenes) <- cancers
# initialize cancer_by_group data frame denoting number of top genes in that group
cancer_by_group <- as.data.frame(matrix(0,nrow=length(cancers), ncol=length(unique(group_df$group))))
rownames(cancer_by_group) <- cancers
colnames(cancer_by_group) <- unique(group_df$group)
# initialize df denoting total number of epigenes in the group for each cancer
total_by_group <- cancer_by_group
# initialize fisher's exact p-value dataframe
fisher_by_group <- as.data.frame(matrix(NA,nrow=length(cancers), ncol=length(unique(group_df$group))))
rownames(fisher_by_group) <- cancers
colnames(fisher_by_group) <- unique(group_df$group)
# initialize permutation p-value dataframe
perm_pval_df <- fisher_by_group


# function to find number of genes in each group for a given gene set
calculate_overlap = function(genes, group_df){
  group_freq <- table((group_df %>% filter(gene%in%genes))[,"group"])
  return(group_freq)
}

##### Fisher's exact test for each cancer and complex. 
##### The number of genes in an epigenetic complex is compared between top NMF epigenes and not top NMF epigenes.
# function to create contingency table and fisher test
cont_table <- function(cgroup, top_group_overlap, total_group_overlap, num_top_genes, num_epigenes, verbose=F){
  # cgroup: complex group to calculate fisher's on
  # top_group_overlap: vector of number of top nmf genes overlapping with groups
  # total_group_overlap: vector of number of total analysis epigenes overlapping with groups
  # num_top_genes: number of top NMF genes
  # num_epigenes: number of analysis epigenes
  
  # make contingency table
  ct <- matrix(0,nrow=2,ncol=2)
  rownames(ct) <- c('In Group','Not In Group')
  colnames(ct) <- c('Top NMF', 'Not Top NMF')
  
  ct['In Group','Top NMF'] <- ifelse(cgroup %in% names(top_group_overlap),top_group_overlap[cgroup],0)
  ct['In Group','Not Top NMF'] <- ifelse(cgroup%in% names(total_group_overlap),total_group_overlap[cgroup]-ct['In Group','Top NMF'],0)
  ct["Not In Group","Top NMF"] <- num_top_genes-ct['In Group','Top NMF']
  ct["Not In Group","Not Top NMF"] <- num_epigenes-num_top_genes-ct["In Group","Not Top NMF"]
  if(verbose){
    print(ct)
  }
  # fishers
  fish <- fisher.test(ct)
  
  return(as.numeric(fish$estimate))
  # pval <- format(signif(fish$p.value,5),scientific = T)
  # return(fish$p.value)
}

# function to calculate permutation p-value
calc_pval <- function(group_name, top_group_overlap, permutation_top_gene_group){
  # group_name: group to calculate pvalue for
  # top_group_overlap: vector containing number of overlapping top NMF genes in each group
  # permutation_top_gene_group: overlap of permuted random top NMF genes for all permutations
  
  # pvalue calculated as proportion of permutations with overlap number >= overlap of true top NMF genes
  
  pval <- sum(permutation_top_gene_group[,group_name] >= 
                ifelse(group_name %in% names(top_group_overlap),top_group_overlap[group_name],0))/nrow(permutation_top_gene_group)
  return(pval)
}
# for each cancer,
num_perm=10000
set.seed(0)
for (c in cancers) {
  # number of top nmf genes
  num_top_genes <- table(analysis_genes$Cancer_Type)[c]
  # total epigenes
  epigenes_directory <- paste0('../../',project,'_',c,'/RNA-seq_datasets/')
  epigenes_file <- list.files(epigenes_directory,pattern = paste0(c,'_epigenes_log2norm_sd'))
  epigenes_df <- read.csv(paste0(epigenes_directory,epigenes_file),row.names=1)
  epigenes <- rownames(epigenes_df)
  total_epigenes[c] <- length(epigenes)
  # overlap of analysis epigenes and groups' genes
  epi_group_overlap <- calculate_overlap(epigenes, group_df)
  # initialize dataframes to record numbers for each permutation
  # permutation_top_gene_group <- as.data.frame(matrix(0,nrow=num_perm, ncol=length(unique(group_df$group))))
  # colnames(permutation_top_gene_group) <- unique(group_df$group)
  permutation_fisher_df <- as.data.frame(matrix(1,nrow=num_perm, ncol=length(unique(group_df$group))))
  colnames(permutation_fisher_df) <- unique(group_df$group)
  # permutation test
  for(p in seq(num_perm)){
    # randomly select top NMF genes from analysis epigenes
    top_genes <- sample(epigenes,size = num_top_genes, replace = F)
    top_overlap <- calculate_overlap(top_genes, group_df)
    # permutation_top_gene_group[p,names(top_overlap)] <- top_overlap
    
    # print(top_genes)
    # print(top_overlap)
    
    # # calculate fisher's exact test score for each complex group
    fish_groups <- sapply(colnames(permutation_fisher_df), cont_table, top_overlap, epi_group_overlap, num_top_genes, length(epigenes))
    permutation_fisher_df[p,names(fish_groups)] <- fish_groups
    
    
  }
  # actual results
  top_genes <- (analysis_genes%>%filter(Cancer_Type==c))[,"Gene"]
  top_overlap <- calculate_overlap(top_genes, group_df)
  cancer_by_group[c,names(top_overlap)] <- top_overlap
  total_by_group[c,names(epi_group_overlap)] <- epi_group_overlap
  fish_groups <- sapply(names(top_overlap), cont_table, top_overlap, epi_group_overlap, num_top_genes, length(epigenes)) 
  fisher_by_group[c,names(fish_groups)] <- fish_groups
  
  # IMPLEMENT PVALUE CALCULATION FUNCTION (PROPORTION OF TESTS WITH GREATER OR EQUAL TO TRUE OVERLAP)
  # pvals <- sapply(names(permutation_top_gene_group),calc_pval, top_overlap, permutation_top_gene_group)
  pvals <- sapply(names(permutation_fisher_df), calc_pval, top_overlap, permutation_fisher_df)
  # UPDATE DATA FRAME
  perm_pval_df[c,names(pvals)] <- pvals
  # break
}

write.csv(cancer_by_group,paste0('../../Summary_tables/',project,'_',analysis,'_top_genes_epigenetic_groups_overlap.csv'))


fisher_by_group[is.na(fisher_by_group)] <- 0
fisher_by_group[sapply(fisher_by_group, is.infinite)] <- max(fisher_by_group[!sapply(fisher_by_group, is.infinite)])+100
# # normalize by total number of top epigenes in cancer (i.e. row-sums = 1)
# total_cancer_genes <- table(analysis_genes$Cancer_Type)
# normalized_cancer_by_group <- cancer_by_group
# normalize <- function(cancer){
#   normalized_cancer_by_group[cancer,] <<- normalized_cancer_by_group[cancer,]/total_cancer_genes[cancer]
# }
# invisible(sapply(names(total_cancer_genes),normalize))
# rownames(normalized_cancer_by_group) <- mapply(paste0, rownames(normalized_cancer_by_group),'(', total_cancer_genes[rownames(normalized_cancer_by_group)],')')
# colnames(normalized_cancer_by_group) <- mapply(paste0, colnames(normalized_cancer_by_group),'(', group_size[colnames(normalized_cancer_by_group)],')')

# heatmap
# color_border <- brewer.pal(n=9,"Reds")[c(8,1)]
# col_fun <- colorRamp2(c(0,0.049999999,0.05),c(color_border,'grey'))
color_border <- brewer.pal(n=9,"Reds")
png(paste0("../../Summary_tables/",project,"_",analysis,"_top_genes_epigenetic_groups_fishers_OR.png"), units="in", width=15, height=10, res=500)
h <- Heatmap(fisher_by_group,
             name = 'OR',
             # heatmap_legend_param = list(title='p-value', at=c(0:5)/100, labels=c(c(0:4)/100,'>=0.05'),
             #                             labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
             #                             title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
             #                             legend_height = unit(2.5, "in")),
             column_title = 'Epifactor Group (# of genes)',
             column_title_gp = gpar(fontsize=20,fontfamily='Helvetica'),
             column_title_side = "bottom",
             row_title = 'Cancer Type (# of genes)',
             row_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
             row_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
             column_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
             column_names_rot = 60,
             col=color_border)
# add more titles to heatmap
ht <- draw(h,
           column_title=paste0("Top ",a_title," Genes in Each Complex Group Fisher's Odds-Ratio"), column_title_gp=gpar(fontsize=25, fontfamily='Helvetica')
)
ht
dev.off()
write.csv(fisher_by_group[row_order(ht),column_order(ht)],paste0('../../Summary_tables/',project,'_',analysis,'_top_genes_epigenetic_groups_fishers_OR.csv'))



# heatmap for permutation test
color_border <- brewer.pal(n=9,"Reds")[c(8,1)]
col_fun <- colorRamp2(c(0,0.049999999,0.05),c(color_border,'grey'))
png(paste0("../../Summary_tables/",project,"_",analysis,"_top_genes_epigenetic_groups_perm_overlap_OR.png"), units="in", width=15, height=10, res=500)
h <- Heatmap(perm_pval_df,
             name = 'pvalue',
             heatmap_legend_param = list(title='p-value', at=c(0:5)/100, labels=c(c(0:4)/100,'>=0.05'),
                                         labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                         title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                         legend_height = unit(2.5, "in")),
             column_title = 'Epifactor Group',
             column_title_gp = gpar(fontsize=20,fontfamily='Helvetica'),
             column_title_side = "bottom",
             row_title = 'Cancer Type',
             row_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
             row_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
             column_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
             column_names_rot = 60,
             col=col_fun)
# add more titles to heatmap
ht <- draw(h,
           column_title=paste0("Top ",a_title," Genes Permutation Odds Ratio P-value"), column_title_gp=gpar(fontsize=25, fontfamily='Helvetica')
)
ht
dev.off()
write.csv(perm_pval_df[row_order(ht),column_order(ht)],paste0('../../Summary_tables/',project,'_',analysis,'_top_genes_epigenetic_groups_perm_overlap_OR.csv'))

# heatmap for p adjusted permutation test
fisher_heatmap <- t(apply(perm_pval_df,1,p.adjust,method="BH")) # BH correction within each cancer type
png(paste0("../../Summary_tables/",project,"_",analysis,"_top_genes_epigenetic_groups_perm_overlap_OR_BH.png"), units="in", width=15, height=10, res=500)
h <- Heatmap(fisher_heatmap,
             name = 'pvalue',
             heatmap_legend_param = list(title='p-value', at=c(0:5)/100, labels=c(c(0:4)/100,'>=0.05'),
                                         labels_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                         title_gp=gpar(fontsize=15,fontfamily="Helvetica"),
                                         legend_height = unit(2.5, "in")),
             column_title = 'Epifactor Group',
             column_title_gp = gpar(fontsize=20,fontfamily='Helvetica'),
             column_title_side = "bottom",
             row_title = 'Cancer Type',
             row_title_gp=gpar(fontsize=20, fontfamily='Helvetica'),
             row_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
             column_names_gp = gpar(fontsize=15, fontfamily='Helvetica'),
             column_names_rot = 60,
             col=col_fun)
# add more titles to heatmap
ht <- draw(h,
           column_title=paste0("Top ",a_title," Genes Permutation Odds Ratio P-value"), column_title_gp=gpar(fontsize=25, fontfamily='Helvetica')
)
ht
dev.off()
write.csv(perm_pval_df[row_order(ht),column_order(ht)],paste0('../../Summary_tables/',project,'_',analysis,'_top_genes_epigenetic_groups_perm_overlap_OR_BH.csv'))
