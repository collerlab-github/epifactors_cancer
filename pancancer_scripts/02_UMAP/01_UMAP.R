library(umap)
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(dplyr)
library(stringr)

project <- 'TCGA'
dir.create(paste0('../../',project,'_PAN/02_UMAP'),showWarnings = F)

setwd(paste0('../../',project,'_PAN/02_UMAP'))

set.seed(123)
# load epigene data
filename <- list.files('../RNA-seq_datasets',pattern='log2norm_counts.csv',full.names = T)
epigenes <- read.csv(filename,row.names=1)
# patients <- gsub('(.*)\\w$','\\1',colnames(epigenes))
patients <- gsub('\\.','-',colnames(epigenes))
epigenes <- rownames(epigenes)

# load batch corrected pancancer counts (XENA)
filename <- list.files('../RNA-seq_datasets',pattern='rsem',full.names=T)
exp_df <- read.table(filename)
colnames(exp_df) <- exp_df[1,]
names(patients) <- gsub('(.*)\\w$','\\1',patients)
# colnames(exp_df) <- paste0(gsub("-","\\.", exp_df[1,]))
exp_df_filtered <- exp_df[which(exp_df$sample %in% epigenes), which(colnames(exp_df) %in% names(patients))]
colnames(exp_df_filtered) <- patients[colnames(exp_df_filtered)]
rnames <- exp_df[rownames(exp_df_filtered),'sample']
exp_df_filtered <- apply(exp_df_filtered, 2, as.numeric)
rownames(exp_df_filtered) <- rnames
exp_df_filtered <- t(exp_df_filtered)
exp_df_filtered <- na.omit(exp_df_filtered)
exp_df_filtered <- exp_df_filtered[!is.na(rownames(exp_df_filtered)),]
write.csv(exp_df_filtered, '../RNA-seq_datasets/filtered_pancan.csv')

# rm(exp_df)

# load patient info
filename <- list.files('../01_filtering',pattern='patient_data.csv', full.names=T)
patient_data <- read.csv(filename,row.names=1)
rownames(patient_data) <- gsub('\\.','-',rownames(patient_data))
patient_data <- na.omit(patient_data[rownames(exp_df_filtered),,drop=F])

# umap
epigenes.umap <- umap(as.data.frame(exp_df_filtered),random_state=123)

umap_coords <- data.frame(epigenes.umap$layout)
umap_coords$cancer <- patient_data[rownames(umap_coords),1]
umap_coords <- na.omit(umap_coords)
umap_coords$label <- ifelse(duplicated(umap_coords$cancer),'No',umap_coords$cancer)

ggplot(umap_coords, aes(x=X1,y=X2,label=cancer))+
  geom_point(size=1,aes(color=cancer),alpha=0.4)+
  guides(color = guide_legend(override.aes = list(size=5, alpha=0.7), title='Cancer'))+
  theme_bw()+
  geom_text_repel(data = subset(umap_coords, !(label == "No")),box.padding=0.5,size=3)+
  geom_text_repel(data = subset(umap_coords%>%filter(X1>5,X1<6,X2>-5,X2< -2.5),!duplicated(umap_coords%>%filter(X1>5,X1<6,X2>-5,X2< -2.5)%>%select(cancer))),box.padding=0.5,size=3)+
  geom_text_repel(data = subset(umap_coords%>%filter(X1>2,X1<2.5,X2>-6,X2< -5),!duplicated(umap_coords%>%filter(X1>2,X1<2.5,X2>-6,X2< -5)%>%select(cancer))),box.padding=0.5,size=3)+
  geom_text_repel(data = subset(umap_coords%>%filter(X1>-7,X1< -5,X2>-7,X2< -5),!duplicated(umap_coords%>%filter(X1>-7,X1< -5,X2>-7,X2< -5)%>%select(cancer))),box.padding=0.5,size=3)+
  # geom_text_repel(data = subset(umap_coords%>%filter(X1>-5,X1< 0,X2>0,X2< 4),!duplicated(umap_coords%>%filter(X1>-5,X1< 0,X2>0,X2< 4)%>%select(cancer))),box.padding=0.5,size=3)+
  xlab('UMAP1')+
  ylab('UMAP2')+
  # geom_text_repel(data = subset(umap_coords%>%filter(X1>-12.5,X1< -10,X2>-10,X2< -5),!duplicated(umap_coords%>%filter(X1>-12.5,X1< -10,X2>-10,X2< -5)%>%select(cancer))),box.padding=0.5,size=3)+
  # geom_text_repel(data = subset(umap_coords%>%filter(X1>-10,X1< -5,X2>-1,X2< 5),!duplicated(umap_coords%>%filter(X1>-10,X1< -5,X2>-1,X2< 5)%>%select(cancer))),box.padding=0.5,size=3)+
  # geom_text_repel(data = subset(umap_coords%>%filter(X1>-3,X1< 0,X2>-1,X2< 2.5),!duplicated(umap_coords%>%filter(X1>-3,X1< 0,X2>-1,X2< 2.5)%>%select(cancer))),box.padding=0.15,size=3)+
  # geom_text_repel(data = subset(umap_coords%>%filter(X1>-3,X1< 0,X2>3,X2< 10),!duplicated(umap_coords%>%filter(X1>-3,X1< 0,X2>3,X2< 10)%>%select(cancer))),box.padding=0.5,size=3)+
  # geom_text_repel(data = subset(umap_coords%>%filter(X1>0,X1< 6,X2>-1,X2< 6),!duplicated(umap_coords%>%filter(X1>0,X1< 6,X2>1,X2< 6)%>%select(cancer))),box.padding=0.5,size=3)+
  # 
  theme(text=element_text(family="Helvetica", size=10),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10))+
  ggtitle('UMAP of All Cancer Epigenes')
  
ggsave('umap_bc_all_epigenes_final.png',width=7,height=7)

# umap of top nmf epigenes
nmf_genes <- c()
cancers <- na.omit(unique(patient_data$type))
# extract top nmf genes
for (cancer in cancers){
  top_genes <- c(rownames(read.csv(paste0('../../',project,'_',cancer,'/03_nmf/Rank_2/nmf_lee_rank2_feature1_genes.csv'),row.names=1)),
                 rownames(read.csv(paste0('../../',project,'_',cancer,'/03_nmf/Rank_2/nmf_lee_rank2_feature1_genes.csv'),row.names=1)))
  nmf_genes <- unique(c(nmf_genes,top_genes))
}
# subset top nmf epigenes
nmf_epigenes <- exp_df_filtered[,nmf_genes[which(nmf_genes%in%colnames(exp_df_filtered))]]
# umap
nmf_epigenes.umap <- umap(nmf_epigenes,random_state=123)
nmf_umap_coords <- data.frame(nmf_epigenes.umap$layout)
nmf_umap_coords$Cancer <- patient_data[rownames(nmf_umap_coords),1]
nmf_umap_coords <- na.omit(nmf_umap_coords)
nmf_umap_coords$label <- ifelse(duplicated(nmf_umap_coords$Cancer),'No',nmf_umap_coords$Cancer)

# png('umap_top_nmf_epigenes_0.5_legend.png',width=1200,height=1200, res=300, units = 'px')
ggplot(nmf_umap_coords, aes(x=X1,y=X2,label=Cancer))+
  geom_point(size=0.5,aes(fill=Cancer, color=Cancer),alpha=0.4)+
  geom_text_repel(data = subset(nmf_umap_coords, !(label == "No")),box.padding=0.5,max.overlaps = 20)+
  # geom_text_repel(data = subset(nmf_umap_coords%>%filter(X1>-10,X1< -9,X2>-7,X2< -5),!duplicated(nmf_umap_coords%>%filter(X1>-10,X1< -9,X2>-7,X2< -5)%>%select(Cancer))),box.padding=0.5)+
  geom_text_repel(data = subset(nmf_umap_coords%>%filter(X1>-9,X1< -8,X2>0,X2< 2.5),!duplicated(nmf_umap_coords%>%filter(X1>-9,X1< -8,X2>0,X2< 2.5)%>%select(Cancer))),box.padding=0.5)+
  theme_bw()+
  xlab('UMAP1')+
  ylab('UMAP2')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank(),
        text=element_text(family="Helvetica",size=30),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), title=element_text(size=15),
        legend.text = element_text(family='Helvetica',size=10), legend.title=element_text(family="Helvetica",size=10),
        axis.text = element_text(size=15))+
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border=element_blank(),
  #       panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position='none')+
  ggtitle('UMAP of Top NMF Cancer Epigenes')+
  guides(color = guide_legend(override.aes=list(size=5,color='black',pch=21,alpha=0.7)))
ggsave('umap_bc_top_nmf_epigenes_0.5.png',width=10,height=7)
# dev.off()
