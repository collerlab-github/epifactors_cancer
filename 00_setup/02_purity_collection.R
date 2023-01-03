
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# if ('TCGAbiolinks' %in% installed.packages() ==F)
#   BiocManager::install("TCGAbiolinks")
# library(TCGAbiolinks)
library(readxl)
library(stringr)

# 1. read in CPE and CHAT purity tables
# 2. read in cancer patient list
# 3. filter cancer patients for CPE and CHAT cutoffs
# 4. decide on cutoff and save patient list

### only need to manipulate these ###
# cancer
cancer <- 'ACC'
# name of cpe file
cpe_file <- "../../metadata/TCGA_Tumor_Purity_estimates.xlsx"
# name of chat file
chat_file <- paste0("../../metadata/AGPall/AGP-",str_to_lower(cancer),".txt")
# supplementary patient file
sup_file <- '../supplemental/TCGA_ACC_supp_3.xlsx'
# sheet name in sup_file
sheet_name <- 'Subtypes'
# range for excel patients
cell_range <- anchored("A7",dim=c(NA,1))
# best purity cutoff (default NA) ***DECIDE AFTER YOU LOOK AT THE PURITY CUTOFF TABLE**
cutoff <- NA
######################################

# read in cpe purity data
cpe_df <- read_excel(cpe_file,sheet = "Supp Data 1", range = "A4:G9368",col_names=T)
cpe_df <- cpe_df[which(cpe_df$`Cancer type`==cancer),]
class(cpe_df) # data frame
dim(cpe_df) # 9364  7
cpe <- as.data.frame(cpe_df[,c(1,2,7)])
rownames(cpe) <- cpe$`Sample ID`
head(cpe)

# read in CHAT purity data
chat_df <- read.table(chat_file, sep='\t', header=T)
class(chat_df) # data frame
dim(chat_df) # 87 11
chat <- chat_df[,c(1,2)]
chat$sampleid <- gsub('\\.','-',chat$sampleid)
rownames(chat) <- chat$sampleid
head(chat)

# read in cancer patient list
cancer_df <- read_excel(sup_file,sheet=sheet_name,range=cell_range)[,1,drop=T]
cancer_df <- paste0(cancer_df,'-01A')

# find patients in common among cpe, chat, and tcga datasets
cancer_patients <- Reduce(intersect, list(cancer_df,rownames(chat),rownames(cpe)))

# merge purity values into 1 df
# stringsAsFactors = False keeps character vector (string) entries as characters
purity_df <- data.frame(row.names=cancer_patients, CPE = cpe[cancer_patients,3], CHAT = chat[cancer_patients,2], stringsAsFactors = FALSE)

# analyze purity cutoffs
purity_df[which(purity_df$CPE=="NaN"),1] <- "0"
purity_cutoffs <- c()
purity_cutoffs <- append(purity_cutoffs,nrow(purity_df)) 
purity_cutoffs <- append(purity_cutoffs,nrow(purity_df[which((purity_df$CPE>=0.6)|purity_df$CHAT>=0.6),,drop=F]))
purity_cutoffs <- append(purity_cutoffs,nrow(purity_df[which((purity_df$CPE>=0.7)|purity_df$CHAT>=0.7),,drop=F]))
purity_cutoffs <- append(purity_cutoffs,nrow(purity_df[which((purity_df$CPE>=0.8)|purity_df$CHAT>=0.8),,drop=F]))
purity_cutoffs <- append(purity_cutoffs,nrow(purity_df[which((purity_df$CPE>=0.85)|purity_df$CHAT>=0.85),,drop=F]))
purity_cutoffs <- append(purity_cutoffs,nrow(purity_df[which((purity_df$CPE>=0.9)|purity_df$CHAT>=0.9),,drop=F]))
purity_cutoffs <- t(data.frame(purity_cutoffs))
rownames(purity_cutoffs) <- cancer
colnames(purity_cutoffs) <- c(0,0.6,0.7,0.8,0.85,0.9)

# save purity cutoffs
write.csv(purity_cutoffs, file='purity_cutoffs.csv', quote = F)
if(!is.na(cutoff)){
  # save purity values
  write.csv(purity_df[which((purity_df$CPE>=cutoff)|purity_df$CHAT>=cutoff),,drop=F], file=paste0(cancer,'_purity_patients_',cutoff,'.csv'),quote=F)
}

