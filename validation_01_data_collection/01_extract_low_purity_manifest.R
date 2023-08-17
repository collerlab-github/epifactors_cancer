library(argparse)
library(dplyr)
# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("project",nargs=1,help="project name (TCGA or TARGET)")
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
# parse the arguments
args <- parser$parse_args()
print(args)
project <- args$project
cancer <- args$cancer

setwd(paste0('../../',project,'_',cancer,'/01_data_collection'))
print(getwd())
# define filenames

# inputs
purity_patients_file <- paste0(cancer,'_low_purity_patients.csv')
manifest_file <- paste0('gdc_manifest_',cancer,'.txt')
samples_file <- paste0('gdc_sample_sheet_',cancer,'.tsv')
# outputs
purity_sample_sheet_file <- paste0('purity_gdc_sample_sheet_',cancer,'.txt')
purity_manifest_file <- paste0('purity_gdc_manifest_',cancer,'.txt')

# read purity cutoff patients
purity_patients <- read.csv(purity_patients_file,row.names=1,header=T)
print(paste0('Number of patients: ', nrow(purity_patients)))

# read in samples and manifest file
print('Reading manifest and sample data')
manifest <- read.delim(manifest_file,header=T,sep='\t')
samples <- read.delim(samples_file,header=T,sep='\t') %>% filter(Sample.Type=='Primary Tumor')


# filter samples by given patient list, then filter manifest by the filename of the samples
print('Extracting patients from manifest')
purity_samples <- samples[which(samples$Case.ID %in% rownames(purity_patients)),]
purity_samples <- purity_samples[grepl('01A$',purity_samples$Sample.ID),]
purity_samples <- purity_samples[!duplicated(purity_samples$Case.ID),]
# rownames(purity_patients)[which(rownames(purity_patients)%in%samples$Sample.ID==F)]
# remove purity sample duplicates
print('saving outputs')
purity_samples <- purity_samples[duplicated(purity_samples$Sample.ID)==F,]
print(paste0('Number of patients in manifest: ', nrow(purity_samples)))
write.table(purity_samples, file=purity_sample_sheet_file, quote=F,sep='\t',row.names=F)
purity_manifest <- manifest[which(manifest$id%in%purity_samples$File.ID),]
write.table(purity_manifest, file=purity_manifest_file,quote=F,sep="\t",row.names=F)
