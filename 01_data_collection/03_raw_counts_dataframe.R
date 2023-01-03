library(argparse)

# create parser object to add command line arguments
parser <- ArgumentParser()
# specify desired options
# by default ArgumentParser will add a help option
parser$add_argument("project",nargs=1, help="project name (TCGA or TARGET)")
parser$add_argument("cancer", nargs=1, help="type of cancer (e.g. BRCA)")
# parse the arguments
args <- parser$parse_args()
print(args)
project <- args$project
cancer <- args$cancer


setwd(paste0('../../',project,'_',cancer,'/01_data_collection/raw_data'))

# number of raw counts files
print(paste0('Number of patients: ',length(list.files())))

# read in data frame that matches sample to file name
samples <- read.delim(paste0("../purity_gdc_sample_sheet_",cancer,".txt"),header=T,sep="\t")
# set rownames
rownames(samples) <- gsub('.gz','',samples$File.Name)
# initialize list of raw counts files
files <- list.files()

# write all samples' counts into one dataframe
# initialize counter
counter <- 0

print('Curating patient raw counts dataframe')
# this function writes the counts of one sample into a larger counts dataframe of all samples
get_counts <- function(filename,samples){
  # filename is the name of the file
  # samples is the sample_sheet datframe that matches sample to filename
  # No output, the function only writes to an existing dataframe
  
  # increment counter
  counter <<- counter + 1
  if (counter%%25==1) {
    print(paste0('Iteration: ',counter))
  }
  # read in the sample's count file
  f <- read.delim(filename, header=F, sep='\t')
  # set count file column to name of the sample
  colnames(f)[2] <- as.character(samples[filename,'Sample.ID'])
  dim(f)
  # if counter is 1 (aka first sample being read), initialize a GLOBAL counts dataframe
  # global is important because we want to the dataframe to exist outside the scope of the function
  if (counter==1){
    counts <<- f # '<<-' means global assignment
  }
  else {
    counts <<- merge(counts,f,by=1,sort=F)
  }
  return()
}

# apply get_counts function to every file in the files object
l <- lapply(files,FUN=get_counts,samples=samples)

# eliminate last 5 rows
counts2 <- counts[1:(nrow(counts)-5),-1]
rownames(counts2) <- counts[1:(nrow(counts)-5),1]


print('Writing raw counts data frame')
# make rnaseq counts directory
dir.create("../../RNA-seq_datasets")
# write raw counts table
write.csv(counts2,file=paste0("../../RNA-seq_datasets/",cancer,"_raw_counts.csv"),quote=F)
