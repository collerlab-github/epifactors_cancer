setwd("~/Desktop/SCAD_Dseq/MEGENA_hv")


library(readxl)
coexp_modules <- read_excel("ontology_output.xlsx",col_names=TRUE)

colnames(coexp_modules)

# write info file (module, source, description)
dseq_megena_info <- coexp_modules[,c(1,3,4)]
# msd$MODULE <- gsub(pattern="WGCNA_",replacement="",x=msd$MODULE) # remove WGCNA_
#write.table(msd, file="dseq_megena_info.txt", sep='\t', row.names=F, col.names=T, quote=F)

#concatenate functional annotations of pathway colors
library(dplyr)
#dseq_megena_info <- read.delim("dseq_megena_info.txt", header=TRUE)

filter(.data=dseq_megena_info, MODULE=='c1_7', SOURCE=='canonical')
# find unique module and source combos
dseq_groups <- dseq_megena_info[!duplicated(dseq_megena_info[,c(1,2)]),c(1,2)]
# combine into one column
dseq_groups$combined <- paste(dseq_groups$MODULE,dseq_groups$SOURCE)

# loop through each module-source combo in the combined column
# make a result dataframe for the results
result_dseq <- data.frame(MODULE=character(), SOURCE=character(), DESCR=character())
for (modsource in dseq_groups$combined){
  # what we know: modsource
  # first, identify the module and the source
  mod <- filter(dseq_groups, combined==modsource)$MODULE
  source <- filter(dseq_groups, combined==modsource)$SOURCE
  
  # second, get a list of pathway descriptions corresponding to the mod and source
  descriptions <- filter(dseq_megena_info, MODULE==mod, SOURCE==source)$DESCR
  
  # third, combine the pathway descriptions using paste
  combined_descriptions <- paste(descriptions, collapse='; ')
  
  # fourth, add new entry into result_dseq dataframe
  result_list <- setNames(as.list(c(mod, source, combined_descriptions)),c('MODULE','SOURCE','DESCR'))
  result_dseq <- rbind(result_dseq,result_list)
}

write.table(result_dseq, file="dseq_megena_modules_combineDescr_hv.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

###########################3
# write module and gene file. use individual excel files
library(stringr)
# store list of module files
color_files <- list.files(path="dseq_all_WGCNA/")
# create big dataframe
big_modgen <- data.frame(MODULE=character(), GENE=character())

# function to load module genes and combine it with big dataframe
create_modgen <- function(color_file){
  modgen <- read.csv(paste0("dseq_all_WGCNA/",color_file), header=T, row.names=1) # read in 1 file 
  # get color name (remove ME and .csv from name)
  color_name <- str_sub(string=color_file, start=3, end=str_length(color_file)-4)
  # add color column to modgen
  modgen$MODULE <- color_name
  colnames(modgen)[1] <- 'GENE'
  
  # reorder columns
  modgen <- modgen[, c(2,1)]
  
  # add modgen to big_modgen in global scope
  big_modgen <<- rbind(big_modgen, modgen)
  return('file added')
}


# combine genes of all 71 modules into 1 dataframe
sapply(X=color_files, FUN=create_modgen) # (running 71x)


write.table(big_modgen, file="dseq_all_WGCNA_modules.txt", sep='\t', row.names=F, col.names=T, quote=F)
