setwd("~/Desktop/SCAD_Dseq/MEGENA/MEGENA_output.2")

modules <- read.csv("modules_dseq.csv",header=TRUE,row.names=1)
modules[1:6,1:6]

modules.annot <- data.frame(MODULE=character(),GENE=character())

# what we want to do for each module
# step 1: add module name to the MODULE column of modules.annot for as many rows as there are genes
# step 2: add gene name to the GENE column of modules.annot for each gene in the module (not including NA)
    # make a mini data frame that combines step 1 and step 2
    # rbind() to combine modules.annot with the mini dataframe


structure <- function(mod_name){
  # step 1: add module name to the MODULE column of modules.annot for as many rows as there are genes
  gene_list <- modules[mod_name,] # list of genes
  gene_list_isna <- is.na(gene_list) # boolean list of whether entry is NA
  final_gene_list <- gene_list[,gene_list_isna==F] # final list of only genes
  final_gene_list <- as.character(final_gene_list) # transform from class dataframe to class character
  
  # create data frame with 2 columns (module, gene). repeat module name for length of final gene list (NAs not present)
  mini_mod.annot <- data.frame(MODULE=rep(mod_name,length(final_gene_list)),GENE=final_gene_list)
  
  # add specific data frame to overall (step 2)
  # <<- assigns a global variable so that the change to the datframe applies outside the function's scope
  modules.annot <<- rbind(modules.annot, mini_mod.annot)
  }

# iterate through the module names and apply the structure function (this can be done with sapply)
# note: in sapply, the first argument of the applied function is the iterated element
# in this case: for iteration i,  mod_name is the ith element in rownames(modules)
sapply(rownames(modules),FUN=structure)

write.table(modules.annot, file="megena_modgen_dseq.txt",quote=FALSE,sep='\t',row.names=F,col.names=T)

