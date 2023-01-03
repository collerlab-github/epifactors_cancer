#****************************************************************************************
#
#                       Gene co-Expression Network Analysis Package
#                               - GO Term Enrichment Analyis
#
#        Objective: 
#                   Automatic construction of gene co-expression networks &
#                   Decomposition of a network into tightly co-expressed modules
#
#        Input:     DNA microarray data file
#
#        Output:    Topological overlap matrix (TOM) and gene co-expressed modules
#
#        Authors:   Bin Zhang                    
#
#                   Array Data Analysis Group (ADAG)
#                   Departments of Human Genetics and Biostatistics
#                   University of California at Los Angeles                
#
#        Contact:   binzhang@mednet.ucla.edu
#
#
#
#        Date:      May 25, 2005
#
#        
#        REV            DATE            BY           DESCRIPTION
#
#        
#****************************************************************************************

#----------------------------------------------------------------------------------------
#                            Description of major variables
#
# INPUT:
#
#     minputfname   ~ module assignment table with the last column as module label
#     geneSymbolCol ~ which column contains gene symbols to be maSharedhed to GO term tables
#     ontologyfname ~ GO term database
#
# Output:
#     GO terms enriched in individual modules
#
#----------------------------------------------------------------------------------------
rm(list=ls(all=TRUE))
library(class)
library(WGCNA)
library(cluster)
library(survival)
library(rpart)
library(lattice) # require is design for use inside functions 
library(Hmisc)

############################## Input #########################################
# 1) set up location of the file subject to GO term enrichment analysis:
setwd("/u/home/j/jennytch/project-xyang123/SCAD_dseq_WGCNA/") # USE HOFFMAN2

# Source function file
source("R-tomfunctions.R")
collect_garbage()
#memory.size(TRUE)   # check the maximum memory that can be allocated
# increase the available memory to 4GB
#memory.limit(size=4095)

# 2) module information file (tab-delimited text file with the last column as module assignment)
################################## IT HAS TO BE "Symbol \t module" CASE SENSITIVE #######
##################################
##### switch MODULE GENE -> Symbol module ###########
## purpose of this is to feed into ontology script ##

#a = read.table("GTEX_lung_combined_080417_less1000.txt", header=T, sep='\t', check.names=F)
#a = read.table("combined_coex_canon_v7_logtransformed.txt", header=T, sep='\t', check.names=F)
#a = read.table("multiscale_significant.modules.split.txt", header=T, sep='\t', check.names=F)
#a = read.table("plos_gtexv7_eqtl.txt", header=T, sep='\t', check.names=F)
a = read.table("megena_modgen_dseq.txt", header=T, sep='\t', check.names=F) #input megena output file 

b = a[,c("GENE", "MODULE")]
colnames(b)= c("Symbol", "module")
write.table(b, file="megena.coexpr.modules.txt", sep = "\t", row.names = F, col.names = T, quote = F)
#minputfname= "MEGENAmodulesLiverMale.txt"  ## in the 1st column we have the gene symbols, in the 2nd col we have the module names
minputfname= "megena.coexpr.modules.txt"
##################################
##################################
# which column includes gene symbol
geneSymbolCol = 1 ## we put the gene symbols in the 1st column!!
# 3) gene ontology database for human V3 chip
#ontologyfname = "/Users/jasonhong/Downloads/MSigDB_Canonical_Pathways.txt" ## this includes reactome, kegg, biocarta pathways
#ontologyfname = "/Users/jasonhong/Box Sync/Mergeomics/GTEx_V7/WGCNA/MSigDB_Canonical_Pathways.txt"
ontologyfname = "msigdb.pathway+gobp.txt"

# specify where to hold the analysis results
##################################
##################################
OutputDirectory= "~/"  
##################################
##################################
############################## END of Input #########################################
############################## Output #########################################
if(OutputDirectory!="") {
 dir.create(OutputDirectory)
}
#######+++++++Parameters may be subject to change+++++++########
# signifLevel: the pvalue cutoff for enriched GO terms
# ctopNrows:   top N most significantly enriched GO terms in each module will be compiled into a single file
# background:  total number of genes on the chip; if you don't know the number, you can simply set it as 0
#  then, the program will use the annotated genes as background
#
OntologyAnalysisDull(inputfname=minputfname, minPopHits=3, ## min number of the gene overlap is 3
                     identifier_col=geneSymbolCol,gene_symbolcol=geneSymbolCol,
                     outdir=OutputDirectory, 
                     ontologyfnlist=c(ontologyfname), useEASEannotation=T, 
                     signifLevel=0.005, ctopNrows=5, background=0) ## show top 5 GO terms with P<5e-3

############################## END of Output #########################################
