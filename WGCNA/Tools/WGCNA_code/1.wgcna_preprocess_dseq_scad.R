# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
patientData <- read.csv("mapped_normcount.csv",header = TRUE);
#patientData = unique(patientData)
# Take a quick look at what is in the data set:
dim(patientData);
names(patientData);

#Store the patient diagnosis 
#diagnosis = patientData[1,]
# remove the diagnosis row
#patientData<- patientData[-1,]
# # #=====================================================================================
# # #
# # #  Code chunk 2
# # #
# # #=====================================================================================
# # 
# # 

# call out healthy samples
patientData = patientData[,c(1,31:80)]

  datExpr0 = as.data.frame((patientData));

 
  # rownames(datExpr0) <-patientData[,1]
  # datExpr0 <- datExpr0[,-1]
  # datExpr0 = t(datExpr0)
  # names(datExpr0) = patientData$BIOCHEMICAL #column that has metabolite names
  datExpr0 = as.data.frame(t(patientData[, -c(1:1)]));
  names(datExpr0) = patientData$X;
  rownames(datExpr0) = names(patientData[,-c(1:1)]);
  
#
#
# #=====================================================================================
# #
# #  Code chunk 3
# #
# #=====================================================================================
#
#

 gsg = goodSamplesGenes(datExpr0, verbose = 3);
 gsg$allOK
#
#
# #=====================================================================================
# #
# #  Code chunk 4
# #
# #=====================================================================================
#
#
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
#
#
# #=====================================================================================
# #
# #  Code chunk 5
# #
# #=====================================================================================
#
#
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Gene Clustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#dev.off()


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 1400000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 1400000, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
#names(datExpr) = patientData$BIOCHEMICAL
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================





#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================




#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================

print(names(datExpr))
save(datExpr,datExpr0, file = "dseq_input_scad.RData")

