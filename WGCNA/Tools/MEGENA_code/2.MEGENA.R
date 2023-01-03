library(MEGENA)
#patientData <- read.csv("mapped_normcount.csv",header = TRUE);
#patientData = unique(patientData)
# datExpr0 = as.data.frame((patientData));
# datExpr0 = as.data.frame(t(patientData[, -c(1:1)]));
# names(datExpr0) = patientData$X;
# rownames(datExpr0) = names(patientData[,-c(1:1)]);
# datExpr0 = as.data.frame(t(datExpr0))
# save(datExpr0, file="datExpr0.RData")
 
lnames <-load(file="datExpr1.RData")
lnames

# input parameters
method = "spearman"
FDR = 0.05 
cor.pval = 0.01
n.perm = 100
hub.perm = 100
doPar = FALSE
num.cores = 4
mod.pval = 0.05
hub.pval = 0.05
min.size = 10
max.size = 2500
do.MSigDB = FALSE
do.survival = FALSE
do.DEG.subtype = FALSE
do.DEG.TumorNormal = FALSE 
TCGA.barcode = TRUE
plot.modules = TRUE
annot.table=NULL

# find the corr between gene pairs:
sig.corList <- calculate.correlation(as.matrix(datExpr1), doPerm = n.perm, 
                                     doPar = TRUE, num.cores = num.cores, FDR.cutoff = FDR, 
                                     n.increment = 100, is.signed = FALSE, output.permFDR = TRUE, 
                                     output.corTable = TRUE, saveto = NULL, method = method)
save(sig.corList, file="sig.corListLVRml.ing.RData")

system.time(PFN.ijw <- calculate.PFN(edgelist = sig.corList, doPar = TRUE, num.cores = num.cores))
colnames(PFN.ijw) <- c("TAIL", "HEAD", "WEIGHT")
net.file <- "network_dseq.RData" 
save(PFN.ijw, file="net.file")
# save the PFN - resulting network - into following file
net.file <- "network_lvr_ml.ing.txt"
write.table(PFN.ijw, file = net.file, sep = "\t", row.names = F, col.names = T, quote = F)
colnames(PFN.ijw) <- c("row", "col", "weight")

g <- graph.data.frame(PFN.ijw,directed = FALSE)

##### perform MCA clustering.
MEGENA.output <- do.MEGENA(g,
                           mod.pval = mod.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                           min.size = 50,max.size = vcount(g)/2,
                           doPar = doPar,num.cores = num.cores,n.perm = hub.perm,
                           save.output = FALSE)


summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pval = mod.pval,hub.pval = hub.pval,
                                       min.size = 50,max.size = vcount(g)/2,
                                       annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                       output.sig = TRUE)
write.csv(summary.output$module.table, "moduleTable_ing.csv")

modules <- matrix(NA, nrow = 130, ncol = 2884) #nrow=number of modules, ncol=length of the module that has most genes
modules <- as.data.frame(modules)

for(i in 1:length(summary.output$modules)){
  rownames(modules)[i] <- names(summary.output$modules)[i]
  modules[i,1:length(summary.output$modules[[i]])] <- summary.output$modules[[i]]
}

write.csv(modules, "modules_ing.csv")
modules <- as.data.frame(t(modules), stringsAsFactors = F)
for_GO <- gather(modules)
write.csv(for_GO, "modules_for_GO_ing.csv")

