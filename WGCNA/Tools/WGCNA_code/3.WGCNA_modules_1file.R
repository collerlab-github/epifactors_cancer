###     Collect the WGCNA modules into one file ###
###################################################
intModules=names(table(moduleColors))
all.genes = colnames(datExpr)
pathways <- c()
for (module in intModules){
  mod.name <- paste("WGCNA_", module, sep='')
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their genes
  genes = all.genes[modGenes];
  add.this <- as.character(rep(mod.name, length(genes)))
  add.this <- cbind(add.this, genes)
  if (module != 0)
    pathways <- rbind(pathways, add.this)
}
colnames(pathways) <- c("MODULE", "GENE")
write.table(pathways, file=paste0("Wmodules", ".txt"), sep = "\t", row.names = F, col.names = T, quote = F)

