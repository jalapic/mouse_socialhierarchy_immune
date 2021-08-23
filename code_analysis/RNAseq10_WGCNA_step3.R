# continue from step2.R



# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
# dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = my_power);
# saveRDS(dissTOM,glue("results_RNAseqRDS/dissTom_{my_tissue}.RDS"))

dissTOM <- readRDS(glue("results_RNAseqRDS/dissTom_{my_tissue}.RDS"))


restGenes= (moduleColors != "grey")
dissTOM = 1-TOMsimilarityFromExpr(datExpr[,restGenes],
                                  power = my_power,
                                  networkType = "signed hybrid");
plotTOM = dissTOM^7;
# look at the network w/out "grey" modules
colnames(dissTOM) = rownames(dissTOM) = colnames(datExpr[restGenes])
hier1=flashClust(as.dist(dissTOM), method = "average" )



plotDendroAndColors(hier1, colors = data.frame(moduleColors[restGenes]),
                    c("Module Tree Cut"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05, main = "Gene dendrogram and module colors")

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
diag(dissTOM) = NA






png(glue("results_figures/TOMplot_{my_tissue}_{my_networkType}_{my_power}.png"),   width = 12, height = 12, units = "cm", res = 600)

TOMplot(dissim = 1-dissTOM^7, 
        hier1, 
        as.character(moduleColors[restGenes]), 
        min = "Network heatmap plot, removed unassigned genes")
dev.off()



# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, status))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,9);
par(cex = 0.9)

png(glue("results_figures/dendoAndcolors_{my_tissue}.png"),
    width = 12, height = 18, units = "cm", res = 600)

plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

dev.off()
