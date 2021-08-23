my_tissue = 'Liver'
my_tissue = 'Spleen'

dds <- readRDS(glue("results_RNAseqRDS/dds_{my_tissue}.RDS"))
## Perform a variance-stabilizing transformation
vsd <- getVarianceStabilizedData(dds)
## Many functions expect the matrix to be transposed
datExpr <- t(vsd)
## check rows/cols
nrow(datExpr)
ncol(datExpr)
rownames(datExpr)



Samples = rownames(datExpr);
collectGarbage()
# ===========================================================================

# Run WGCNA for each alpha, beta, and gamma 



WGCNA_part_1 <- function(my_tissue = "Liver"){
  dds <- readRDS(glue("results_RNAseqRDS/dds_{my_tissue}.RDS"))
  ## Perform a variance-stabilizing transformation
  vsd <- getVarianceStabilizedData(dds)
  ## Many functions expect the matrix to be transposed
  datExpr_tissue_temp <- t(vsd)
  nrow(datExpr_tissue_temp)
  rownames(datExpr_tissue_temp)
  saveRDS(datExpr_tissue_temp, glue("results_RNAseqRDS/datExpr_{my_tissue}.RDS"))
  
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr_tissue_temp, powerVector = powers, verbose = 5, )
  # Plot the results:
  print(sft)
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = glue("{my_tissue}: Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  
  abline(h=0.85,col="red")  #changed to 0.8 (originally 0.9 in the tutorial)
  abline(h=0.90,col="blue")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = glue("{my_tissue}:Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  print(sft)
  
}

WGCNA_part_1("Liver") # 5
WGCNA_part_1("Spleen") # 16 WTF

# TOMType = "unsigned", networkType = "unsigned",
# TOMType = "signed", networkType = "signed hybrid",
# TOMType = "signed", networkType = "signed",


WGCNA_get_net <- function(my_tissue = "Liver",
                          my_power =  5, 
                          my_TOMType ="unsigned", 
                          my_networkType = "unsigned"){
  
  x <- readRDS(glue("results_RNAseqRDS/datExpr_{my_tissue}.RDS"))  
  net = blockwiseModules(x, 
                         power = my_power,
                         TOMType = my_TOMType, 
                         networkType = my_networkType,
                         minModuleSize = 30,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.25,
                         numericLabels = TRUE, 
                         pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)
  
  
  saveRDS(net, glue("results_RNAseqRDS/net_{my_networkType}_{my_tissue}_{my_power}.RDS"))
  
  
}

# WGCNA_get_net("Liver", 4, "unsigned", "unsigned")
# WGCNA_get_net("Liver", 5, "unsigned", "unsigned")
# WGCNA_get_net("Liver", 6, "unsigned", "unsigned")
# WGCNA_get_net("Liver", 3, "signed", "signed hybrid")
# WGCNA_get_net("Liver", 4, "signed", "signed hybrid")
WGCNA_get_net("Liver", 5, "signed", "signed hybrid")
# WGCNA_get_net("Liver", 6, "signed", "signed hybrid")
# 
# WGCNA_get_net("Spleen", 12, "unsigned", "unsigned")
# WGCNA_get_net("Spleen", 14, "unsigned", "unsigned")
# WGCNA_get_net("Spleen", 16, "unsigned", "unsigned")
# WGCNA_get_net("Spleen", 18, "unsigned", "unsigned")
# WGCNA_get_net("Spleen", 12, "signed", "signed hybrid")
WGCNA_get_net("Spleen", 14, "signed", "signed hybrid")
# WGCNA_get_net("Spleen", 16, "signed", "signed hybrid")
# WGCNA_get_net("Spleen", 18, "signed", "signed hybrid")

# WGCNA_get_net("Spleen", 1, "signed", "signed hybrid")





# ========================================================================================

files <- list.files(path = "results_RNAseqRDS", pattern = "net_*")
files <- files[9]
for (i in 1:length(files)){
  filename <- files[i]
  
  net <- readRDS(glue('results_RNAseqRDS/{filename}'))
  
  # open a graphics window
  sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  
  # Plot the dendrogram and the module colors underneath
  
  dev.off() # make sure you do this before AND after 
  png(file = glue("results_figures/cluster_dendo_{str_sub(filename,5,-5)}.png"),
      width=600, height=350)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, 
                      main = str_sub(filename, 5, -5))
  
  dev.off()
  
  
}


# moduleLabels = net$colors
# moduleColors = labels2colors(net$colors)
# MEs = net$MEs;
# geneTree = net$dendrograms[[1]];

