library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("~/Desktop/data_analysis/elife/D2/all_0.1_local")

LG180_integrated <- readRDS("elife_microglial_cells_only.rds")

LG180_markers <- read.csv(file = "LG180_markers.csv", header=T,row.names =1)

top10 <- LG180_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# Function to change x labels of individual violin plots
VlnPlot_2 <- function(object, features, pt.size, ncol, xlab) {
  
  # Main function
  main_function <- function(object = object, features = features,pt.size = pt.size, ncol = ncol, xlab = xlab) {
    VlnPlot(object = object, features = features, pt.size = pt.size, ncol = ncol) + 
      labs(x = xlab)
  }
  
  # Apply main function on all features
  p <- lapply(X = features, object = object, xlab = xlab, pt.size = pt.size, ncol = ncol,
              FUN = main_function)
  
  # Arrange all plots using cowplot
  # Adapted from Seurat
  # https://github.com/satijalab/seurat/blob/master/R/plotting.R#L1100
  # ncol argument adapted from Josh O'Brien
  # https://stackoverflow.com/questions/10706753/how-do-i-arrange-a-variable-list-of-plots-using-grid-arrange
  cowplot::plot_grid(plotlist = p, ncol = ncol)
}

pdf("Cluster_all_DEGs_test.pdf", width=24, height=36)
VlnPlot_2(object = LG180_integrated, features = top10$gene, pt.size = 0, ncol =5, xlab = "Cluster")
dev.off()

#For Cluster-3, skip lncRNAs (Gm42418, Gm26917), AY036118 and the mitochondrially encoded genes (mt-*)

pdf("Cluster3_DEGs.pdf", width=24, height=6)
VlnPlot_2(object = LG180_integrated, features = c("Lars2", "Xist", "Dst","Macf1","Atp5g1","Kcnq1ot1","Rps26","Dbi","Rps8","Cox6c"),
          pt.size = 0, ncol =5, xlab = "Cluster")
dev.off()
