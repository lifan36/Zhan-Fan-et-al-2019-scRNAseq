library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(patchwork)

setwd("")

LG180_integrated <- readRDS("elife_microglial_cells_only.rds")

#Figure 3 a-f
Cluster_2 <- subset(LG180_integrated, idents = "2")
VlnPlot(Cluster_1_2, features = c("Tmem119","P2ry12","P2ry13", "Selplg", "Cx3cr1", "Csf1r"), group.by = "Condition", pt.size = 0, combine = T)& 
  theme( plot.title = element_text( face = "italic") )
