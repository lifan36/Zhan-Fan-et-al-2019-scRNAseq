library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("~/Desktop/data_analysis/elife/D2/all_0.1_local")

LG180_integrated <- readRDS("elife_microglial_cells_only.rds")
#for figure 3

MG_homeostasis <- c("Trem2", "Tmem119","P2ry12","P2ry13", "Selplg", "Aif1", "Cx3cr1", "Csf1r", "Hexb")
DoHeatmap(LG180_integrated, features = MG_homeostasis, size = 2, draw.lines = T) 

#Figure 3 b-g
Cluster_1_2 <- subset(LG180_integrated, idents = c("1", "2"))
VlnPlot(Cluster_1_2, features = c("Tmem119","P2ry12","P2ry13", "Selplg", "Cx3cr1", "Csf1r"), group.by = "Condition", pt.size = 0, combine = T)
