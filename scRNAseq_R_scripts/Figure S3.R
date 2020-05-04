library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("~/Desktop/data_analysis/elife/D2/all_0.1_local")

LG180_integrated <- readRDS("~/Desktop/data_analysis/elife/D2/elife_Ctrl_D0_D2.rds")

DefaultAssay(LG180_integrated) <- 'RNA'
LG180_integrated <- ScaleData(LG180_integrated)

#Figure S3a
VlnPlot(object = LG180_integrated, features = "Ly6c1", pt.size = 0.1)
VlnPlot(object = LG180_integrated, features = "Ly6a", pt.size = 0.1)

#Figure S3b
VlnPlot(object = LG180_integrated, features = "Acta2", pt.size = 0.1) 
VlnPlot(object = LG180_integrated, features = "Vtn", pt.size = 0.1)
VlnPlot(object = LG180_integrated, features = "Tagln", pt.size = 0.1) 
VlnPlot(object = LG180_integrated, features = "Myl9", pt.size = 0.1)

#Figure S3c
VlnPlot(object = LG180_integrated, features = "Slc1a2", pt.size = 0.1) 
VlnPlot(object = LG180_integrated, features = "Mt3", pt.size = 0.1)
VlnPlot(object = LG180_integrated, features = "Clu", pt.size = 0.1) 
VlnPlot(object = LG180_integrated, features = "Aldoc", pt.size = 0.1)

