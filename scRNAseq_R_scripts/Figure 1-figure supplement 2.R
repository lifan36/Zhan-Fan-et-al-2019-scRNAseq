library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("~/Desktop/data_analysis/elife/D2/all_0.1_local")

LG180_integrated <- readRDS("~/Desktop/data_analysis/elife/D2/elife_Ctrl_D0_D2.rds")

#Figure S2a-c
VlnPlot(object = LG180_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample_name", ncol = 3, pt.size=0, idents=NULL)

#Figure S2d
FeatureScatter(object = LG180_integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample_name", pt.size=0.5)

#Figure S2e
Idents(LG180_integrated) <- "sample_name"
DimPlot(LG180_integrated, reduction = "pca")

#Figure S2f
ElbowPlot(LG180_integrated, ndims = 50)

