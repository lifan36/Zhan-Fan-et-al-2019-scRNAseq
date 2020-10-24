library(ggplot2)
library(Seurat)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
library(patchwork)
# Figure 6a and 6b plots
# Subset Mac2+ cells based on Mac2 expression level

# Load in Seurat object
setwd("")
data<-readRDS("elife_microglial_cells_only.rds")
data@active.ident
# Extract raw data and subset "Lgals3" expression data
raw <- data@assays$RNA@scale.data

lgals3.data <- raw[rownames(raw)=="Lgals3"]
names(lgals3.data) <- colnames(raw)

data <- AddMetaData(
  object = data,
  metadata = lgals3.data,
  col.name = "lgals3.data")
head(data@meta.data)

# Mac2+ cells considered those with > 1 SD Mac2+ expression
lgals3.cutoff <- mean(lgals3.data) + 1*sd(lgals3.data)
table(data@meta.data$lgals3.data>lgals3.cutoff)

cell_index <- which(data@meta.data$lgals3.data>lgals3.cutoff) 
cell_name <- rownames(data@meta.data)[cell_index]

# add a level to ident 
data$mac2 <- "no"
data$mac2[row.names(data@meta.data) %in% cell_name] <- "mac2"
data$cluster_mac2 <- as.numeric(data$seurat_clusters)
data$cluster_mac2[row.names(data@meta.data) %in% cell_name] <- "mac2"
data$seurat_clusters_new <- as.numeric(data$seurat_clusters)

# UMAP and visualizations of Mac2+ cells
Idents(data)<-"mac2"
DimPlot(data,cols=c("grey","orange"))
ggsave("umap mac2+.pdf", plot = last_plot(), device = "pdf", path = "",
       scale = 0.8, width = 5, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)

RidgePlot(data, features=c("Lgals3"), col=c("grey","orange"))
ggsave("ridge plot mac2+ microglia.pdf", plot = last_plot(), device = "pdf", path = "",
       scale = 0.8, width = 5.5, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)

