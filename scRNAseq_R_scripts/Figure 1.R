library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)
library(patchwork)
setwd("")

LG180_integrated <- readRDS("elife_Ctrl_D0_D2.rds")

DefaultAssay(LG180_integrated) <- 'RNA'
LG180_integrated <- ScaleData(LG180_integrated, features = rownames(LG180_integrated))

#Fig. 1c
#top20 HeatMap
LG180_markers <- FindAllMarkers(LG180_integrated, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
top20 <- LG180_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
pdf("HeatMapTop20_all_0.1.pdf", width=16, height=8)
DoHeatmap(LG180_integrated, features = top20$gene) + NoLegend()
dev.off()

#Fig. 1d
pdf("UMAP_all_0.1.pdf", width=8, height=6)
DimPlot(LG180_integrated, reduction = 'umap', label = T)
dev.off()

#Fig. 1e
VlnPlot(object = LG180_integrated, features = "Itgam", pt.size = 0) & 
  theme( plot.title = element_text( face = "italic") )

VlnPlot(object = LG180_integrated, features = "Aif1", pt.size = 0) & 
  theme( plot.title = element_text( face = "italic") )

VlnPlot(object = LG180_integrated, features = "Cx3cr1", pt.size = 0) & 
  theme( plot.title = element_text( face = "italic") )

VlnPlot(object = LG180_integrated, features = "Csf1r", pt.size = 0) & 
  theme( plot.title = element_text( face = "italic") )

#Remove clusters 6, 7, 8, 9 with low Itgam/CD11d, Aif1/Iba1, Cx3cr1 and Csf1r expression
table(LG180_integrated@active.ident)
LG180_integrated <- subset(LG180_integrated, idents = c("6","7", "8", "9"), invert = TRUE)
table(LG180_integrated@active.ident)

####Reclustering
LG180_integrated <- FindVariableFeatures(LG180_integrated, selection.method = "vst", nfeatures = 2000)
LG180_integrated <- ScaleData(object = LG180_integrated)
LG180_integrated <- RunPCA(object = LG180_integrated, features = rownames(x = LG180_integrated), verbose = FALSE)
ElbowPlot(LG180_integrated)
LG180_integrated <- FindNeighbors(LG180_integrated, dims = 1:15)
LG180_integrated <- FindClusters(LG180_integrated, resolution = 0.075)
LG180_integrated <- RunUMAP(LG180_integrated, dims = 1:15)

DimPlot(LG180_integrated, reduction = "umap", label = T)

# rename clusters from 0-5 to 1-6.
n <- dim(table(LG180_integrated@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
LG180_integrated@active.ident <- plyr::mapvalues(x = LG180_integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)

LG180_integrated@active.ident <- factor(LG180_integrated@active.ident, levels=1:n)

DimPlot(LG180_integrated, reduction = 'umap', label = T)

DefaultAssay(LG180_integrated) <- 'RNA'
VlnPlot(object = LG180_integrated, features = c("Itgam"), pt.size = 0)

#remove cluster 6 because of lack of Itgam/Cd11b
LG180_integrated <- subset(LG180_integrated, idents = c("6"), invert = TRUE)
#Fig. 1i
DimPlot(LG180_integrated, reduction = 'umap', label = T)
#Fig. 1j
DimPlot(LG180_integrated, reduction = 'umap', label = T, split.by = "sample_name", ncol = 3)

DefaultAssay(LG180_integrated) <- 'RNA'
LG180_integrated <- ScaleData(LG180_integrated, features = rownames(LG180_integrated))
saveRDS(LG180_integrated, file = 'elife_microglial_cells_only.rds')

