library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)
library(patchwork)
setwd("")

LG180_integrated <- readRDS("~/Desktop/data_analysis/elife/D2/elife_Ctrl_D0_D2.rds")

DefaultAssay(LG180_integrated) <- 'RNA'
LG180_integrated <- ScaleData(LG180_integrated)

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


#Remove clusters 6, 7, 8, 9 with low Itgam/CD11b and Aif1/Iba1 expression
table(LG180_integrated@active.ident)
LG180_integrated <- subset(LG180_integrated, idents = c("6","7", "8", "9"), invert = TRUE)
table(LG180_integrated@active.ident)
saveRDS(LG180_integrated, file = 'elife_microglial_cells_only.rds')

