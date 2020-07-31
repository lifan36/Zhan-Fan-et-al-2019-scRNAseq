library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)
library(patchwork)

setwd("")

LG180_integrated <- readRDS("elife_microglial_cells_only.rds")

#Figure 2a
LG180_markers <- FindAllMarkers(LG180_integrated, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(LG180_markers, "LG180_markers.csv")
top5 <- LG180_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pdf("HeatMapTop5_all_0.1.pdf", width=16, height=8)
DoHeatmap(LG180_integrated, features = top5$gene) + NoLegend()
dev.off()

#Figure 2b, c, d
Cluster_3 <- subset(LG180_integrated, idents = "3") 
Cluster_4 <- subset(LG180_integrated, idents = "4")
Cluster_5 <- subset(LG180_integrated, idents = "5")

###Cluster-3:Rps8 and Rps26
VlnPlot(Cluster_3, features = c("Rps8","Rps26"), group.by = "Condition", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )
###Cluster-4:Mki67 and Top2a
VlnPlot(Cluster_4, features = c("Mki67","Top2a"), group.by = "Condition", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )
###Cluster-5:Cd74 and H2-Aa
VlnPlot(Cluster_5, features = c("Cd74","H2-Aa"), group.by = "Condition", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )

#Figure 2e
RidgePlot(LG180_integrated, features=c("Cd74"))

