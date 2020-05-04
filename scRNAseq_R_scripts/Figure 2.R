library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("~/Desktop/data_analysis/elife/D2/all_0.1_local")

LG180_integrated <- readRDS("elife_microglial_cells_only.rds")

#Figure 2a
LG180_markers <- FindAllMarkers(LG180_integrated, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(LG180_markers, "LG180_markers.csv")
top5 <- LG180_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
pdf("HeatMapTop5_all_0.1.pdf", width=10, height=8)
DoHeatmap(LG180_integrated, features = top5$gene) + NoLegend()
dev.off()

#Figure 2b
Cluster_3 <- subset(LG180_integrated, idents = "3")
Cluster_4 <- subset(LG180_integrated, idents = "4")
Cluster_5 <- subset(LG180_integrated, idents = "5")
Cluster_6 <- subset(LG180_integrated, idents = "6")

###Cluster-3:Rps8 and Rps26
VlnPlot(Cluster_3, features = c("Rps8","Rps26"), group.by = "Condition", pt.size = 0.1, combine = FALSE)
###Cluster-4:Mki67 and Top2a
VlnPlot(Cluster_3, features = c("Mki67","Top2a"), group.by = "Condition", pt.size = 0.1, combine = FALSE)
###Cluster-5:Cd74 and H2-Aa
VlnPlot(Cluster_3, features = c("Cd74","H2-Aa"), group.by = "Condition", pt.size = 0.1, combine = FALSE)
###Cluster-6:Ngp and Camp
VlnPlot(Cluster_3, features = c("Ngp","Camp"), group.by = "Condition", pt.size = 0.1, combine = FALSE)
