# Figure 4a

library(ggplot2)
library(Seurat)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
library(patchwork)
# Load in Seurat object
setwd("")
data<-readRDS("elife_microglial_cells_only.rds")
data@active.ident


# Find markers for Cluster 5:
marker <- FindMarkers(data,5, logfc.threshold = 0, test.use = "MAST")
write.csv(marker, "Cluster5_markers.csv")
# Volcano plot of marker genes =====

# Identify DEGs for cluster 5 markers
marker$colours <- c("NC")
marker$colours[marker$avg_logFC >= 0.5 & marker$p_val_adj <= 0.05] <- c("UP")
marker$colours[marker$avg_logFC <= -0.5 & marker$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_mature <- c("Cx3cr1", "Csf1r", "Mafb", "Tmem119")
genes_to_plot_mature <- marker[row.names(marker) %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- "Mature"

genes_select_immature <- c("Ifitm1", "Lgals3", "Il1b")
genes_to_plot_immature <- marker[row.names(marker)  %in% genes_select_immature, ]
genes_to_plot_immature$Cluster <- c("Immature")

genes_to_plot <- rbind(genes_to_plot_mature, genes_to_plot_immature)

# Set color palette
my_color <- c("#2B8CBE", "#D7301F", "skyblue","seashell3", "plum1")
my_color_1 <- c("#2B8CBE","Grey", "#D7301F")


# Plot volcano plot
dev.off()
ggplot() + 
  geom_point(data=marker, aes(x=avg_logFC, y=-log10(p_val_adj), colour=colours),
             shape=19, alpha=1, size=1) +
  scale_color_manual(values = my_color_1,
                     name="DEGs",
                     breaks=rev(names(table(marker$colours))),
                     labels=rev(names(table(marker$colours)))) +
  geom_point(data=genes_to_plot,
             aes(x=avg_logFC, y=-log10(p_val_adj)),
             shape=19, alpha=1, size=3) +
  geom_text_repel(data=genes_to_plot,
                  aes(x=avg_logFC, y=-log10(p_val_adj), label = row.names(genes_to_plot)), 
                  color="black", fontface = 'bold',size = 5, box.padding = 0.5,
                  point.padding = 0.5, segment.size=0.25, segment.colour="black") +
  ylab("-Log10[FDR]") + xlab("Log2FC") +
  ggtitle("Cluster 5 vs all other clusters")+
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=15),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=15),
        axis.title.y = element_text(colour = "black", size=15),
        plot.title = element_text(size = 15, face = "bold")) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks=seq(-4, 4, 2), limits=c(-4, 4))+
  NoLegend()

ggsave("Volcano_cluster_5.pdf", plot = last_plot(), device = "pdf", path = "",
       scale = 0.8, width = 7, height = 7, units = c("in"),
       dpi = 600, limitsize = FALSE)

# #Figure 4b
#change the color code
p1 <- FeaturePlot(LG180_integrated, features = c("Lgals3")) & 
  theme( plot.title = element_text( face = "italic") )
scale <- scale_color_gradientn(colours = brewer.pal(11,"RdYlBu")[11:1],limits=c(0,4))
p2 <- lapply(p1, function (x) x + scale)
CombinePlots(p2)

#Figure 4c
Cluster_5 <- subset(data, idents = "5")
VlnPlot(Cluster_5, features = "Lgals3", group.by = "Condition", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )

