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
Idents(data) <- "Condition"


# Find markers for D0 vs Ctrl:
marker <- FindMarkers(data, ident.1 = "D0", ident.2 = "Ctrl", logfc.threshold = 0, test.use = "MAST")
write.csv(marker, "D0vsCtrl_markers.csv")
# Volcano plot of marker genes =====

# Identify DEGs for cluster 5 markers
marker$colours <- c("NC")
marker$colours[marker$avg_logFC >= 0.5 & marker$p_val_adj <= 0.05] <- c("UP")
marker$colours[marker$avg_logFC <= -0.5 & marker$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_mature <- c("Cx3cr1", "Csf1r", "Mafb", "Tmem119","P2ry12","ApoE")
genes_to_plot_mature <- marker[row.names(marker) %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- "Mature"

genes_select_immature <- c("Lgals3", "Il1b","Cd52","Cd74","Lyz2")
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
  ggtitle("D0 vs Ctrl")+
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
  scale_x_continuous(breaks=seq(-2, 2, 1), limits=c(-2, 2))+
  NoLegend()

ggsave("Volcano_D0vsCtrl.pdf", plot = last_plot(), device = "pdf", path = "",
       scale = 0.8, width = 7, height = 7, units = c("in"),
       dpi = 600, limitsize = FALSE)

# #Figure 4b
#change the color code
FeaturePlot(LG180_integrated, features = c("Lgals3")) & 
  theme( plot.title = element_text( face = "italic") )

#Figure 4c
Cluster_4 <- subset(data, idents = "4")
VlnPlot(Cluster_4, features = "Lgals3", group.by = "Condition", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )

