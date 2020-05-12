# Figure 6c
# Compare Mac2+ cell transcriptome with homeostatic microglia cluster (MG1 and MG2) from control condition

# Load Seurat object
setwd("")
data<-readRDS("elife_Ctrl_D0_D2.rds")

# Remove clusters 7-9
Idents(data) <- "seurat_clusters_new"
data <- subset(data,idents = c(1:6))

# Label clusters 1 and 2 as homeostatic clusters
# Extract contol conditions from this homeostatic cluster
# Compare with Mac2+ cells

data$cluster_mac2 <- data$seurat_clusters_new
data$cluster_mac2[data$mac2 == "mac2"] <- "mac2"
data$cluster_mac2[data$Condition == "Ctrl" & data$seurat_clusters_new == 1] <- "Ctrl_homeo"
data$cluster_mac2[data$Condition == "Ctrl" & data$seurat_clusters_new == 2] <- "Ctrl_homeo"

Idents(data)<-"cluster_mac2"
marker <- FindMarkers(data,ident.1="mac2",ident.2="Ctrl_homeo", logfc.threshold = 0,
                      test.use = "MAST")

library(ggplot2)
library(ggrepel)

# Volcano plot of DEGs =====

# Identify DEGs for Mac2+ vs homeostatic control microglia
marker$colours <- c("NC")
marker$colours[marker$avg_logFC >= 0.5 & marker$p_val_adj <= 0.05] <- c("UP")
marker$colours[marker$avg_logFC <= -0.5 & marker$p_val_adj <= 0.05] <- c("DN")

# Select genes of interest to highlight
genes_select_mature <- c("Cx3cr1", "Csf1r", "Mafb", "Tmem119")
genes_to_plot_mature <- marker[row.names(marker) %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- "Mature"

genes_select_immature <- c("Lyz2", "Mmp8", "Mmp9", "Pf4", "Ifit3", "Cdk1", "Lgals3")
genes_to_plot_immature <- marker[row.names(marker)  %in% genes_select_immature, ]
genes_to_plot_immature$Cluster <- c("Immature")

genes_to_plot <- rbind(genes_to_plot_mature, genes_to_plot_immature)

# Set color palette
my_color <- c("#2B8CBE", "#D7301F", "skyblue","seashell3", "plum1")
my_color_1 <- c("#D7301F","Grey", "#2B8CBE")

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
  ggtitle("Mac2+ cells vs Homeostatic MG") +
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=20),
        axis.text.y = element_text(colour = "black", size=20),
        axis.title.x = element_text(colour = "black", size=20),
        axis.title.y = element_text(colour = "black", size=20),
        plot.title= element_text(size = 20, face = "bold")) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks=seq(-4, 4, 2), limits=c(-4, 4))

ggsave("Volcano_Mac2_pos_genes.pdf", plot = last_plot(), device = "pdf", path = "",
       scale = 0.8, width = 7, height = 7, units = c("in"),
       dpi = 600, limitsize = FALSE)
