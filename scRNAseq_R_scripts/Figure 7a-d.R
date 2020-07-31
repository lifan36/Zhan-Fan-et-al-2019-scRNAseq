# Generating Figure 7a-d plots
# Compare Ctrl homeostatic microglia vs Ctrl Mac2+ cells

# Load Seurat object
setwd("")
data<-readRDS("elife_Ctrl_D0_D2.rds")

# Remove clusters 7-9
#Idents(data) <- "seurat_clusters_new"
#data <- subset(data,idents = c(1:6))

# Extract raw lgals3 data
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

# Label clusters 1 and 2 as homeostatic clusters
# Extract contol conditions from this homeostatic cluster
# Compare with Mac2+ cells

data$cluster_mac2 <- as.numeric(data$seurat_clusters)
data$cluster_mac2[data$Condition == "Ctrl" & data$mac2 == "mac2"] <- "Ctrl_mac2"
data$cluster_mac2[data$Condition == "Ctrl" & data$seurat_clusters_new == 1] <- "Ctrl_homeo"
data$cluster_mac2[data$Condition == "Ctrl" & data$seurat_clusters_new == 2] <- "Ctrl_homeo"
table(data$cluster_mac2)

# Comparing Mac2+ cells from control samples vs Homeostatic microglia from control samples
Idents(data)<-"cluster_mac2"
marker <- FindMarkers(data,ident.1="Ctrl_mac2",ident.2="Ctrl_homeo", logfc.threshold = 0,
                      test.use = "MAST")

write.csv(marker,"mac2_ctrl_vs_homeostatic.csv")

# Bar plot of number of Mac2+ cells per condition
library(ggplot2)
mac2.number <- as.data.frame(table(data$mac2, data$Condition))
colnames(mac2.number) <- c("Mac2","Condition","Cell.number")
total.number <- as.data.frame(table(data$Condition))
colnames(total.number) <- c("Condition","Total.no")
mac2.only <- subset(mac2.number, Mac2 == "mac2")
mac2.only <- merge(mac2.only,total.number)
mac2.only$percent <- paste(round(mac2.only$Cell.number / mac2.only$Total.no * 100,digits=1), "%", sep="")

ggplot(mac2.only,aes(x=Condition,y=Cell.number,fill=Condition,label=percent))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("Ctrl" = "grey", "D0" = "maroon1", "D2" = "dodgerblue"))+
  theme_classic()+
  geom_text(vjust = -1)+
  ylim(0,1500)+
  ylab("Number of Mac2+ cells")+
  NoLegend()

ggsave("Mac2+ cells by condition.pdf",plot = last_plot(),
       path = "~/Desktop",
       width = 3, height = 3.5, units = "in")

# Volcano plot of DEGs from Mac2+ vs homeostatic MG from control samples =====
marker$colours <- c("NC")
marker$colours[marker$avg_logFC >= 0.5 & marker$p_val_adj <= 0.05] <- c("UP")
marker$colours[marker$avg_logFC <= -0.5 & marker$p_val_adj <= 0.05] <- c("DN")

# selected genes of interest to highlight
genes_select_mature <- c("Cx3cr1", "Csf1r", "Mafb", "Tmem119")
genes_to_plot_mature <- marker[row.names(marker) %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- "Mature"

genes_select_immature <- c("Lyz2", "Mmp8", "Mmp9", "Pf4",  "Cdk1", "Lgals3")
genes_to_plot_immature <- marker[row.names(marker)  %in% genes_select_immature, ]
genes_to_plot_immature$Cluster <- c("Immature")

genes_to_plot <- rbind(genes_to_plot_mature, genes_to_plot_immature)

# set color palette
my_color <- c("#2B8CBE", "#D7301F", "skyblue","seashell3", "plum1")
my_color_1 <- c("#D7301F","Grey", "#2B8CBE")

# plot volcano plot
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
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=20),
        axis.text.y = element_text(colour = "black", size=20),
        axis.title.x = element_text(colour = "black", size=20),
        axis.title.y = element_text(colour = "black", size=20),
        plot.title = element_text(size = 20, face = "bold")) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks=seq(-4, 4, 2), limits=c(-4, 4))+
  NoLegend()

ggsave("Volcano_CtrlMac2pos_v_homeostatic.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop",
       scale = 0.8, width = 7, height = 7, units = c("in"),
       dpi = 600, limitsize = FALSE)

# compare DEG numbers between Mac2+ control vs MG1 and Mac2+ all vs MG1
all.marker<-read.csv("mac2_vs_homeostatic microglia.csv",header=T)
all.marker.sig <-subset(all.marker,p_val_adj<0.05)

library(VennDiagram)

all.marker.up <- all.marker.sig$X[all.marker.sig$avg_logFC>0]
all.marker.dn <- all.marker.sig$X[all.marker.sig$avg_logFC<0]

marker.sig <- subset(marker, p_val_adj < 0.05)
marker.up<-row.names(marker.sig)[marker.sig$avg_logFC>0]
marker.dn<-row.names(marker.sig)[marker.sig$avg_logFC<0]

overlap.up <- all.marker.up[all.marker.up %in% marker.up]
overlap.dn <- all.marker.dn[all.marker.dn %in% marker.dn]


## cat.pos - position of category titles, represented by degree from the
## middle of the circle

## cat.dist - distance of the category titles from the edge of the circle


dev.off()
# draw venn with upregulated genes
g <- draw.pairwise.venn(area1 = length(all.marker.up), area2 = length(marker.up), cross.area = length(overlap.up), 
                        category = c("Mac2+ all", "Mac2+ ctrl"),fontfamily="sans", lwd=c(5,5),
                        fill = c("maroon1", "maroon1"), col = c("white", "white"), 
                        alpha = c(0.5, 0.5), cat.pos = c(0, 0), cat.dist = c(0.02, 0.02))

ggsave("venn_up_mac2_mg.pdf", plot = g, device = "pdf", path = "~/Desktop",
       scale = 1, width = 10, height = 10, units = c("in"),
       dpi = 300, limitsize = FALSE)

#############

adjustcolor( "dodgerblue", alpha.f = 1)
dev.off()
g <- draw.pairwise.venn(area1 = length(all.marker.dn), area2 = length(marker.dn), cross.area = length(overlap.dn), 
                        category = c("Mac2+ all", "Mac2+ ctrl"),fontfamily="sans", lwd=c(5,5),
                        fill = c("dodgerblue", "dodgerblue"), col = c("white", "white"), 
                        alpha = c(0.5, 0.5), cat.pos = c(0, 0), cat.dist = c(0.02, 0.02))


ggsave("venn_dn_mac2_mg.pdf", plot = g, device = "pdf", path = "~/Desktop",
       scale = 1, width = 10, height = 10, units = c("in"),
       dpi = 300, limitsize = FALSE)
