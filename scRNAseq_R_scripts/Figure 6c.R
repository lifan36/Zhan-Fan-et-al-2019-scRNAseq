# Figure 6c

# Load Seurat object
setwd("")
data<-readRDS("elife_microglial_cells_only.rds")

# Saving the active ident with the 5 clusters used for this revision in the metadata
data$RevisedClusters <- Idents(data)

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
data$cluster_mac2 <- as.numeric(data$RevisedClusters)
data$cluster_mac2[row.names(data@meta.data) %in% cell_name] <- "mac2"
data$seurat_clusters_new <- as.numeric(data$RevisedClusters)
#### Lines 34 and 36 were changed from seurat_clusters to RevisedClusters after as.numeric()

# UMAP and visualizations of Mac2+ cells
Idents(data)<-"mac2"
DimPlot(data,cols=c("grey","orange"))
ggsave("umap mac2+.pdf", plot = last_plot(), device = "pdf", 
       scale = 0.8, width = 5, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE, path = "~/Desktop")

RidgePlot(data, features=c("Lgals3"), col=c("grey","orange"))
ggsave("ridge plot mac2+ microglia.pdf", plot = last_plot(), device = "pdf",
       scale = 0.8, width = 5.5, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE, path = "~/Desktop")
