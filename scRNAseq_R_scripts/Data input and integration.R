library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)

setwd("/athena/ganlab/scratch/lif4001/Fan_microglia/data_analysis")

lib_LG180_A19.data <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Fan_microglia/LG180_A19/outs/filtered_feature_bc_matrix")
lib_LG180_A19.data_filter <- CreateSeuratObject(counts = lib_LG180_A19.data, min.cells = 10, min.features = 200, project = "LG180_A19")
lib_LG180_A19.data_filter[["Condition"]] = c('Ctrl')
lib_LG180_A19.data_filter[["Batch"]] = c('Batch-1')
lib_LG180_A19.data_filter[["sample_name"]] = c('Ctrl-1')

lib_LG180_A31.data <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Fan_microglia/cellranger_count_LG180_A31/outs/filtered_feature_bc_matrix")
lib_LG180_A31.data_filter <- CreateSeuratObject(counts = lib_LG180_A31.data, min.cells = 10, min.features = 200, project = "LG180_A31")
lib_LG180_A31.data_filter[["Condition"]] = c('Ctrl')
lib_LG180_A31.data_filter[["Batch"]] = c('Batch-2')
lib_LG180_A31.data_filter[["sample_name"]] = c('Ctrl-2')

lib_LG180_A32.data <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Fan_microglia/cellranger_count_LG180_A32/outs/filtered_feature_bc_matrix")
lib_LG180_A32.data_filter <- CreateSeuratObject(counts = lib_LG180_A32.data, min.cells = 10, min.features = 200, project = "LG180_A32")
lib_LG180_A32.data_filter[["Condition"]] = c('Ctrl')
lib_LG180_A32.data_filter[["Batch"]] = c('Batch-2')
lib_LG180_A32.data_filter[["sample_name"]] = c('Ctrl-3')

lib_LG180_A22A23.data <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Fan_microglia/LG180_A22A23/outs/filtered_feature_bc_matrix")
lib_LG180_A22A23.data_filter <- CreateSeuratObject(counts = lib_LG180_A22A23.data, min.cells = 10, min.features = 200, project = "LG180_A22A23")
lib_LG180_A22A23.data_filter[["Condition"]] = c('D0')
lib_LG180_A22A23.data_filter[["Batch"]] = c('Batch-1')
lib_LG180_A22A23.data_filter[["sample_name"]] = c('D0-1')

lib_LG180_A35A36.data <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Fan_microglia/LG180_A35A36/outs/filtered_feature_bc_matrix")
lib_LG180_A35A36.data_filter <- CreateSeuratObject(counts = lib_LG180_A35A36.data, min.cells = 10, min.features = 200, project = "LG180_A35A36")
lib_LG180_A35A36.data_filter[["Condition"]] = c('D0')
lib_LG180_A35A36.data_filter[["Batch"]] = c('Batch-2')
lib_LG180_A35A36.data_filter[["sample_name"]] = c('D0-2')

lib_LG180_A37A38.data <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Fan_microglia/LG180_A37A38/outs/filtered_feature_bc_matrix")
lib_LG180_A37A38.data_filter <- CreateSeuratObject(counts = lib_LG180_A37A38.data, min.cells = 10, min.features = 200, project = "LG180_A37A38")
lib_LG180_A37A38.data_filter[["Condition"]] = c('D0')
lib_LG180_A37A38.data_filter[["Batch"]] = c('Batch-2')
lib_LG180_A37A38.data_filter[["sample_name"]] = c('D0-3')

lib_LG180_A20A21.data <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Fan_microglia/LG180_A20A21/outs/filtered_feature_bc_matrix")
lib_LG180_A20A21.data_filter <- CreateSeuratObject(counts = lib_LG180_A20A21.data, min.cells = 10, min.features = 200, project = "LG180_A20A21")
lib_LG180_A20A21.data_filter[["Condition"]] = c('D2')
lib_LG180_A20A21.data_filter[["Batch"]] = c('Batch-1')
lib_LG180_A20A21.data_filter[["sample_name"]] = c('D2-1')

lib_LG180_A24A25.data <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Fan_microglia/LG180_A24A25/outs/filtered_feature_bc_matrix")
lib_LG180_A24A25.data_filter <- CreateSeuratObject(counts = lib_LG180_A24A25.data, min.cells = 10, min.features = 200, project = "LG180_A24A25")
lib_LG180_A24A25.data_filter[["Condition"]] = c('D2')
lib_LG180_A24A25.data_filter[["Batch"]] = c('Batch-1')
lib_LG180_A24A25.data_filter[["sample_name"]] = c('D2-2')

lib_LG180_A33A34.data <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Fan_microglia/cellranger_count_LG180_A33A34/outs/filtered_feature_bc_matrix")
lib_LG180_A33A34.data_filter <- CreateSeuratObject(counts = lib_LG180_A33A34.data, min.cells = 10, min.features = 200, project = "LG180_A33A34")
lib_LG180_A33A34.data_filter[["Condition"]] = c('D2')
lib_LG180_A33A34.data_filter[["Batch"]] = c('Batch-2')
lib_LG180_A33A34.data_filter[["sample_name"]] = c('D2-3')

rm(lib_LG180_A19.data, lib_LG180_A31.data, lib_LG180_A32.data, lib_LG180_A22A23.data, lib_LG180_A35A36.data, lib_LG180_A37A38.data, lib_LG180_A20A21.data, lib_LG180_A24A25.data, lib_LG180_A33A34.data)                                            

project_list <- c(lib_LG180_A19.data_filter, lib_LG180_A31.data_filter, lib_LG180_A32.data_filter, lib_LG180_A22A23.data_filter, lib_LG180_A35A36.data_filter, lib_LG180_A37A38.data_filter, lib_LG180_A20A21.data_filter, lib_LG180_A24A25.data_filter, lib_LG180_A33A34.data_filter)

for (i in 1:length(project_list)){
  project_list[[i]][['percent.mt']] <- PercentageFeatureSet(project_list[[i]], pattern = "^mt.")
  project_list[[i]] <- subset(x = project_list[[i]], subset = nCount_RNA < 20000 & percent.mt < 10)
  project_list[[i]] <- NormalizeData(project_list[[i]], verbose = FALSE)
  project_list[[i]] <- FindVariableFeatures(project_list[[i]], selection.method = 'vst', nfeatures = 3000)
}

LG180_A19 <- project_list[[1]]
LG180_A31 <- project_list[[2]]
LG180_A32 <- project_list[[3]]
LG180_A22A23 <- project_list[[4]]
LG180_A35A36 <- project_list[[5]]
LG180_A37A38 <- project_list[[6]]
LG180_A20A21 <- project_list[[7]]
LG180_A24A25 <- project_list[[8]]
LG180_A33A34 <- project_list[[9]]

rm(lib_LG180_A19.data_filter, lib_LG180_A31.data_filter, lib_LG180_A32.data_filter, lib_LG180_A22A23.data_filter, lib_LG180_A35A36.data_filter, lib_LG180_A37A38.data_filter, lib_LG180_A20A21.data_filter, lib_LG180_A24A25.data_filter, lib_LG180_A33A34.data_filter)
rm(project_list)

#1 by 1 by 1 by 1 Integration
project_list <- c(LG180_A19, LG180_A31, LG180_A32, LG180_A22A23, LG180_A35A36, LG180_A37A38, LG180_A20A21, LG180_A24A25, LG180_A33A34)
anchors <- FindIntegrationAnchors(object.list = project_list, dims = 1:30)
LG180_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
rm(LG180_A19, LG180_A31, LG180_A32, LG180_A22A23, LG180_A35A36, LG180_A37A38, LG180_A20A21, LG180_A24A25, LG180_A33A34, anchors)
saveRDS(LG180_integrated, file = 'integrated_LG180_1by1by1by1.rds')

#Scale integrated data
DefaultAssay(LG180_integrated) <- 'integrated'

all.genes <- rownames(LG180_integrated)
LG180_integrated <- ScaleData(LG180_integrated, features = all.genes)
LG180_integrated <- RunPCA(LG180_integrated, features = VariableFeatures(object = LG180_integrated))
ElbowPlot(LG180_integrated)

DimPlot(LG180_integrated, reduction = 'pca')
LG180_integrated <- FindNeighbors(LG180_integrated, dims = 1:10)
LG180_integrated <- FindClusters(LG180_integrated, resolution = 0.1)
LG180_integrated <- RunUMAP(LG180_integrated, dims = 1: 10, perplexity = 20)
DimPlot(LG180_integrated, reduction = 'umap', label = T)

DimPlot(LG180_integrated, reduction = "umap", split.by = "Condition", label = T)
DimPlot(LG180_integrated, reduction = "umap", split.by = "Batch", label = T)
DimPlot(LG180_integrated, reduction = "umap", split.by = "sample_name", label = T)

# rename clusters
n <- dim(table(LG180_integrated@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
LG180_integrated@active.ident <- plyr::mapvalues(x = LG180_integrated@active.ident, from = current.cluster.ids, to = new.cluster.ids)

LG180_integrated@active.ident <- factor(LG180_integrated@active.ident, levels=1:n)

saveRDS(LG180_integrated, file = 'elife_Ctrl_D0_D2.rds')



