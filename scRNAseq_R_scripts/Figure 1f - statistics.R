# Figure 1f
# Calculating differential abundance for cluster distributions using EdgeR 

library(edgeR)
library(Seurat)

# Load in Seurat Object
setwd("~/Documents/Lab/Lihong revision/D2_included")
data <- readRDS("elife_Ctrl_D0_D2.rds")

# Subsetting Seurat data for each pair of comparisons
all.data <- list()

Idents(data) <- "Condition"
D0 <- subset(data,Condition == "D2", invert = T)
D2 <- subset(data,Condition == "D1", invert = T)
D0_D2 <- subset(data, Condition == "Ctrl", invert = T)

all.data[[1]] <- D0
all.data[[2]] <- D2
all.data[[3]] <- D0_D2

# Use GLM to analyze statistics of condition for cell abundance changes

comparisons <- c("D0vsCtrl","D2vsCtrl","D0vsD2")

for (i in 1:3){
  a <- all.data[[i]]
  
  # Table of abundances by Condition
  abundances <- table(a$seurat_clusters_new,a$sample_name)
  abundances <- unclass(abundances)
  
  # Extract meta data from Seurat object
  metadata <- a@meta.data
  extra.info <- metadata[match(colnames(abundances), a$sample_name),]
  
  # Create EdgeR object
  y.ab <- DGEList(abundances, samples=extra.info)
  y.ab
  
  # Design matrix
  design <- model.matrix(~factor(Batch) + factor(Condition) , y.ab$samples)
  
  # Estimate dispersion
  y.ab <- estimateDisp(y.ab, design, trend="none")
  summary(y.ab$common.dispersion)
  plotBCV(y.ab, cex=1)
  
  # QL dispersion
  fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
  summary(fit.ab$var.prior)
  summary(fit.ab$df.prior)
  plotQLDisp(fit.ab, cex=1)
  
  # Calculate GLM for Condition of sample
  res <- glmQLFTest(fit.ab, coef=ncol(design))
  summary(decideTests(res))
  results <- topTags(res)
  
  # Write data
  write.csv(results,paste("DA_",comparisons[i],".csv",sep=""))
}





