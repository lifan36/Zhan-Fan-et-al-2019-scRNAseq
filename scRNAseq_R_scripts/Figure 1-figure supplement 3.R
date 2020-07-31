library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(patchwork)
setwd("")

LG180_integrated <- readRDS("elife_Ctrl_D0_D2.rds")

DefaultAssay(LG180_integrated) <- 'RNA'
LG180_integrated <- ScaleData(LG180_integrated)

#Figure 1-figure supplement 3a
VlnPlot(object = LG180_integrated, features = "Retnlg", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )
VlnPlot(object = LG180_integrated, features = "S100a9", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )

#Figure 1-figure supplement 3b
VlnPlot(object = LG180_integrated, features = "Ly6c1", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )
VlnPlot(object = LG180_integrated, features = "Ly6a", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )

#Figure 1-figure supplement 3c
VlnPlot(object = LG180_integrated, features = "Acta2", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") ) 
VlnPlot(object = LG180_integrated, features = "Vtn", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )
VlnPlot(object = LG180_integrated, features = "Tagln", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") ) 
VlnPlot(object = LG180_integrated, features = "Myl9", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )

#Figure 1-figure supplement 3d
VlnPlot(object = LG180_integrated, features = "Slc1a2", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )
VlnPlot(object = LG180_integrated, features = "Mt3", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )
VlnPlot(object = LG180_integrated, features = "Clu", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") ) 
VlnPlot(object = LG180_integrated, features = "Aldoc", pt.size = 0.1) & 
  theme( plot.title = element_text( face = "italic") )

