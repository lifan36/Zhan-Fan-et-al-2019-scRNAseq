# Figure 6de

library(dplyr)
library(Seurat)
library(MAST)
library(RColorBrewer)

### the microglial progenitor gene signature file "MG_dev.csv" was obtained from
### the 2016 science paper (supplementary table) by Matcovitch-Natan and Winter et al, 
### title: "Microglia development follows a stepwise program to regulate brain homeostasis" 

# Load DEG list for Mac2+ vs homeostatic microglia from Figure 6c
setwd("")
mac2<-read.csv("mac2_vs_homeostatic microglia.csv",header=T)
colnames(mac2)<-as.character(c("gene","p_val","avg_logFC","pct.1","pct.2","p_val_adj"))
mac2_up <- mac2$gene[mac2$avg_logFC>0&mac2$p_val_adj<0.05]
mac2_dn <- mac2$gene[mac2$avg_logFC<0&mac2$p_val_adj<0.05]

# Load developmental genes identified in Matcovitch-Natan and Winter et al
setwd("")
df <- read.csv("Dev_Genes.csv",header=T)
YS_genes <- df$Gene[df$Cluster==1]
E1_genes <- df$Gene[df$Cluster==2]
E2_genes <- df$Gene[df$Cluster==3]
P1_genes <- df$Gene[df$Cluster==4]
P2_genes <- df$Gene[df$Cluster==5]
A1_genes <- df$Gene[df$Cluster==6]
A2_genes <- df$Gene[df$Cluster==7]

##################
# find overlapping genes in Mac2+ upregulated and downregulated that belong to each dev_stage
dev_stage <- c("YS", "E1", "E2", "P1", "P2", "A1", "A2")
dev.up.overlap <- list()
for (i in 1:length(dev_stage)){
  dev_gene <- df$Gene[df$Cluster == i]
  dev_gene_up <- mac2_up[mac2_up %in% dev_gene]
  dev.up.overlap[[i]] <- dev_gene_up
}

dev.dn.overlap <- list()
for (i in 1:length(dev_stage)){
  dev_gene <- df$Gene[df$Cluster == i]
  dev_gene_dn <- mac2_dn[mac2_dn %in% dev_gene]
  dev.dn.overlap[[i]] <- dev_gene_dn
}

# calculates % of genes in mac2+ DEGs that belong to each dev_stage
percent_summary_mac2 <- data.frame()
percent_each_mac2 <- numeric()
for (i in 1:7) {
  dev_gene_up <- dev.up.overlap[[i]]
  percent_up <- length(dev_gene_up)/sum(unlist(lapply(dev.up.overlap,length)))*100
  dev_gene_dn <- dev.dn.overlap[[i]]
  percent_dn <- length(dev_gene_dn)/sum(unlist(lapply(dev.dn.overlap,length)))*100
  percent_each_mac2 <- c(percent_up, percent_dn)
  percent_summary_mac2 <- rbind(percent_summary_mac2,percent_each_mac2)
}


colnames(percent_summary_mac2) <- c("Percent_UP_MAC2_pos", "Percent_DN_MAC2_pos")
percent_summary_mac2$dev_stage <- c(dev_stage)
percent_summary_mac2$dev_stage <- factor(percent_summary_mac2$dev_stage, levels = percent_summary_mac2$dev_stage)

# Set color palette
my_color <- rev(brewer.pal(nrow(percent_summary_mac2),"OrRd"))[1:nrow(percent_summary_mac2)]
my_labels <- paste(percent_summary_mac2$dev_stage, sprintf("%.2f",round(percent_summary_mac2$Percent_UP_MAC2_pos, 2)), sep=": ")
my_labels <- paste(my_labels, "%", sep="")

# Plot donut chart
ggplot(percent_summary_mac2, aes(x=2, y=Percent_UP_MAC2_pos, fill=dev_stage)) +
  geom_bar(width = 1, stat = "identity", color="black", size=0.25) + coord_polar("y", start=0, direction = -1) +
  xlim(0.5, 2.5) + scale_fill_manual(values = my_color,
                                     name="MG Dev Genes", labels=my_labels) +
  geom_text(aes(x=0.5, y=0), label=paste("overlap",sum(unlist(lapply(dev.up.overlap,length))),sep="\n"), size=8) +
  ylab(NULL) + xlab(NULL) +
  theme(aspect.ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=14),
        legend.position='right',
        legend.direction = "vertical",
        axis.text.x = element_text(colour = "black", size=0),
        panel.grid.major  = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank()) 

ggsave("Mac2microglia_UP_devstages_overlap_notitle.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop",
       scale = 0.9, width = 8, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)


#################
percent_summary_mac2_reorder <- percent_summary_mac2[c(7:1),]
percent_summary_mac2_reorder$dev_stage <- factor(percent_summary_mac2_reorder$dev_stage, levels=percent_summary_mac2_reorder$dev_stage)
my_color <- rev(brewer.pal(nrow(percent_summary_mac2),"GnBu"))

my_labels <- paste(percent_summary_mac2_reorder$dev_stage, sprintf("%.2f",round(percent_summary_mac2_reorder$Percent_DN_MAC2_pos, 2)), sep=": ")
my_labels <- paste(my_labels, "%", sep="")

ggplot(percent_summary_mac2_reorder, aes(x=2, y=Percent_DN_MAC2_pos, fill=dev_stage)) +
  geom_bar(width = 1, stat = "identity", color="black", size=0.25) + coord_polar("y", start=0, direction = -1) +
  xlim(0.5, 2.5) + scale_fill_manual(values = my_color, breaks=percent_summary_mac2_reorder$dev_stage,
                                     name="MG Dev Genes", label=my_labels) +
  
  geom_text(aes(x=0.5, y=0), label=paste("overlap",sum(unlist(lapply(dev.dn.overlap,length))),sep="\n"), size=8) +
  ylab(NULL) + xlab(NULL) +
  theme(aspect.ratio = 1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        legend.title = element_text(colour="black", size=14, face="bold"),
        legend.text = element_text(colour="black", size=14),
        legend.position='right',
        legend.direction = "vertical",
        axis.text.x = element_text(colour = "black", size=0),
        panel.grid.major  = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())

ggsave("Mac2microglia_DN_devstages_overlap_notitle.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop",
       scale = 0.9, width = 8, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)
