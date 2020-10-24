# Supplementary Figure 6
# Comparing DEGs from Mac2+ vs homeostatic microglia with
# Hammond et al 2019 Supplementary Table 1 microglial clusters during development + aging
# Title: Single-Cell RNA Sequencing of Microglia throughout the Mouse Lifespan and in the Injured Brain Reveals Complex Cell-State Changes

library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# Only focus on genes that are overlapping with the dataset 
#######
# Load in DEGs from Mac2+ vs homeostatic microglia comparison
setwd("")
mac2<-read.csv("mac2_vs_homeostatic microglia.csv",header=T)
colnames(mac2)<-as.character(c("gene","p_val","avg_logFC","pct.1","pct.2","p_val_adj"))
mac2_up <- mac2$gene[mac2$avg_logFC>0&mac2$p_val_adj<0.05]
mac2_dn <- mac2$gene[mac2$avg_logFC<0&mac2$p_val_adj<0.05]

##################

# Use edited data set from Hammond et al including  fold change values
df <- read.csv("Hammond_direction.csv",header=T)

df$cluster[df$cluster==1]<-2 #"P4"
df$cluster[df$cluster==2]<-1 #"E/P"
df$cluster[df$cluster==3]<-1 #"E14.5"
df$cluster[df$cluster==4]<-2 #"P4"
df$cluster[df$cluster==5]<-1 #"E14.5"
df$cluster[df$cluster==6]<-1 #"E14.5"
df$cluster[df$cluster==7]<-3 #"Older"
df$cluster[df$cluster==8]<-4 #"Aged"
df$cluster[df$cluster==9]<-5 #"White"
df$cluster[df$cluster==10]<-4 #"Aged"
df$cluster[df$cluster==11]<-6 #"Mono"


# generate list of genes that overlap between data set and Mac2+ DEGs
dev_stage <- c("E14.5","P4", "P30","P540","White","Monocyte")
# calculate how many of the mac2 upregulated genes are overlapping with upregulated genes of dev genes
dev.up.overlap <- list()
for (i in 1:length(dev_stage)){
  dev_gene <- df$gene[df$cluster == i & df$foldchange > 1]
  dev_gene_up <- mac2_up[mac2_up %in% dev_gene]
  dev.up.overlap[[i]] <- dev_gene_up
}

dev.dn.overlap <- list() 
# calculate how many of the mac2 downregulated genes are overlapping with downregulated genes of dev genes
for (i in 1:length(dev_stage)){
  dev_gene <- df$gene[df$cluster == i & df$foldchange < 1]
  dev_gene_dn <- mac2_dn[mac2_dn %in% dev_gene]
  dev.dn.overlap[[i]] <- dev_gene_dn
}

# calculates % of genes in Mac2+ DEGs that belong to each dev_stage
percent_summary_mac2 <- data.frame()
percent_each_mac2 <- numeric()
for (i in 1:6) {
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

# Donut plot for upregulated DEGs
ggplot(percent_summary_mac2, aes(x=2, y=Percent_UP_MAC2_pos, fill=dev_stage)) +
  geom_bar(width = 1, stat = "identity", color="black", size=0.25) + coord_polar("y", start=0, direction = -1) +
  xlim(0.5, 2.5) + scale_fill_manual(values = my_color,
                                     name="MG Dev Genes", labels=my_labels)+ 
  
  geom_text(aes(x=0.5, y=0), label=paste("overlap",sum(unlist(lapply(dev.up.overlap,length))),sep="\n"), size=5) +
  ylab(NULL) + xlab(NULL) +
  theme(aspect.ratio = 1) +
  theme_minimal()+ 
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

ggsave("Mac2microglia_UP_devstages_Hammond_overlap.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop",
       scale = 0.8, width = 8, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)


#################
# Donut plot for downregulated DEGs
my_color <- rev(brewer.pal(nrow(percent_summary_mac2),"GnBu"))
my_labels <- paste(percent_summary_mac2$dev_stage, sprintf("%.2f",round(percent_summary_mac2$Percent_DN_MAC2_pos, 2)), sep=": ")
my_labels <- paste(my_labels, "%", sep="")

ggplot(percent_summary_mac2, aes(x=2, y=Percent_DN_MAC2_pos, fill=dev_stage)) +
  geom_bar(width = 1, stat = "identity", color="black", size=0.25) + coord_polar("y", start=0, direction = -1) +
  xlim(0.5, 2.5) + scale_fill_manual(values = my_color, breaks=percent_summary_mac2$dev_stage,
                                     name="MG Dev Genes", label=my_labels) +
  
  geom_text(aes(x=0.5, y=0), label=paste("overlap",sum(unlist(lapply(dev.dn.overlap,length))),sep="\n"), size=5) +
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

ggsave("Mac2microglia_DN_devstages_Hammond_overlap.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop",
       scale = 0.8, width = 8, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)
