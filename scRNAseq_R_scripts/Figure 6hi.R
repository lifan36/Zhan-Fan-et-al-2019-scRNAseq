# Figure 6hi 
# Bargraphs for GSEA data and scatterplot for IPA data
# Input Mac2+ vs homeostatic microglia comparison DEGs

library(stringr)
library(dplyr)
library(ggplot2)

setwd("")
# input for GSEA analysis
mac2<-read.csv("mac2_vs_homeostatic microglia.csv",header=T)
mac2_dn<-mac2[mac2$avg_logFC<0&mac2$p_val_adj<0.05,]
mac2_dn<-mac2_dn[1:1000,]
mac2_up<-mac2[mac2$avg_logFC>0&mac2$p_val_adj<0.05,]
mac2_up<-mac2_up[1:1000,]
mac2_DE<-rbind(mac2_up,mac2_dn)
write.csv(mac2_DE,"mac2_DEGs.csv")

# set the directory where your csv file is located
setwd("")

df <- read.csv("GSEA_overlap.csv",header=T)
colnames(df) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR")

df$FDR <- as.numeric(df$FDR)

BP_1 <- dplyr::mutate(df, minus_logp=0-log10(FDR))
BP_1 <- dplyr::mutate(BP_1, enrich=paste(BP_1$genes_in_data, "/", BP_1$genes_in_gsea, sep=" ") )
BP_2 <- BP_1[order(-BP_1$minus_logp),]
BP_2$Name <-  gsub("HALLMARK_", "", BP_2$name)
BP_2$Name <- gsub("*_", " ", BP_2$Name)
BP_2$Name <-  factor(BP_2$Name, levels=rev(BP_2$Name))

df_summary <- data.frame(name=BP_2$Name, cluster=files[i])
data_conslidate <- rbind(data_conslidate, df_summary)

ggplot(data=BP_2, aes(x=Name  , y=minus_logp)) +
  theme_bw()+
  theme(panel.grid.major.x  = element_line(size =0.2, colour="grey88"),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(colour = "black", size=12),
        axis.text.y = element_text(colour = "black", size=12),
        axis.title.x = element_text(colour = "black", size=14),
        axis.title.y = element_text(colour = "black", size=14)) +
  ylab("-Log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity",  width=0.7, fill="grey", alpha=0.8) +
  # maroon1 dodgerblue
  geom_text(aes(label=enrich), vjust=0.4, 
            hjust=0, size=4, color="black", 
            stat="Identity", y=0.05*max(BP_2$minus_logp)) +
  coord_flip() + 
  theme(aspect.ratio = 1.5)

ggsave(paste("GSEA_all", ".pdf", sep=""), plot = last_plot(), device = "pdf", path = "",
       scale = 0.5, width = 12, height = 6, units = c("in"),
       dpi = 600, limitsize = FALSE)

#IPA results
ipa<-read.csv("IPA_upstream.csv",header=T)
ipa <- subset(ipa, p.value.of.overlap < 0.05)
ipa$colors <- "grey"
ipa$colors[ipa$Activation.z.score>0]<-"red"
ipa$colors[ipa$Activation.z.score<0]<-"blue"

ggplot(data=ipa,aes(x=-log(p.value.of.overlap),y=Activation.z.score,label=Upstream.Regulator,color=colors))+
  geom_point(show.legend=F,size=3)+theme_classic()+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=-log(0.05),linetype="dotted")+
  geom_text_repel(color="black")+
  scale_color_manual(values=c("red"="red","blue"="blue","grey"="grey"))+
  ylab("Predicted Z score")+ xlab("-Log10(p-value)")

ggsave(paste("IPA z scores", ".pdf", sep=""), plot = last_plot(), device = "pdf", path = "",
       scale = 1, width = 4.5, height = 4, units = c("in"),
       dpi = 600, limitsize = FALSE)


