library(ggplot2)
setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")

Bubble.data <- read.csv("Bubble_Plot.txt", header = TRUE, sep = "\t", check.names = FALSE)
Bubble.data$`Total dm-CpG` <- gsub("Hypermethylated", "Hyper", Bubble.data$`Total dm-CpG`)
Bubble.data$`Total dm-CpG` <- gsub("Hypomethylated", "Hypo", Bubble.data$`Total dm-CpG`)
Bubble.data$JitCoOr <- jitter(as.numeric(factor(Bubble.data$`Genomic Region`)))
Bubble.data$JitCoOrPow <- jitter(as.numeric(factor(Bubble.data$`Total dm-CpG`)))

ggplot(Bubble.data, aes(x = `Genomic Region`, y = `Total dm-CpG`)) +
  scale_fill_manual("Status",values=c("red4","blue4"))+
  geom_point(data=Bubble.data,aes(x=JitCoOr, y=JitCoOrPow,size = `dm-CpGs`, colour = `Total dm-CpG`), alpha=.5)+
  labs(title="Differential hyper/hypomethylation distribution") +
  xlab("Genomic Regions")+
  ylab("Total number of dm-CpGs")+
  theme(legend.position = c(0.5, 0.5))+
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust = 0.0), 
        text=element_text(size=10, color = "black", face="bold"),
        axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"),
        axis.ticks = element_blank()) + 
  geom_text(data=Bubble.data,aes(x=JitCoOr, y=JitCoOrPow,label=`Genomic Region`)) + 
  scale_size(range = c(1,60), guide = FALSE) +
  scale_y_discrete(breaks =1:3 , labels=c("Hyper","Hypo"," "), limits = c(1,2))+
  scale_x_discrete(breaks =1:12 , labels=c("TSS200", "TSS1500", "5'UTR","1st Exon","Gene Body","3'UTR", "Island", "S Shore", "N Shore", "S shelf", "N Self", "Promoter"))+ 
  theme_bw()

ggsave("Bubbleplot_CpG.pdf",  dpi = 600, width = 10, height = 6, units = c("in"))

save.image("Bubble_plot.RData")

### Bubble plot based on Hypermethylation/Hypomethylation ratio ###########
# 4705/7348 = 0.64
Bubble.data$Ratio <- ifelse(Bubble.data$`Total dm-CpG`=="Hypo", Bubble.data$`dm-CpGs`*1.25, Bubble.data$`dm-CpGs`)
ggplot(Bubble.data, aes(x = `Genomic Region`, y = Ratio)) +
  scale_fill_manual("Status",values=c("red4","blue4"))+
  geom_point(data=Bubble.data,aes(x=JitCoOr, y=JitCoOrPow,size = Ratio, colour = `Total dm-CpG`), alpha=.5)+
  labs(title="Differential hyper/hypomethylation distribution") +
  xlab("Genomic Regions")+
  ylab("Total number of dm-CpGs")+
  theme(legend.position = c(0.5, 0.5))+
  theme(plot.title = element_text(color="black", size=14, face="bold", hjust = 0.0), 
        text=element_text(size=10, color = "black", face="bold"),
        axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"),
        axis.ticks = element_blank()) + 
  geom_text(data=Bubble.data,aes(x=JitCoOr, y=JitCoOrPow,label=`Genomic Region`)) + 
  scale_size(range = c(1,60), guide = FALSE) +
  scale_y_discrete(breaks =1:3 , labels=c("Hyper","Hypo"," "), limits = c(1,2))+
  scale_x_discrete(breaks =1:12 , labels=c("TSS200", "TSS1500", "5'UTR","1st Exon","Gene Body","3'UTR", "Island", "S Shore", "N Shore", "S shelf", "N Self", "Promoter"))+ 
  theme_bw()

ggsave("Bubbleplot_CpG_Ratio_1.pdf",  dpi = 600, width = 10, height = 6, units = c("in"))
save.image("Bubble_plot.RData")


#####################################################################
###################### Bubble plot revised ##########################
ggplot(tmp, aes(x = GenomicRegion, y = Ratio, show.legend=TRUE)) +
  scale_colour_manual(values=c("red4","blue4"))+
  geom_point(data=tmp,aes(x=GenomicRegion, y=JitCoOrPow, size = Ratio, colour = Status), alpha=.8)+
  labs(title="Differential hyper/hypomethylation distribution") +
  xlab("GenomicRegions")+
  ylab("Total number of dm-CpGs")+
  theme(legend.position = c(0.5, 0.5))+
  geom_text(data=tmp,aes(x=GenomicRegion, y=JitCoOrPow,label=`GenomicRegion`))+
  scale_size(range = c(1,60), guide = FALSE, trans="log10") +
  theme(legend.position="bottomright")+
  scale_y_discrete(breaks =1:3 , labels=c("Hyper","Hypo"," "), limits = c(1,2))+
  scale_x_discrete(breaks =1:12 ) +
  theme_bw()

ggsave("Bubbleplot_CpG_Ratio_1.pdf",  dpi = 600, width = 15, height = 6, units = c("in"))

