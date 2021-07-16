library(ggplot2)
setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")


Bubble.data <- read.csv("Bubble_Plot - Copy.txt", header = TRUE, sep = "\t", check.names = FALSE)
Bubble.data$`Total dm-CpG` <- gsub("Hypermethylated", "Hyper", Bubble.data$`Total dm-CpG`)
Bubble.data$`Total dm-CpG` <- gsub("Hypomethylated", "Hypo", Bubble.data$`Total dm-CpG`)
Bubble.data$JitCoOr <- jitter(as.numeric(factor(Bubble.data$GenomicRegion)))
Bubble.data$JitCoOrPow <- jitter(as.numeric(factor(Bubble.data$`Total dm-CpG`)))
 
tmp <- Bubble.data
tmp$Ratio <- ifelse(Bubble.data$`Total dm-CpG`=="Hypo", Bubble.data$`dm-CpGs`*1.25, Bubble.data$`dm-CpGs`)

## Fix order
tmp$GenomicRegion <- factor(tmp$GenomicRegion, levels = c("TSS200","TSS1500","Promoter","5'UTR","1st Exon","Gene Body",
                                                                            "3'UTR","Island","S Shore","N Shore","S shelf","N Self"))



ggplot(tmp, aes(x = GenomicRegion, y = Ratio, show.legend=FALSE)) +
  scale_colour_manual(values=c("red","blue"))+
  geom_point(data=tmp,aes(x=GenomicRegion, y=JitCoOrPow, size = Ratio, colour = `Total dm-CpG`), alpha=.5)+
  labs(title="Differential hyper/hypomethylation distribution") +
  xlab("GenomicRegions")+
  ylab("Total number of dm-CpGs")+
  theme(legend.position = c(0.5, 0.5))+
  geom_text(data=tmp,aes(x=GenomicRegion, y=JitCoOrPow,label=GenomicRegion))+
  scale_size(range = c(1,60), guide = FALSE, trans="log10") + 
  #theme(legend.position="bottomright")+ 

  scale_y_discrete(breaks =1:3 , labels=c("Hyper","Hypo"," "), limits = c(1,2))+
  scale_x_discrete(breaks =1:12, labels=c("TSS200", "TSS1500", "5'UTR","1st Exon","Gene Body","3'UTR", "Island", "S Shore", "N Shore", "S shelf", "N Self", "Promoter") ) +
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal") 



ggsave("foo300_1.pdf",width = 320, height = 160, dpi = 800, units = "mm")
ggsave(filename = "foo300.jpg", width = 400, height = 160, dpi = 800, units = "mm", device='jpg')
