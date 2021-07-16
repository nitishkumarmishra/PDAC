library(ggplot2)
library(ggthemes)
options(scipen = 999)  # turns of scientific notations like 1e+40

setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
# Read data

tmp <- stacked.plot
tmp$chr <- rownames(tmp)
tmp$total <- tmp$Hyper.per.Mb + tmp$Hypo.per.Mb
tmp1 <- tmp[,c("Hyper.per.Mb", "chr", "total")]
tmp2 <- tmp[,c("Hypo.per.Mb", "chr", "total")]
colnames(tmp2) <- gsub("Hypo.per.Mb", "Hyper.per.Mb", colnames(tmp2))
tmp2$Hyper.per.Mb <- as.numeric(paste0("-", tmp2$Hyper.per.Mb))

tmp1$Direction <- rep("Hypermethylated", 22)
tmp2$Direction <- rep("Hypomethylated", 22)
tmp.plot <- rbind(tmp1, tmp2)

colnames(tmp.plot) <- gsub("Hyper.per.Mb", "Mb", colnames(tmp.plot))
# X Axis Breaks and Labels 
brks <- seq(-10, 10, 5)
lbls = paste0(as.character(c(seq(10, 0, -5), seq(5, 10, 5))))
#paste0(as.character(c(seq(15, 0, -5), seq(5, 15, 5))), "m")
# Plot
ggplot(tmp.plot, aes(x = reorder(chr, -total), y = Mb, fill = Direction)) +   # Fill column
  xlab("Chromosome")+
  ylab("dm-CpG/Mb (Hyper/Hypo))")+
  geom_bar(stat = "identity", width = .6) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  coord_flip() +  # Flip axes
  labs(title="Differential hyper/hypomethylation frequency") +
  theme_tufte() +  # Tufte theme from ggfortify
  theme_bw()+
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks = element_blank()) +   # Centre plot title
  #scale_fill_brewer(palette = "Dark2")  # Color palette
  scale_fill_manual(values = c("red4",  "blue4"))+
  theme(legend.position = c(0.86, 0.9))+
  theme(text=element_text(size=13, color = "black", face="bold"),
    plot.title = element_text(color="black", size=14, face="bold"),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size=0.7),
    legend.background = element_rect(color = "black", fill = "white", size = 0.5, linetype = "solid"))
ggsave("Stacked_ggplo2.pdf",  dpi = 500, width = 7, height = 6, units = c("in"))
#dev.print(pdf, 'cg26597242 Boxplot.pdf ', width = 6, height = 6)
save.image("Stacked_ggplo2.RData")

Stacked <- read.csv("Stacked_CpG_data.txt", header = TRUE, sep = ",", row.names = 1)
Stacked$HyperVsHypo <-  (Stacked$Hyper.per.Mb/Stacked$Hypo.per.Mb)
Stacked$HypoVsHyper <- (Stacked$Hypo.per.Mb/Stacked$Hyper.per.Mb)
Stacked <- Stacked[order(Stacked$HyperVsHypo, decreasing = TRUE),]


write.csv(Stacked, "Stacked_CpG_data.txt")