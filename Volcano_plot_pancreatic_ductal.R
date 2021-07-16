library(ggplot2)
library(ggrepel)
#genes <- as.data.frame(resSig)
#genes$Gene <- sapply(strsplit(rownames(genes),"\\|"),'[[',1)


TCGA.Harmonized <- read.table("TCGA_GDC_Harmonided_Uniq_Status_GTF.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)



diffExp.DESeq2 <- read.csv("diffExp_DESeq_ALL_Genes_10_19_2018.txt",header = TRUE)
rownames(diffExp.DESeq2) <- diffExp.DESeq2$X


diffExp.DESeq2$Gene_Symbol <- TCGA.Harmonized$Gene_Symbol[(match(rownames(diffExp.DESeq2),TCGA.Harmonized$ENSB_ID))]
diffExp.DESeq2$Gene_status <- TCGA.Harmonized$Gene_status[(match(rownames(diffExp.DESeq2),TCGA.Harmonized$ENSB_ID))]


all_diff <- na.omit(diffExp.DESeq2) # remove genes with no annotation


genes <- all_diff  # diff_Exp_genes 

Gene_Type <- "protein_coding"

genes <- all_diff[all_diff$Gene_status== Gene_Type,]

#for(Gene_Type in c("protein_coding", "lincRNA")){

  
  #genes <- all_diff[all_diff$Gene_status== Gene_Type,]


x.cut=2;y.cut=0.05

Significance <- ifelse(genes$log2FoldChange >= x.cut & genes$padj < y.cut, "Upregulated", ifelse(genes$log2FoldChange <= -x.cut & genes$padj < y.cut, "Downregulated", "Not significant"))

ggplot(genes, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), size = 2) +
  scale_color_manual(values = c( "blue4","gray41" ,"red4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  geom_text_repel(data = subset(genes, (padj < 0.03 & abs(log2FoldChange) >=2)),
                  aes(label = Gene_Symbol)) +
  
  theme(legend.position="right")+
  ggtitle("Volcano plot")+
  xlab("log2(FoldChange)") + ylab("-log10(padj") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 


ggsave(filename = paste0(Gene_Type,"VolcanoPlot_10_20_2019.pdf"), width = 12, height = 8, dpi = 800)


#}




### Run this for lincRNA - as there is no downregulated lincRNA - change in color

genes <- all_diff  # diff_Exp_genes 

Gene_Type <- "lincRNA"

genes <- all_diff[all_diff$Gene_status== Gene_Type,]

x.cut=2;y.cut=0.05

Significance <- ifelse(genes$log2FoldChange >= x.cut & genes$padj < y.cut, "Upregulated", ifelse(genes$log2FoldChange <= -x.cut & genes$padj < y.cut, "Downregulated", "Not significant"))

ggplot(genes, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), size = 2) +
  scale_color_manual(values = c("gray41" ,"red4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  geom_text_repel(data = subset(genes, (padj < 0.05 & abs(log2FoldChange) >=2)),
                  aes(label = Gene_Symbol),segment.color = 'transparent') +
  
  theme(legend.position="right")+
  ggtitle("Volcano plot")+
  xlab("log2(FoldChange)") + ylab("-log10(padj") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 


ggsave(filename = paste0(Gene_Type,"VolcanoPlot_10_20_2019.pdf"), width = 12, height = 8, dpi = 800)





###################################################

####miRNA



diffExp.DESeq2_miRNA <- read.csv("diffExp_DESeq_miRNA_ALL_10_19_2018.txt",header = TRUE)
rownames(diffExp.DESeq2_miRNA) <- diffExp.DESeq2_miRNA$X


#diffExp.DESeq2$Gene_Symbol <- TCGA.Harmonized$Gene_Symbol[(match(rownames(diffExp.DESeq2),TCGA.Harmonized$ENSB_ID))]
#diffExp.DESeq2$Gene_status <- TCGA.Harmonized$Gene_status[(match(rownames(diffExp.DESeq2),TCGA.Harmonized$ENSB_ID))]

#all_diff <- na.omit(diffExp.DESeq2) # remove genes with no annotation


genes <- diffExp.DESeq2_miRNA[,-9]  # diff_Exp_genes 

genes <- na.omit(genes)

#Gene_Type <- "protein_coding"
Gene_Type <- "miRNA"
#genes <- all_diff[all_diff$Gene_status== Gene_Type,]

x.cut=2;y.cut=0.05

Significance <- ifelse(genes$log2FoldChange >= x.cut & genes$padj < y.cut, "Upregulated", ifelse(genes$log2FoldChange <= -x.cut & genes$padj < y.cut, "Downregulated", "Not significant"))

ggplot(genes, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), size = 2) +
  scale_color_manual(values = c("gray41" ,"red4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  geom_text_repel(data = subset(genes, (padj < 0.05 & abs(log2FoldChange) >=2)),
                  aes(label = X),segment.color = 'transparent') +
  
  theme(legend.position="right")+
  ggtitle("Volcano plot")+
  xlab("log2(FoldChange)") + ylab("-log10(padj") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 


ggsave(filename = paste0(Gene_Type,"VolcanoPlot_10_20_2019.pdf"), width = 12, height = 8, dpi = 800)


#########################

# Diff Methylation volcano plot


meth <- read.csv("DiffMeth_Total_CpG_No_chrXYM.txt",header = TRUE)


Gene_Type <- "methylation"


#genes <- meth
meth1 <- meth[meth$adj.P.Val < 0.07,]
genes <-meth1


list_cpg <- c("cg03544320","cg03692651","cg07915921","cg15811515","cg19717586","cg22620090","cg22674699","cg22784954",
          "cg22797031","cg00446123","cg02772121","cg07388969","cg11201447","cg13860281","cg14931884", 
          "cg16389901","cg20765408","cg20852851")

genes_subset <- genes[genes$X %in% list_cpg,] 


x.cut=0.2;y.cut=0.005

Significance <- ifelse(genes$MeanDiff >= x.cut & genes$adj.P.Val < y.cut, "Upregulated", ifelse(genes$MeanDiff <= -x.cut & genes$adj.P.Val < y.cut, "Downregulated", "Not significant"))

ggplot(genes, aes(x = MeanDiff, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significance), size = 1) +
  scale_color_manual(values = c( "blue4","gray41" ,"red4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  
  geom_text_repel(data = genes_subset,
                  aes(label = X)) +
  
  
  theme(legend.position="right")+
  ggtitle("Volcano plot")+
  xlab("Delta Beta Value") + ylab("-log10(padj") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 


ggsave(filename = paste0(Gene_Type,"VolcanoPlot_10_20_2019.pdf"), width = 12, height = 8, dpi = 800)











