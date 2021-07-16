library(dplyr)
library(pacman)
p_load('DESeq2',"edgeR","limma")



load("PAAD.RData")

TCGA_PAAD <- read.csv("study_view_clinical_data.txt", header = TRUE, sep = "\t")
TCGA_PAAD_select <- subset(TCGA_PAAD, select = c(Patient.ID, Sample.ID, Neoplasm.Histologic.Type.Name, Person.Neoplasm.Status, Patient.s.Vital.Status, Disease.Free.Status,Diagnosis.Age,Overall.Survival..Months.,Sex,Ethnicity.Category,Race.Category ))
TCGA_PAAD_select.ductal <- TCGA_PAAD_select[grep(pattern = "Pancreas-Adenocarcinoma Ductal", TCGA_PAAD_select$Neoplasm.Histologic.Type.Name),]




TCGA.Harmonized <- read.table("TCGA_GDC_Harmonided_Uniq_GTF.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

colnames(TCGA.Harmonized) <- c("Chr", "Reference", "Transcript", "Start", "End", "Strand", "Transcript_ID", "Gene_Symbol")



geneExp <- PAAD.HTSeq.Counts # RNA
#geneExp <- PAAD.miRNA.Counts # miRNA
##geneExp <- PAAD.miRNA.RPM


################ DIFF GENE EXPRESSION ######################################

# Select subset of Pancreatic adenocarcinoma ductal
geneExp <- geneExp[,substr(colnames(geneExp),1,12) %in% TCGA_PAAD_select.ductal$Patient.ID]
dim(geneExp)

cancerID <- grep("01A", colnames(geneExp))
normalID <- grep("11A", colnames(geneExp))

length(cancerID)
length(normalID)

cancerExp <- geneExp[,cancerID]
cnts <- cancerExp[rowSums(cancerExp==0)< ncol(cancerExp)*0.20,] ## Remove all gene which have 25% zero's
keep <- rowSums(cpm(cnts)>1) >= ncol(cnts)*0.20 #### Select only genes which have have CPM > 1 for >=50% samples
cnts <- cnts[keep,]
gene <- rownames(cnts)
cnts <- cbind(geneExp[gene,cancerID], geneExp[gene,normalID])
design <- model.matrix(~0 + factor(c(rep(2, length(cancerID)), rep(1, length(normalID)))))
colnames(design) <- c("Normal", "Tumor")
cont.matrix <- makeContrasts("Tumor-Normal", levels = design)
factors <- factor(c(rep("Tumor", length(cancerID)), rep("Normal", length(normalID))))
cnts.dgelist <- DGEList(cnts, group=factors)
cnf <- calcNormFactors(cnts.dgelist, method = "TMM")

######################## DESeq ###########################################
dds <- DESeqDataSetFromMatrix(cnf,DataFrame(factors), ~ factors)
dds <- DESeq(dds)
resmiRNA <- results(dds, pAdjustMethod = "BH")
resmiRNAOrdered <- resmiRNA[order(resmiRNA$log2FoldChange,decreasing = TRUE),]
#resmiRNASig <- subset(resmiRNAOrdered,(padj < 0.05 & abs(resmiRNAOrdered$log2FoldChange) >= 2))

resmiRNASig <- subset(resmiRNAOrdered,(padj < 0.05))

diffExp.DESeq2 <- as.data.frame(resmiRNASig)
diffExp.DESeq2$foldChange <- 2^diffExp.DESeq2$log2FoldChange
diffExp.DESeq2$Gene_Symbol <- TCGA.Harmonized$Gene_Symbol[(match(rownames(diffExp.DESeq2),TCGA.Harmonized$Transcript_ID))]

dim(diffExp.DESeq2)

write.csv(diffExp.DESeq2,file="diffExp_DESeq_ALL_Genes_10_19_2018.txt",row.names = TRUE)


# # ###### Limma #############
# # 
# #eneExp <- log10(geneExp + 1) # log2 values
# fit <- lmFit(geneExp,design)
# fit2 <- contrasts.fit(fit,cont.matrix)
# fit2 <- eBayes(fit2)
# adj.method = "BH"; adj.pval = 0.05; raw.pval = 0.05; logFC = 1
# aradeger <- topTable(fit2, adjust.method = adj.method, number = length(fit2), sort.by='logFC')
# aradeger1 <- data.frame(aradeger[aradeger$adj.P.Val < adj.pval & aradeger$P.Value < raw.pval & abs(aradeger$logFC) >=logFC,])
# diffExp.voom <- aradeger1
# tmp <- decideTests(fit2, p.value = 0.05, lfc = 1)
# summary(tmp)## We will get number of up and downregulated genes
# 
# 
# #####################################################
# #########################GSEA###########################
# 
# 
# 
# v <- voom(cnf, design, plot = FALSE)
# fit <- lmFit(v, design)
# fit2 <- contrasts.fit(fit, cont.matrix)
# fit2 <- eBayes(fit2)
# adj.method = "BH"; adj.pval = 0.05; raw.pval = 0.05; logFC = 1
# aradeger <- topTable(fit2, adjust.method = adj.method, number = length(fit2), sort.by='logFC')
# aradeger1 <- data.frame(aradeger[aradeger$adj.P.Val < adj.pval & aradeger$P.Value < raw.pval & abs(aradeger$logFC) >=logFC,])
# diffExp.voom <- aradeger1
# tmp <- decideTests(fit2, p.value = 0.05, lfc = 1)
# summary(tmp)## We will get number of up and downregulated genes
# 
# 
# 






library(ggplot2)
library(fgsea)


###############################
#rank <- read.csv("gsea.rnk.txt", header = FALSE, sep = "\t")

# Use all the genes without cut off
all_diff <- diffExp.DESeq2
all_diff <- na.omit(all_diff)
tmp <- all_diff$log2FoldChange
names(tmp) <- all_diff$Gene_Symbol
tmp1 <- tmp[order(tmp, decreasing = TRUE)]
##################################################
c2 <- qusage::read.gmt("c2.cp.v6.2.symbols.gmt")

#tmp1 <- tmp[order(tmp, decreasing = TRUE)]

fgseaRes <- fgsea(pathways = c2, 
                  stats = tmp1,
                  minSize=5,
                  maxSize=500,
                  nperm=1000)

plotEnrichment(c2[["KEGG_GLYCOLYSIS_GLUCONEOGENESIS"]],tmp)+
  labs(title="KEGG_GLYCOLYSIS_GLUCONEOGENESIS") 




plotEnrichment(c2[["KEGG_PENTOSE_PHOSPHATE_PATHWAY"]],tmp)+
  labs(title="KEGG_CITRATE_CYCLE_TCA_CYCL") + stat_smooth(se = FALSE)



### KEGG pathway:

c5<- clusterProfiler::read.gmt("c2.cp.v6.2.symbols.gmt")
c2 <- qusage::read.gmt("c2.cp.kegg.v6.2.symbols.gmt")





# library(PPInfer)
# GSEA.barplot(fgseaRes, category = 'pathway', score = 'NES',
#              pvalue = 'pval', sort = 'NES', decreasing = TRUE) +
#   scale_fill_continuous(low = 'red', high = 'green')



dev.print(pdf, 'gsea_barplot.pdf', width = 10, height = 10)

#####################


library("pathview")
hsa04110 <- pathview(gene.data  = tmp1,
                     pathway.id = "hsa00010",
                     species    = "hsa",
                     limit      = list(gene=max(abs(tmp1)), cpd=1))


eg = bitr(names(tmp1), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")




