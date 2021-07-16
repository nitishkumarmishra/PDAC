library(data.table)
library(cowplot)
library(ggplot2)
library(ggpubr)

# exp_meth_protein_coding <- read.csv("exp_meth_protein_coding_10_19_2018.txt", header = TRUE, check.names = FALSE)
# 
# 
# 
# 
# df <- data.frame(matrix(ncol = 3, nrow = 137))  ## 137*3
# colnames(df) <- c("exp","meth","Probe")
# datalist =list()
# 
# 
# gene_name <- c("cg20911165")
#   
#   Probe = which(gene_name == exp_meth_protein_coding$ProbeID)
#   
#   
#   
#   datalist[[1]] <- exp_meth_protein_coding[Probe,]
#   
# 
#   
#   
# 
# dataset <- rbind(datalist[[1]])
# 
# df$exp <- as.numeric(as.character(c(dataset[1,6:142])))
# df$meth <- as.numeric(as.character(c(dataset[1,143:279])))
# df$Probe <- as.factor(c(rep("cg20911165",137)))
# 
# 
# # Use box plot as marginal plots
# ggscatterhist(
#   df, x = "exp", y = "meth", size = 3, alpha = 0.6,
#   color = "red4",
#   margin.plot = "boxplot",
#   ggtheme = theme_bw()
# )
# 
# 
# 
# 
# ## 
# 
# 
# 
# 
# 
# 
# 
# TCGA_PAAD <- read.csv("study_view_clinical_data.txt", header = TRUE, sep = "\t")
# TCGA_PAAD_select <- subset(TCGA_PAAD, select = c(Patient.ID, Sample.ID, Neoplasm.Histologic.Type.Name, Person.Neoplasm.Status, Patient.s.Vital.Status, Disease.Free.Status, Tumor.Other.Histologic.Subtype, Overall.Survival..Months.))
# TCGA_PAAD_select.ductal <- TCGA_PAAD_select[grep(pattern = "Pancreas-Adenocarcinoma Ductal", TCGA_PAAD_select$Neoplasm.Histologic.Type.Name),]
# 
# 
# load("PAAD.RData")
# 
# 
# 
# 
# # Expression data
# geneExp <- PAAD.FPKM.UQ 
# geneExp <- log2(geneExp + 1)
# 
# 
# 
# # Methylation data
# 
# All_BMIQ<-read.csv("Methylation_All_CpGs.csv",header = TRUE,check.names = FALSE)
# 
# # Select subset of Pancreatic adenocarcinoma ductal
# geneExp <- geneExp[,substr(colnames(geneExp),1,12) %in% TCGA_PAAD_select.ductal$Patient.ID]
# dim(geneExp)
# 
# 
# cancerID <- grep("01A", colnames(geneExp))
# normalID <- grep("11A", colnames(geneExp))
# 
# length(cancerID)
# length(normalID)
# 
# cancerExp <- geneExp[,normalID]
# cnts <- cancerExp[rowSums(cancerExp==0)< ncol(cancerExp)*0.20,] ## Remove all gene which have 25% zero's
# keep <- rowSums(cpm(cnts)>1) >= ncol(cnts)*0.20 #### Select only genes which have have CPM > 1 for >=50% samples
# cnts <- cnts[keep,]
# gene <- rownames(cnts)
# 
# 
# ################ MERGE METHYLATION AND EXPRESSION DATA ################################
# 
# hm450.anno <-read.table("hm450.hg38.manifest.tsv.gz",header=T,sep="\t",na.strings="",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)
# 
# rownames(hm450.anno) <- hm450.anno$probeID
# 
# #hm450.anno.unmasked <- hm450.anno[hm450.anno$MASK.general == 'FALSE',]
# hm450.anno.unmasked <- hm450.anno %>% filter(MASK_general %in% "FALSE")
# 
# rm(hm450.anno)
# 
# hm450.anno.epig <-read.table("hm450.manifest.EpigeneticallySilence_Sort.tsv",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)
# 
# ## hm450.manifest.EpigeneticallySilence_Sort.tsv is the file for CpG sites with distance from TSS. I generated this file from GDC harmonized DNA methylation file for genCode22
# 
# hm450.anno.epig <- hm450.anno.epig[which(hm450.anno.epig$Position_to_TSS!='.'),]
# 
# hm450.anno.epig.tss1500 <- hm450.anno.epig[abs(as.numeric(hm450.anno.epig$Position_to_TSS)) <= 1500,]
# 
# hm450.anno.epig.tss1500.unmasked <- hm450.anno.epig.tss1500[hm450.anno.epig.tss1500$ProbeID%in%hm450.anno.unmasked$probeID,]
# 
# hm450.tss1500 <- hm450.anno.epig.tss1500.unmasked[,c("ProbeID", "Transcript_ID", "Gene_Symbol","Position_to_TSS", "Chromosome", "Gene_Type" )]
# 
# colnames(hm450.tss1500) <- gsub("Chromosome","Chr", colnames(hm450.tss1500))
# 
# #TCGA.Harmonized <- read.table("TCGA_GDC_Harmonided_RNASeq-GFT.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)
# 
# TCGA.Harmonized <- read.table("TCGA_GDC_Harmonided_Uniq_GTF.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)
# 
# colnames(TCGA.Harmonized) <- c("Chr", "Reference", "Transcript", "Start", "End", "Strand", "Transcript_ID", "Gene_Symbol")
# 
# masked <- which(TCGA.Harmonized$Chr%in% c("chrX", "chrY", "chrM"))
# 
# TCGA.Har.Masked <- TCGA.Harmonized[-masked,]
# 
# TCGA.Har.Masked <- TCGA.Har.Masked[,c("Chr", "Transcript_ID", "Gene_Symbol")]
# 
# masked <- which(hm450.tss1500$Chr%in% c("chrX", "chrY", "chrM"))
# 
# hm450.tss1500.Masked <- hm450.tss1500[-masked,]
# 
# hm450.tss1500.Masked.Uniq <- unique(hm450.tss1500.Masked[,c("ProbeID", "Gene_Symbol", "Chr","Gene_Type")]) ## Unique gene and probe within +/- 1500 from TSS
# 
# # Add Ensembl Gene ID
# hm450.tss1500.Masked.Uniq$gene_id <- TCGA.Harmonized$Transcript_ID[match(hm450.tss1500.Masked.Uniq$Gene_Symbol,TCGA.Harmonized$Gene_Symbol)]
# 
# 
# # Tumor expression data
# cnts_normal<- geneExp[gene,normalID]
# #cnts_normal = cnts_normal[,substr(colnames(cnts_normal),1,12) %in%  TCGA_PAAD_select.ductal$Sample.ID]
# colnames(cnts_normal) <- substr(colnames(cnts_normal),1,15)
# 
# # Methylation data
# #Tumor_BMIQ <- read.csv("Tumor_BMIQ.txt",header = TRUE, check.names = FALSE)
# rownames(All_BMIQ) <- All_BMIQ[,1]
# All_BMIQ <- All_BMIQ[,-1]
# normalID2 <- grep("11A", colnames(All_BMIQ)) # Samples in Methylation Data
# meth_data <- All_BMIQ[,normalID2]
# 
# #meth_data = meth_data[,substr(colnames(meth_data),1,15) %in% TCGA_PAAD_select.ductal$Sample.ID]
# colnames(meth_data) <- substr(colnames(meth_data),1,15)
# 
# 
# common_patients <-Reduce(intersect,list(colnames(cnts_normal),colnames(meth_data)))
# 
# 
# # Take forward only tumor samples that are common in both
# cnts_normal <- cnts_normal[,common_patients]
# meth_data <- meth_data[,common_patients]
# 
# 
# merge1 <-cbind(hm450.tss1500.Masked.Uniq, cnts_normal[match(hm450.tss1500.Masked.Uniq$gene_id,rownames(cnts_normal)),])
# merge1 <- na.omit(merge1)
# dim(merge1)
# 
# # Merge methylation data with probe_id
# merge2 <-cbind(merge1, meth_data[match(merge1$ProbeID,rownames(meth_data)),])
# merge2<- na.omit(merge2)
# dim(merge2)
# #write.csv(merge2,file="exp_meth_ductal_data.txt", row.names = FALSE)
# 
# 
# #Protein Coding genes 
# merge2_protein_coding <- merge2[merge2$Gene_Type=="protein_coding",]
# merge2_protein_coding <- na.omit(merge2_protein_coding)
# dim(merge2_protein_coding)
# write.csv(merge2_protein_coding,file = "Normal_exp_meth_protein_coding.txt",row.names = FALSE)
# 
# 
# 
# 

exp_meth_protein_coding <- read.csv("exp_meth_protein_coding_11_6_2018.txt", header = TRUE, check.names = FALSE)
normal_exp_meth <- read.csv("Normal_exp_meth_protein_coding.txt")






for (gene_name in c("cg20911165","cg03609102")) {
  
  df <- data.frame(matrix(ncol = 3, nrow = 149))  ## 137+ 3  146 + 3
  colnames(df) <- c("exp","meth","Type")
  datalist =list()



Probe_tumor = which(gene_name == exp_meth_protein_coding$ProbeID)
Probe_normal = which(gene_name == normal_exp_meth$ProbeID)


datalist[[1]] <- exp_meth_protein_coding[Probe_tumor,]
normal_data <- normal_exp_meth[Probe_normal,]




dataset <- rbind(datalist[[1]])

df$exp[1:146] <- as.numeric(as.character(c(dataset[1,6:151])))
df$meth[1:146] <- as.numeric(as.character(c(dataset[1,152:297])))
#df$Probe[1:137] <- as.factor(c(rep("cg20911165",137)))
df$Type[1:146] <- (c(rep("Tumor",146)))


# add Normal data # 3 normal samples
df$exp[147:149] <- as.numeric(as.character(c(normal_data[1,6:8])))
df$meth[147:149] <- as.numeric(as.character(c(normal_data[1,9:11])))
df$Type[147:149] <- (c(rep("Normal",3)))

df$Type <- as.factor(df$Type)
  
  
# Use box plot as marginal plots
ggscatterhist(
  df, x = "exp", y = "meth", size = 2, alpha = 0.6,
  color = "Type",
  margin.plot = "boxplot", margin.params = list(fill = "Type"),
  palette = c("blue4","red4"),cor.coef = TRUE, conf.int = TRUE, add = "reg.line", add.params = list(color = "gray4",fill = NA), # Customize reg. line
  cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"), xlab = "log2(FPKM-UQ)", ylab = "DNA Methylation (Beta Value)", title = gene_name,
  ggtheme = theme_bw()
) 


ggsave(filename = paste0("muc5b_",gene_name,".pdf"), width = 10, height = 8, dpi = 800)

}
