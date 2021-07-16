library(dplyr)

TCGA_PAAD <- read.csv("study_view_clinical_data.txt", header = TRUE, sep = "\t")
TCGA_PAAD_select <- subset(TCGA_PAAD, select = c(Patient.ID, Sample.ID, Neoplasm.Histologic.Type.Name, Person.Neoplasm.Status, Patient.s.Vital.Status, Disease.Free.Status, Tumor.Other.Histologic.Subtype, Overall.Survival..Months.))
TCGA_PAAD_select.ductal <- TCGA_PAAD_select[grep(pattern = "Pancreas-Adenocarcinoma Ductal", TCGA_PAAD_select$Neoplasm.Histologic.Type.Name),]


load("PAAD.RData")

# Expression data
geneExp <- PAAD.FPKM.UQ 
geneExp <- log2(geneExp + 1)

# Methylation data 
Tumor_BMIQ <- read.csv("Tumor_BMIQ.txt",header = TRUE, check.names = FALSE)


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


################ MERGE METHYLATION AND EXPRESSION DATA ################################

hm450.anno <-read.table("hm450.hg38.manifest.tsv.gz",header=T,sep="\t",na.strings="",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

rownames(hm450.anno) <- hm450.anno$probeID

#hm450.anno.unmasked <- hm450.anno[hm450.anno$MASK.general == 'FALSE',]
hm450.anno.unmasked <- hm450.anno %>% filter(MASK_general %in% "FALSE")

rm(hm450.anno)

hm450.anno.epig <-read.table("hm450.manifest.EpigeneticallySilence_Sort.tsv",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

## hm450.manifest.EpigeneticallySilence_Sort.tsv is the file for CpG sites with distance from TSS. I generated this file from GDC harmonized DNA methylation file for genCode22

hm450.anno.epig <- hm450.anno.epig[which(hm450.anno.epig$Position_to_TSS!='.'),]

hm450.anno.epig.tss1500 <- hm450.anno.epig[abs(as.numeric(hm450.anno.epig$Position_to_TSS)) <= 1500,]

hm450.anno.epig.tss1500.unmasked <- hm450.anno.epig.tss1500[hm450.anno.epig.tss1500$ProbeID%in%hm450.anno.unmasked$probeID,]

hm450.tss1500 <- hm450.anno.epig.tss1500.unmasked[,c("ProbeID", "Transcript_ID", "Gene_Symbol","Position_to_TSS", "Chromosome", "Gene_Type" )]

colnames(hm450.tss1500) <- gsub("Chromosome","Chr", colnames(hm450.tss1500))

#TCGA.Harmonized <- read.table("TCGA_GDC_Harmonided_RNASeq-GFT.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

TCGA.Harmonized <- read.table("TCGA_GDC_Harmonided_Uniq_GTF.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)

colnames(TCGA.Harmonized) <- c("Chr", "Reference", "Transcript", "Start", "End", "Strand", "Transcript_ID", "Gene_Symbol")

masked <- which(TCGA.Harmonized$Chr%in% c("chrX", "chrY", "chrM"))

TCGA.Har.Masked <- TCGA.Harmonized[-masked,]

TCGA.Har.Masked <- TCGA.Har.Masked[,c("Chr", "Transcript_ID", "Gene_Symbol")]

masked <- which(hm450.tss1500$Chr%in% c("chrX", "chrY", "chrM"))

hm450.tss1500.Masked <- hm450.tss1500[-masked,]

hm450.tss1500.Masked.Uniq <- unique(hm450.tss1500.Masked[,c("ProbeID", "Gene_Symbol", "Chr","Gene_Type")]) ## Unique gene and probe within +/- 1500 from TSS

# Add Ensembl Gene ID
hm450.tss1500.Masked.Uniq$gene_id <- TCGA.Harmonized$Transcript_ID[match(hm450.tss1500.Masked.Uniq$Gene_Symbol,TCGA.Harmonized$Gene_Symbol)]



# Tumor expression data
cnts_tumor <- geneExp[gene,cancerID]
cnts_tumor = cnts_tumor[,substr(colnames(cnts_tumor),1,15) %in%  TCGA_PAAD_select.ductal$Sample.ID]
colnames(cnts_tumor) <- substr(colnames(cnts_tumor),1,15)

# Methylation data
#Tumor_BMIQ <- read.csv("Tumor_BMIQ.txt",header = TRUE, check.names = FALSE)
rownames(Tumor_BMIQ) <- Tumor_BMIQ[,1]
Tumor_BMIQ <- Tumor_BMIQ[,-1]
cancerID2 <- grep("01A", colnames(Tumor_BMIQ)) # Samples in Methylation Data
meth_data <- Tumor_BMIQ[,cancerID2]

meth_data = meth_data[,substr(colnames(meth_data),1,15) %in% TCGA_PAAD_select.ductal$Sample.ID]
colnames(meth_data) <- substr(colnames(meth_data),1,15)
dim(meth_data)

common_patients <-Reduce(intersect,list(colnames(cnts_tumor),colnames(meth_data)))
length(common_patients)

# Take forward only tumor samples that are common in both
cnts_tumor <- cnts_tumor[,common_patients]
meth_data <- meth_data[,common_patients]


merge1 <-cbind(hm450.tss1500.Masked.Uniq, cnts_tumor[match(hm450.tss1500.Masked.Uniq$gene_id,rownames(cnts_tumor)),])
merge1 <- na.omit(merge1)
dim(merge1)

# Merge methylation data with probe_id
merge2 <-cbind(merge1, meth_data[match(merge1$ProbeID,rownames(meth_data)),])
merge2<- na.omit(merge2)
dim(merge2)
##columns : 6 to 129 expression data 124 Cancer)
##columms :130:253 methylation data 124 cancer)
write.csv(merge2,file="exp_meth_ductal_data_11_6_2018.txt", row.names = FALSE)


#Protein Coding genes 
merge2_protein_coding <- merge2[merge2$Gene_Type=="protein_coding",]
merge2_protein_coding <- na.omit(merge2_protein_coding)
dim(merge2_protein_coding)
write.csv(merge2_protein_coding,file = "exp_meth_protein_coding_11_6_2018.txt",row.names = FALSE)

#lincRNA
merge2_lincRNA <- merge2[merge2$Gene_Type=="lincRNA",]
merge2_lincRNA <- na.omit(merge2_lincRNA)
dim(merge2_lincRNA)
write.csv(merge2_lincRNA,file = "exp_meth_lincRNA_coding_11_6_2018.txt",row.names = FALSE)

###################################################################################
###################################################################################


### Spearman correlation  
## Run this in parts : first for protein coding then for lincRNA
## Make 3 file, 1: Meth, 2: Expression, 3: Annotation

## Protein Coding
meth_file1 <- merge2_protein_coding[,152:297]
exp_file2 <- merge2_protein_coding[,6:151]
annotation_file3 <- merge2_protein_coding[,1:5]


##lincRNA
meth_file1 <- merge2_lincRNA[,152:297]
exp_file2 <- merge2_lincRNA[,6:151]
annotation_file3 <- merge2_lincRNA[,1:5]



ll <- mapply(function(x,y)cor.test(as.numeric(meth_file1[x,]),as.numeric(exp_file2[y,]), method = "spearman", alternative = "t"),
             1:nrow(meth_file1),
             1:nrow(exp_file2),
             SIMPLIFY=FALSE)

cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
spearman.p <- cbind(cor.value, p.value)
rownames(spearman.p) <- rownames(annotation_file3)
spearman.p <- data.frame(spearman.p)
spearman.p$gene <- annotation_file3$Gene_Symbol
spearman.p$gene_id <- annotation_file3$gene_id
spearman.p$Probe <- annotation_file3$ProbeID
spearman.p$chr <- annotation_file3$Chr
spearman.p <- spearman.p[!as.character(spearman.p$chr) %in% c("chrM","chrX", "chrY"),]
spearman.p$adjP <- p.adjust(spearman.p$p.value, method = c("BH"))
spearman.p1 <- spearman.p[which((spearman.p$p.value <= 0.01)&abs(spearman.p$cor.value) >=0.25),]
spearman.p1 <- spearman.p1[order(spearman.p1$gene),]
spearman.p1.neg <- spearman.p1[which(spearman.p1$cor.value < 0),]
spearman.p1.pos <- spearman.p1[which(spearman.p1$cor.value > 0),]

dim(spearman.p1.pos)

dim(spearman.p1.neg)
################################################################


write.csv(spearman.p1.neg,file="spearman_Negative_Protein_11_6_2018.txt",row.names = FALSE)
write.csv(spearman.p1.pos,file="spearman_Postive_Protein_11_6_2018.txt",row.names = FALSE)

write.csv(spearman.p1.neg,file="spearman_Negative_lincRNA_11_6_2018.txt",row.names = FALSE)
write.csv(spearman.p1.pos,file="spearman_Postive_lincRNA_11_6_2018.txt",row.names = FALSE)
