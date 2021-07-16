library(dplyr)
library(ggpubr)


TCGA_PAAD <- read.csv("study_view_clinical_data.txt", header = TRUE, sep = "\t")
TCGA_PAAD_select <- subset(TCGA_PAAD, select = c(Patient.ID, Sample.ID, Neoplasm.Histologic.Type.Name, Person.Neoplasm.Status, Patient.s.Vital.Status, Disease.Free.Status, Tumor.Other.Histologic.Subtype, Overall.Survival..Months.))
TCGA_PAAD_select.ductal <- TCGA_PAAD_select[grep(pattern = "Pancreas-Adenocarcinoma Ductal", TCGA_PAAD_select$Neoplasm.Histologic.Type.Name),]

# Annotation File
TCGA.Harmonized <- read.table("TCGA_GDC_Harmonided_Uniq_GTF.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)
colnames(TCGA.Harmonized) <- c("Chr", "Reference", "Transcript", "Start", "End", "Strand", "Transcript_ID", "Gene_Symbol")


load("PAAD.RData")

# Expression data
geneExp <- PAAD.FPKM.UQ 
geneExp <- log2(geneExp + 1)

# Methylation data 
meth_data <- PAAD.Meth
 
# miRNA
miRNA_data <- PAAD.miRNA.Counts # miRNA
miRNA_data <- log2(miRNA_data + 1)
 
rm(list= ls()[!(ls() %in% c('TCGA_PAAD_select.ductal','geneExp','meth_data','miRNA_data','TCGA.Harmonized','PAAD.FPKM.UQ', 'PAAD.Meth', 'PAAD.miRNA.Counts', 'PAAD.miRNA.RPM'))])


dim(geneExp) 
dim(meth_data)
dim(miRNA_data)

# for processing seperately
lincRNA_data <- geneExp 

protein_coding <- c("B3GNT3","DMBT1")
miRNA <- c("hsa-mir-196b","hsa-mir-135a-2","hsa-mir-125a","hsa-mir-3200","hsa-mir-196a-2","hsa-mir-196a-1")
cpg <- c("cg03234186","cg02587316")
lincRNA <- c("PVT1","RP11-54H7.4","LINC00941","CASC11", "RP11-1029J19.4","MIR205HG")




############### GeneExpression data ########################
# Select subset of Pancreatic adenocarcinoma ductal
geneExp <- geneExp[,substr(colnames(geneExp),1,12) %in% TCGA_PAAD_select.ductal$Patient.ID]
dim(geneExp)

cancerID_GE <- grep("01A", colnames(geneExp))
normalID_GE <- grep("11A", colnames(geneExp))
length(cancerID_GE)
length(normalID_GE)


geneExp.cancer <- geneExp[,cancerID_GE]
geneExp.normal <- geneExp[,normalID_GE]

geneExp <- cbind(geneExp.cancer, geneExp.normal)

# geneExp[dim(geneExp)[1]+1,cancerID_GE] <- "Tumor"
# geneExp[dim(geneExp)[1],normalID_GE] <- "Normal"
# rownames(geneExp)[3] <- "Type"

geneExp.select <- geneExp[match(protein_coding,(TCGA.Harmonized$Gene_Symbol[(match(rownames(geneExp),TCGA.Harmonized$Transcript_ID))])),]
# replace row names
rownames(geneExp.select) <- TCGA.Harmonized$Gene_Symbol[(match(rownames(geneExp.select),TCGA.Harmonized$Transcript_ID))]

geneExp.select <- as.data.frame(t(geneExp.select))
geneExp.select$Type <- c(rep("Tumor", length(cancerID_GE)), rep("Normal", length(normalID_GE)))



############## Redo for lincRNA ###############################


lncRNA.select <- geneExp[match(lincRNA,(TCGA.Harmonized$Gene_Symbol[(match(rownames(geneExp),TCGA.Harmonized$Transcript_ID))])),]
# replace row names
rownames(lncRNA.select) <- TCGA.Harmonized$Gene_Symbol[(match(rownames(lncRNA.select),TCGA.Harmonized$Transcript_ID))]

lncRNA.select <- as.data.frame(t(lncRNA.select))
lncRNA.select$Type <- c(rep("Tumor", length(cancerID_GE)), rep("Normal", length(normalID_GE)))



# # Select subset of Pancreatic adenocarcinoma ductal
# lincRNA_data <- lincRNA_data[,substr(colnames(lincRNA_data),1,12) %in% TCGA_PAAD_select.ductal$Patient.ID]
# dim(lincRNA_data)
## WORKS ONLY IF UNIQUE GENE NAME 
# lincRNA_data <- lincRNA_data[match(lincRNA,(TCGA.Harmonized$Gene_Symbol[(match(rownames(lincRNA_data),TCGA.Harmonized$Transcript_ID))])),]
# # replace row names
# rownames(lincRNA_data) <- TCGA.Harmonized$Gene_Symbol[(match(rownames(lincRNA_data),TCGA.Harmonized$Transcript_ID))]

# cancerID_linc <- grep("01A", colnames(lincRNA_data))
# normalID_linc <- grep("11A", colnames(lincRNA_data))
# length(cancerID_linc)
# length(normalID_linc)
# 
# lincRNA_data[dim(lincRNA_data)[1]+1,cancerID_linc] <- "Tumor"
# lincRNA_data[dim(lincRNA_data)[1],normalID_linc] <- "Normal"
# rownames(lincRNA_data)[5] <- "Type"





################ Methylation data #############################
# Select subset of Pancreatic adenocarcinoma ductal
meth_data <- meth_data[,substr(colnames(meth_data),1,12) %in% TCGA_PAAD_select.ductal$Patient.ID]
dim(meth_data)


cancerID_MT <- grep("01A", colnames(meth_data))
normalID_MT <- grep("11A", colnames(meth_data))
length(cancerID_MT)
length(normalID_MT)


meth_data.cancer <- meth_data[,cancerID_MT]
meth_data.normal <- meth_data[,normalID_MT]

meth_data <- cbind(meth_data.cancer, meth_data.normal)

# 
# meth_data[dim(meth_data)[1]+1,cancerID_MT] <- "Tumor"
# meth_data[dim(meth_data)[1],normalID_MT] <- "Normal"
# rownames(meth_data)[3] <- "Type"

# ONLY IF UNIQUE NAMES
meth_data <- meth_data[match(cpg,rownames(meth_data)),]

meth_data <- as.data.frame(t(meth_data))
meth_data$Type <- c(rep("Tumor", length(cancerID_MT)), rep("Normal", length(normalID_MT)))



################ miRNA data #############################
# Select subset of Pancreatic adenocarcinoma ductal
miRNA_data <- miRNA_data[,substr(colnames(miRNA_data),1,12) %in% TCGA_PAAD_select.ductal$Patient.ID]
dim(miRNA_data)

cancerID_miRNA <- grep("01A", colnames(miRNA_data))
normalID_miRNA <- grep("11A", colnames(miRNA_data))
length(cancerID_miRNA)
length(normalID_miRNA)

miRNA_data.cancer <- miRNA_data[,cancerID_miRNA]
miRNA_data.normal <- miRNA_data[,normalID_miRNA]

miRNA_data <- cbind(miRNA_data.cancer, miRNA_data.normal)


# ONLY IF UNIQUE NAMES
miRNA_data <- miRNA_data[match(miRNA,rownames(miRNA_data)),]


miRNA_data <- as.data.frame(t(miRNA_data))
miRNA_data$Type <- c(rep("Tumor", length(cancerID_miRNA)), rep("Normal", length(normalID_miRNA)))



# miRNA_data[dim(miRNA_data)[1]+1,cancerID_miRNA] <- "Tumor"
# miRNA_data[dim(miRNA_data)[1],normalID_miRNA] <- "Normal"
# rownames(miRNA_data)[5] <- "Type"


############# PLOTS ###########################
my_comparisons <- list(c("Tumor", "Normal"))



################ Protein coding #####################
# geneExp <- as.data.frame(t(geneExp))
# geneExp$B3GNT3 <- as.numeric(as.character(geneExp$B3GNT3))
# geneExp$DMBT1 <- as.numeric(as.character(geneExp$DMBT1))


## B3GNT3

P <- ggboxplot(geneExp.select, x = "Type", y = "B3GNT3",
               title = "B3GNT3", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))

P + #stat_compare_means(label.y = max(geneExp.select$B3GNT3)+1)
  stat_compare_means(comparisons = my_comparisons,  label.y = max(geneExp.select$B3GNT3)+0.2) 

ggsave("Boxplot_B3GNT3.pdf", width = 4, height = 5)




## DMBT1

P <- ggboxplot(geneExp.select, x = "Type", y = "DMBT1",
               title = "DMBT1", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))

P + #stat_compare_means(label.y = max(geneExp.select$DMBT1)+1)+
  stat_compare_means(comparisons = my_comparisons,  label.y = max(geneExp.select$DMBT1)+0.2) 


ggsave("Boxplot_DMBT1.pdf", width = 4, height = 5)



################## lincRNA #######################################

# lincRNA_data <- as.data.frame(t(lincRNA_data))
# lincRNA_data$PVT1 <- as.numeric(as.character(lincRNA_data$PVT1))
# lincRNA_data$`RP11-54H7.4` <- as.numeric(as.character(lincRNA_data$`RP11-54H7.4`))
# lincRNA_data$LINC00941 <- as.numeric(as.character(lincRNA_data$LINC00941))
# lincRNA_data$CASC11 <- as.numeric(as.character(lincRNA_data$CASC11))
# 

## "PVT1"

P <- ggboxplot(lncRNA.select, x = "Type", y = "PVT1",
               title = "PVT1", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(lncRNA.select$PVT1)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(lncRNA.select$PVT1)+0.2) 

ggsave("Boxplot_PVT1.pdf", width = 4, height = 5)


#"RP11-54H7.4"

P <- ggboxplot(lncRNA.select, x = "Type", y = "`RP11-54H7.4`",
               title = "RP11-54H7.4", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(lncRNA.select$`RP11-54H7.4`)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(lncRNA.select$`RP11-54H7.4`)+0.2) 

ggsave("Boxplot_RP11-54H7_4.pdf", width = 4, height = 5)



## LINC00941

P <- ggboxplot(lncRNA.select, x = "Type", y = "LINC00941",
               title = "LINC00941", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(lncRNA.select$LINC00941)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(lncRNA.select$LINC00941)+0.2) 

ggsave("Boxplot_LINC00941.pdf", width = 4, height = 5)


## CASC11

P <- ggboxplot(lncRNA.select, x = "Type", y = "CASC11",
               title = "CASC11", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(lncRNA.select$CASC11)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(lncRNA.select$CASC11)+0.2) 

ggsave("Boxplot_CASC11.pdf", width = 4, height = 5)


## MIR205HG

P <- ggboxplot(lncRNA.select, x = "Type", y = "MIR205HG",
               title = "MIR205HG", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(lncRNA.select$MIR205HG)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(lncRNA.select$MIR205HG)+0.2) 

ggsave("Boxplot_MIR205HG.pdf", width = 4, height = 5)



## RP11-1029J19.4

P <- ggboxplot(lncRNA.select, x = "Type", y = "`RP11-1029J19.4`",
               title = "RP11-1029J19.4", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(lncRNA.select$`RP11-1029J19.4`)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(lncRNA.select$`RP11-1029J19.4`)+0.2) 

ggsave("Boxplot_RP11-1029J19_4.pdf", width = 4, height = 5)




######################## miRNA ###################################

# miRNA_data <- as.data.frame(t(miRNA_data))
# miRNA_data$`hsa-mir-196b` <- as.numeric(as.character(miRNA_data$`hsa-mir-196b`))
# miRNA_data$`hsa-mir-135a-2` <- as.numeric(as.character(miRNA_data$`hsa-mir-135a-2`))
# miRNA_data$`hsa-mir-125a` <- as.numeric(as.character(miRNA_data$`hsa-mir-125a`))
# miRNA_data$`hsa-mir-3200` <- as.numeric(as.character(miRNA_data$`hsa-mir-3200`))



# `hsa-mir-196b`

P <- ggboxplot(miRNA_data, x = "Type", y = "`hsa-mir-196b`",
               title = "hsa-mir-196b", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(miRNA_data$`hsa-mir-196b`)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(miRNA_data$`hsa-mir-196b`)+0.2) 


ggsave("Boxplot_hsa-mir-196b.pdf", width = 4, height = 5)



# `hsa-mir-135a-2`
P <- ggboxplot(miRNA_data, x = "Type", y = "`hsa-mir-135a-2`",
               title = "hsa-mir-135a-2", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(miRNA_data$`hsa-mir-135a-2`)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(miRNA_data$`hsa-mir-135a-2`)+0.2) 


ggsave("Boxplot_hsa-mir-135a-2.pdf", width = 4, height = 5)


# `hsa-mir-125a`
P <- ggboxplot(miRNA_data, x = "Type", y = "`hsa-mir-125a`",
               title = "hsa-mir-125a", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(miRNA_data$`hsa-mir-125a`)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(miRNA_data$`hsa-mir-125a`)+0.2) 

ggsave("Boxplot_hsa-mir-125a.pdf", width = 4, height = 5)



# `hsa-mir-3200`

P <- ggboxplot(miRNA_data, x = "Type", y = "`hsa-mir-3200`",
               title = "hsa-mir-3200", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(miRNA_data$`hsa-mir-3200`)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(miRNA_data$`hsa-mir-3200`)+0.2) 


ggsave("Boxplot_hsa-mir-3200.pdf", width = 4, height = 5)




# `hsa-mir-196a-2`

P <- ggboxplot(miRNA_data, x = "Type", y = "`hsa-mir-196a-2`",
               title = "hsa-mir-196a-2", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(miRNA_data$`hsa-mir-196a-2`)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(miRNA_data$`hsa-mir-196a-2`)+0.2) 


ggsave("Boxplot_hsa-mir-196a-2.pdf", width = 4, height = 5)




# `hsa-mir-196a-1`

P <- ggboxplot(miRNA_data, x = "Type", y = "`hsa-mir-196a-1`",
               title = "hsa-mir-196a-1", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(miRNA_data$`hsa-mir-196a-1`)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(miRNA_data$`hsa-mir-196a-1`)+0.2) 


ggsave("Boxplot_hsa-mir-196a-1.pdf", width = 4, height = 5)






######################### CPG ##########################
#cpg <- c("cg03234186","cg02587316")


# meth_data <- as.data.frame(t(meth_data))
# meth_data$cg03234186 <- as.numeric(as.character(meth_data$cg03234186))
# meth_data$cg02587316 <- as.numeric(as.character(meth_data$cg02587316))


## cg03234186

P <- ggboxplot(meth_data, x = "Type", y = "cg03234186",
               title = "cg03234186", ylab = "Beta Value", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(meth_data$cg03234186)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(meth_data$cg03234186)+0.2) 


ggsave("Boxplot_cg03234186.pdf", width = 4, height = 5)




## cg02587316

P <- ggboxplot(meth_data, x = "Type", y = "cg02587316",
               title = "cg02587316", ylab = "Beta Value", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + #stat_compare_means(label.y = max(meth_data$cg02587316)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(meth_data$cg02587316)+0.2) 

ggsave("Boxplot_cg02587316.pdf", width = 4, height = 5)



#######################################################################
#######################################################################

























library(dplyr)
library(ggpubr)


TCGA_PAAD <- read.csv("study_view_clinical_data.txt", header = TRUE, sep = "\t")
TCGA_PAAD_select <- subset(TCGA_PAAD, select = c(Patient.ID, Sample.ID, Neoplasm.Histologic.Type.Name, Person.Neoplasm.Status, Patient.s.Vital.Status, Disease.Free.Status, Tumor.Other.Histologic.Subtype, Overall.Survival..Months.))
TCGA_PAAD_select.ductal <- TCGA_PAAD_select[grep(pattern = "Pancreas-Adenocarcinoma Ductal", TCGA_PAAD_select$Neoplasm.Histologic.Type.Name),]

# Annotation File
TCGA.Harmonized <- read.table("TCGA_GDC_Harmonided_Uniq_GTF.txt",header=T,sep="\t",na.strings="NA",row.names=NULL,quote="",comment.char="", stringsAsFactors = FALSE, check.names = FALSE)
colnames(TCGA.Harmonized) <- c("Chr", "Reference", "Transcript", "Start", "End", "Strand", "Transcript_ID", "Gene_Symbol")


load("PAAD.RData")

# Expression data
geneExp <- PAAD.FPKM.UQ 
geneExp <- log2(geneExp + 1)

# Methylation data 
meth_data <- PAAD.Meth

# miRNA
miRNA_data <- PAAD.miRNA.Counts # miRNA
miRNA_data <- log2(miRNA_data + 1)

rm(list= ls()[!(ls() %in% c('TCGA_PAAD_select.ductal','geneExp','meth_data','miRNA_data','TCGA.Harmonized','PAAD.FPKM.UQ', 'PAAD.Meth', 'PAAD.miRNA.Counts', 'PAAD.miRNA.RPM'))])

###########################################################################
lincRNA_data <- geneExp 

protein_coding <- c("B3GNT3","DMBT1")
miRNA <- c("hsa-mir-196b","hsa-mir-135a-2","hsa-mir-125a","hsa-mir-3200", "hsa-mir-196a-2")


cpg <- c("cg03234186","cg02587316")
lincRNA <- c("PVT1","RP11-54H7.4","LINC00941","CASC11", "MIR205HG")

###########################################################################
geneExp <- geneExp[,substr(colnames(geneExp),1,12) %in% TCGA_PAAD_select.ductal$Patient.ID]

cancerID_GE <- grep("01A", colnames(geneExp))
normalID_GE <- grep("11A", colnames(geneExp))

geneExp.cancer <- geneExp[,cancerID_GE]
geneExp.normal <- geneExp[,normalID_GE]

geneExp <- cbind(geneExp.cancer, geneExp.normal)

##########################################################################

geneExp.select <- geneExp[match(protein_coding,(TCGA.Harmonized$Gene_Symbol[(match(rownames(geneExp),TCGA.Harmonized$Transcript_ID))])),]
# replace row names
rownames(geneExp.select) <- TCGA.Harmonized$Gene_Symbol[(match(rownames(geneExp.select),TCGA.Harmonized$Transcript_ID))]

geneExp.select <- as.data.frame(t(geneExp.select))
geneExp.select$Type <- c(rep("Tumor", length(cancerID_GE)), rep("Normal", length(normalID_GE)))

###########################################################################
############# PLOTS ###########################
my_comparisons <- list(c("Tumor", "Normal"))

################ Protein coding #####################
## B3GNT3

P <- ggboxplot(geneExp.select, x = "Type", y = "B3GNT3",
               title = "B3GNT3", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + stat_compare_means(label.y = max(geneExp.select$B3GNT3)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(geneExp.select$B3GNT3)+0.2) 

ggsave("Boxplot_B3GNT3.pdf", width = 4, height = 5)

###################################
lncRNA.select <- geneExp[match(lincRNA,(TCGA.Harmonized$Gene_Symbol[(match(rownames(geneExp),TCGA.Harmonized$Transcript_ID))])),]
# replace row names
rownames(lncRNA.select) <- TCGA.Harmonized$Gene_Symbol[(match(rownames(lncRNA.select),TCGA.Harmonized$Transcript_ID))]

lncRNA.select <- as.data.frame(t(lncRNA.select))
lncRNA.select$Type <- c(rep("Tumor", length(cancerID_GE)), rep("Normal", length(normalID_GE)))

#############################################

P <- ggboxplot(lncRNA.select, x = "Type", y = "MIR205HG",
               title = "MIR205HG", ylab = "log2(FPKM)", xlab = "Sample",
               color = "Type", palette = c("red4","blue4"),add = "jitter",
               add.params = list(size = 0.1, jitter = 0.2))
P + stat_compare_means(label.y = max(lncRNA.select$MIR205HG)+1)+ 
  stat_compare_means(comparisons = my_comparisons,  label.y = max(lncRNA.select$MIR205HG)+0.2) 

ggsave("Boxplot_MIR205HG.pdf", width = 4, height = 5)


