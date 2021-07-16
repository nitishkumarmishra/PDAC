library(limma)
setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
load("PAAD.RData")
rm(list=setdiff(ls(), c("PAAD.Clinical", "PAAD.Meth")))
TCGA_PAAD <- read.csv("study_view_clinical_data.txt", header = TRUE, sep = "\t")
TCGA_PAAD_select <- subset(TCGA_PAAD, select = c(Patient.ID, Sample.ID, Neoplasm.Histologic.Type.Name, Person.Neoplasm.Status, Patient.s.Vital.Status, Disease.Free.Status, Tumor.Other.Histologic.Subtype, Overall.Survival..Months.))
TCGA_PAAD_select.ductal <- TCGA_PAAD_select[grep(pattern = "Pancreas-Adenocarcinoma Ductal", TCGA_PAAD_select$Neoplasm.Histologic.Type.Name),]

PAAD.Meth.ductal <- PAAD.Meth[,substr(colnames(PAAD.Meth), 1,12)%in%TCGA_PAAD_select.ductal$Patient.ID]
#table(substr(colnames(PAAD.Meth[,substr(colnames(PAAD.Meth), 1,12)%in%TCGA_PAAD_select.ductal$Patient.ID]), 14, 16))
cancerID <- grep("01A", substr(colnames(PAAD.Meth.ductal), 14, 16))
normalID <- grep("11A", substr(colnames(PAAD.Meth.ductal), 14, 16))
PAAD.Meth.ductal.cancer <- PAAD.Meth.ductal[,cancerID]
PAAD.Meth.ductal.normal <- PAAD.Meth.ductal[, normalID]


hm450.manifest <- read.csv("../../Pancreatic Endocrine/hm450.hg38.manifest.tsv.gz", header = TRUE, sep = "\t", row.names = 5)
hm450.manifest.select <- subset(hm450.manifest, select = c(designType))
hm450.manifest.select$type <- ifelse(hm450.manifest.select$designType=="II", 2, 1)
hm450.manifest.select.Meth <- hm450.manifest.select[rownames(Meth),]

#annotation <- merge(annotation, hm450.manifest.select, by="ProbeID")
#annotation.mask <- annotation[annotation$MASK_general==FALSE,]

numNAs <- rowSums(is.na(PAAD.Meth.ductal.cancer))
PAAD.Meth.ductal.cancer <- PAAD.Meth.ductal.cancer[!(numNAs > ncol(PAAD.Meth.ductal.cancer)*0.25),]
PAAD.Meth.ductal.normal <- PAAD.Meth.ductal.normal[!(numNAs > ncol(PAAD.Meth.ductal.normal)*0.25),]

PAAD.Meth.ductal.cancer <- PAAD.Meth.ductal.cancer[intersect(rownames(PAAD.Meth.ductal.cancer), rownames(PAAD.Meth.ductal.normal)),]
PAAD.Meth.ductal.normal <- PAAD.Meth.ductal.normal[intersect(rownames(PAAD.Meth.ductal.cancer), rownames(PAAD.Meth.ductal.normal)),]

## Imputation od tumor and normal sample separately doesn't work properly. So I used merged data for imputation.
# impute <- impute::impute.knn(as.matrix(PAAD.Meth.ductal.cancer), k = 15, rowmax = 0.25, colmax = 0.25, rng.seed = 12345)
# PAAD.Meth.ductal.cancer.impute <- impute$data
# 
# impute <- impute::impute.knn(as.matrix(PAAD.Meth.ductal.normal), k = 15, rowmax = 0.25, colmax = 0.25, rng.seed = 12345)
# PAAD.Meth.ductal.normal.impute <- impute$data
# 
# Meth <- cbind(PAAD.Meth.ductal.cancer.impute, PAAD.Meth.ductal.normal.impute)
# #myNorm <- ChAMP::champ.norm(beta = Meth, method =  "BMIQ", plotBMIQ = TRUE, arraytype = "450K", resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
# 
# 
# m1 <- matrix(0, nrow(Meth), ncol(Meth))
# colnames(m1) <- colnames(Meth)
# for (i in 1:ncol(Meth)) {
#   beta <- Meth[,i]
#   Meth.BMIQ <- wateRmelon::BMIQ(beta, design.v = hm450.manifest.select.Meth$type, plots = FALSE)
#   m1[,i] <- Meth.BMIQ$nbeta
#   
# }
# rownames(m1) <- rownames(Meth)
## BMIQ tumor and normal separately is not working great. So used BMIQ for full data (tumor and normal)
PAAD.Meth.ductal.cancer <- PAAD.Meth.ductal.cancer[intersect(rownames(PAAD.Meth.ductal.cancer), rownames(PAAD.Meth.ductal.normal)),]
PAAD.Meth.ductal.normal <- PAAD.Meth.ductal.normal[intersect(rownames(PAAD.Meth.ductal.cancer), rownames(PAAD.Meth.ductal.normal)),]

Meth <- cbind(PAAD.Meth.ductal.cancer, PAAD.Meth.ductal.normal)

impute <- impute::impute.knn(as.matrix(Meth), k = 15, rowmax = 0.25, colmax = 0.25, rng.seed = 12345)
Meth <- impute$data

m1 <- matrix(0, nrow(Meth), ncol(Meth))
colnames(m1) <- colnames(Meth)
for (i in 1:ncol(Meth)) {
  beta <- Meth[,i]
  Meth.BMIQ <- wateRmelon::BMIQ(beta, design.v = hm450.manifest.select.Meth$type, plots = FALSE)
  m1[,i] <- Meth.BMIQ$nbeta
  
}
rownames(m1) <- rownames(Meth)

BMIQ.Meth <- m1
BMIQ.Meth <-  lumi::beta2m(BMIQ.Meth)
cancerID <- grep("01A", substr(colnames(BMIQ.Meth), 14, 16))
normalID <- grep("11A", substr(colnames(BMIQ.Meth), 14, 16))
Tumor.BMIQ <- Meth[,cancerID]
Normal.BMIQ <- Meth[, normalID]
BMIQ.Meth <- cbind(Tumor.BMIQ, Normal.BMIQ)
############ Differential Methylation #################

design <- model.matrix(~0 + factor(c(rep(2, length(cancerID)), rep(1, length(normalID)))))
colnames(design) <- c("Normal", "Tumor")
cont.matrix <- makeContrasts("Tumor-Normal", levels = design)
methFit = lmFit(BMIQ.Meth, design)
#cont.matrix = makeContrasts(TumorVsNormal = Tumor-Normal, levels = design)
methFit2 = contrasts.fit(methFit, cont.matrix)
methFit2 = eBayes(methFit2)
methRresults <- topTable(methFit2, adjust.method = "BH", number = nrow(BMIQ.Meth), sort.by = "p", resort.by="logFC")
t.tumor <- as.data.frame(cbind(rowMeans(Tumor.BMIQ), rowMeans(Normal.BMIQ), rowMeans(Tumor.BMIQ) - rowMeans(Normal.BMIQ), rowMeans(Tumor.BMIQ) / rowMeans(Normal.BMIQ), log2(rowMeans(Tumor.BMIQ) / rowMeans(Normal.BMIQ))))
colnames(t.tumor) <- c("Tumor", "Normal", "MeanDiff", "FoldChange", "log2FC")
limma.results <- merge(methRresults, t.tumor, by="row.names")
rownames(limma.results) <- limma.results$Row.names
limma.results$Row.names <- NULL

results.dms <- limma.results[which(limma.results$adj.P.Val <= 0.01 & abs(limma.results$MeanDiff) >= 0.20) ,]
results.dms <- results.dms[order(results.dms$MeanDiff, decreasing = TRUE),]


hm450.manifest.chromosome <- subset(hm450.manifest, select = c(CpG_chrm, CpG_beg, CpG_end, MASK_general))
results.merge <- merge(hm450.manifest.chromosome, limma.results, by="row.names")
rownames(results.merge) <- results.merge$Row.names; results.merge$Row.names <- NULL
#results.merge[, !results.merge$CpG_chrm %in% c("chrM",  "chrX",  "chrY")]
results.table <- results.merge[!as.character(results.merge$CpG_chrm) %in% c("chrX", "chrY", "chrM"),]
results.table <- results.table[results.table$MASK_general==FALSE,]
mark <- grep("rs", rownames(results.table))
results.table <- results.table[-mark,]
results.dms.selected <- results.table[which(results.table$adj.P.Val <= 0.01 & abs(results.table$MeanDiff) >= 0.20) ,]
write.csv(results.dms.selected, "DiffMeth_No_chrXYM_P.val_0.01.txt")
write.csv(results.table, "DiffMeth_Total_CpG_No_chrXYM.txt")
#######################################################

results.dms.selected <- results.table[which(results.table$adj.P.Val <= 0.01 & abs(results.table$MeanDiff) >= 0.20) ,]

# Source survival plot code 
source("plot_surv_function.R")


unlink("Survival_Plot/KM_Plot_*") 
results <- plot_surv(dir = "Survival_Plot", clinical_patient = PAAD.Clinical, dataMETH = BMIQ.Meth, Genelist=rownames(results.dms.selected), Survresult = FALSE, ThreshTop = 0.5, ThreshDown = 0.3, PercentUp = 10, PercentDown = 10, p.cut = 0.05)

#unlink("Survival_Plot/KM_Plot_*") 
results.DM.P.0.01 <- plot_surv(dir = "Survival_Plot", clinical_patient = PAAD.Clinical, dataMETH = BMIQ.Meth, Genelist=rownames(results), Survresult = TRUE, ThreshTop = 0.5, ThreshDown = 0.3, PercentUp = 10, PercentDown = 10, p.cut = 0.01)

# unlink("Survival_Plot/KM_Plot_*") 
# results.DM.P.0.01 <- plot_surv(dir = "Survival_Plot", clinical_patient = PAAD.Clinical, dataMETH = BMIQ.Meth, Genelist=rownames(results.DM), Survresult = TRUE, ThreshTop = 0.5, ThreshDown = 0.3, PercentUp = 10, PercentDown = 10, p.cut = 0.01)

#######################################################
save.image("PAAD_Methylation_Oct-10.RData")
