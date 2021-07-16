setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
load("PAAD.RData")
rm(list=setdiff(ls(), c("PAAD.Clinical", "PAAD.Meth", "PAAD.FPKM", "PAAD.FPKM.UQ", "PAAD.miRNA.RPM")))

library(dplyr)
library(reshape2)
library(ROCR)

####################### Expression data prepration ##########################
GDC.Harmonized <- read.csv("TCGA_GDC_Harmonided_Uniq_Status_GTF.txt", header = TRUE, sep = "\t", row.names = 1)
GDC.Harmonized.protein <- GDC.Harmonized[GDC.Harmonized$Gene_status=="protein_coding",]
GDC.Harmonized.lincRNA <- GDC.Harmonized[GDC.Harmonized$Gene_status=="lincRNA",]
#GDC.Harmonized.processed_pseudogene <- GDC.Harmonized[GDC.Harmonized$Gene_status=="processed_pseudogene",]

TCGA_PAAD <- read.csv("study_view_clinical_data.txt", header = TRUE, sep = "\t")
TCGA_PAAD_select <- subset(TCGA_PAAD, select = c(Patient.ID, Sample.ID, Neoplasm.Histologic.Type.Name, Person.Neoplasm.Status, Patient.s.Vital.Status, Disease.Free.Status, Tumor.Other.Histologic.Subtype, Overall.Survival..Months.))
TCGA_PAAD_select.ductal <- TCGA_PAAD_select[grep(pattern = "Pancreas-Adenocarcinoma Ductal", TCGA_PAAD_select$Neoplasm.Histologic.Type.Name),]

PAAD.FPKM.UQ.lncRNA <- PAAD.FPKM.UQ[rownames(GDC.Harmonized.lincRNA),]
PAAD.FPKM.UQ.protein <- PAAD.FPKM.UQ[rownames(GDC.Harmonized.protein),]

#cancerID <- grep("01A", substr(colnames(PAAD.FPKM.UQ.lncRNA), 14, 16))
PAAD.FPKM.UQ.lncRNA <- PAAD.FPKM.UQ.lncRNA

#cancerID <- grep("01A", substr(colnames(PAAD.FPKM.UQ.protein), 14, 16))
PAAD.FPKM.UQ.protein <- PAAD.FPKM.UQ.protein

#cancerID <- grep("01A", substr(colnames(PAAD.miRNA.RPM), 14, 16))
PAAD.miRNA.RPM.cancer <- PAAD.miRNA.RPM


colnames(PAAD.FPKM.UQ.lncRNA) <- substr(colnames(PAAD.FPKM.UQ.lncRNA), 1, 15)
colnames(PAAD.FPKM.UQ.protein) <- substr(colnames(PAAD.FPKM.UQ.protein), 1, 15)
colnames(PAAD.miRNA.RPM.cancer) <- substr(colnames(PAAD.miRNA.RPM.cancer), 1, 15)




PAAD.FPKM.UQ.lncRNA.ductal <-  PAAD.FPKM.UQ.lncRNA[,substr(colnames(PAAD.FPKM.UQ.lncRNA), 1, 12)%in% substr(TCGA_PAAD_select.ductal$Sample.ID, 1, 12)]
PAAD.FPKM.UQ.protein.ductal <-  PAAD.FPKM.UQ.protein[,substr(colnames(PAAD.FPKM.UQ.protein), 1, 12)%in% substr(TCGA_PAAD_select.ductal$Sample.ID, 1, 12)]
PAAD.miRNA.RPM.ductal <- PAAD.miRNA.RPM.cancer[,substr(colnames(PAAD.miRNA.RPM.cancer), 1, 12)%in% substr(TCGA_PAAD_select.ductal$Sample.ID, 1, 12)]

## Remove genes which have more than 20% zeros
PAAD.FPKM.UQ.lncRNA.ductal <- PAAD.FPKM.UQ.lncRNA.ductal[apply(PAAD.FPKM.UQ.lncRNA.ductal == 0, 1, sum) <= ncol(PAAD.FPKM.UQ.lncRNA.ductal)*0.20, ]
PAAD.FPKM.UQ.protein.ductal <- PAAD.FPKM.UQ.protein.ductal[apply(PAAD.FPKM.UQ.protein.ductal == 0, 1, sum) <= ncol(PAAD.FPKM.UQ.protein.ductal)*0.20, ]
PAAD.miRNA.RPM.ductal <- PAAD.miRNA.RPM.ductal[apply(PAAD.miRNA.RPM.ductal ==0, 1, sum) <= ncol(PAAD.miRNA.RPM.ductal)*0.20,]

PAAD.FPKM.UQ.lncRNA.ductal <- log2(PAAD.FPKM.UQ.lncRNA.ductal + 1)
PAAD.FPKM.UQ.protein.ductal <- log2(PAAD.FPKM.UQ.protein.ductal + 1)
PAAD.miRNA.RPM.ductal <- log2(PAAD.miRNA.RPM.ductal + 1)

lncRNA.symbol <- subset(GDC.Harmonized.lincRNA, select = c(Gene_Symbol))
protein.symbol <- subset(GDC.Harmonized.protein, select = c(Gene_Symbol))

PAAD.FPKM.UQ.lncRNA.ductal <- merge(PAAD.FPKM.UQ.lncRNA.ductal, lncRNA.symbol, by="row.names")
rownames(PAAD.FPKM.UQ.lncRNA.ductal) <- PAAD.FPKM.UQ.lncRNA.ductal$Gene_Symbol
PAAD.FPKM.UQ.lncRNA.ductal$Row.names <- NULL

PAAD.FPKM.UQ.protein.ductal <- merge(PAAD.FPKM.UQ.protein.ductal, protein.symbol, by="row.names")
rownames(PAAD.FPKM.UQ.protein.ductal) <- PAAD.FPKM.UQ.protein.ductal$Gene_Symbol
PAAD.FPKM.UQ.protein.ductal$Row.names <- NULL

PAAD.FPKM.UQ.lncRNA.ductal$Gene_Symbol <- NULL
PAAD.FPKM.UQ.protein.ductal$Gene_Symbol <- NULL

############################################
rna <- PAAD.FPKM.UQ.lncRNA.ductal
rna <- (rna-min(rna))/(max(rna)-min(rna))
ind_Exp <- as.data.frame(rna)
#ind_Exp <- ind_Exp[common_names,]
ind_Exp$Index <- rownames(ind_Exp)
ind_Exp$TargetID <- rownames(ind_Exp)
#rownames(ind_Exp) <- ind_Exp$TargetID
ind_Exp <- ind_Exp %>% filter(!is.na(Index)) %>% droplevels()
#
dat <- ind_Exp %>% select(Index, TargetID, contains("TCGA")) %>% melt(id = c("Index", "TargetID")) %>% mutate(group = ifelse(grepl("-01", variable), 1, 0)) %>% split(ind_Exp$TargetID)

###############################

se_auc <- function(auc, n1, n2) { #n1 is positive, n2 is negative
  q1 <- auc / (2 - auc)
  q2 <- (2 * auc ^ 2) / (1 + auc)
  sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc ^ 2) + (n2 - 1) * (q2 - auc ^ 2)) / (n1 * n2))
}

ci_auc <- function(auc, n1, n2, level = 0.95) {
  ci <- auc + c(-1, 1) * qnorm((1 + level) / 2) * se_auc(auc, n1, n2)
  ci[1] <- ifelse(ci[1] < 0, 0, ci[1])
  ci[2] <- ifelse(ci[2] > 1, 1, ci[2])
  ci }
save_roc_plot <- function(roc_obj) {}

############ R code for AUC ###########
aucs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "auc")@y.values)
rocs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "tpr", "fpr"))
# Create directory (if it doesn't exist)
if (!dir.exists("lncRNA_ROC_Plot")) dir.create("lncRNA_ROC_Plot")   ### change directory name here
# Loop through ROCs to make plots.
ns <- unname(table(dat[[1]]$group))
auc_df <- data.frame(name = character(length(dat)), 
                     AUC = numeric(length(dat)), 
                     CI_lower = numeric(length(dat)), 
                     CI_upper = numeric(length(dat)), stringsAsFactors = F)
for (i in 1:length(rocs)) {
  roc <- rocs[[i]]
  auc <- aucs[[i]]
  name <- names(aucs[i])
  attributes(roc)$alpha.values[[1]][1] <- 1
  # Flip results if AUC < 0.5
  ## Why he is flipping AUC values
  if (auc < 0.5) {
    auc <- 1 - unlist(auc)
    temp <- roc@x.values
    roc@x.values <- roc@y.values
    roc@y.values <- temp
  }
  ci <- ci_auc(unlist(auc), ns[1], ns[2])
  ci_chr <- paste0(round(ci, 2), collapse = ", ") 
  auc_df[i, "name"] <- name
  auc_df[i, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
  if (auc > 0.5) {
    file_name = paste0("lncRNA_ROC_Plot/roc_", name, ".png")    ### change directory name here (the one you created above.
    png(file_name)
    #plot(roc, main = name, col = "red")
    plot(roc, main = name, colorize=T)  
    abline(0, 1, lty = 2)
    ## We can change position by changing 0.65, 0.20 etc.
    ## These range are based on ROC plot X and Y axis value [1,1] range
    ## AUC start from X aix=0.8, Y axis = 0.13, similarly CI start from position 0.8 on X and 0.08 on Y
    ## rect(0.62, 0.03, 0.98, 0.18); It start on X=0.62 end X=0.98 (go from 0.66 to 0.98 on X)
    ## On Y axis start Y= 0.03 go to 0.18
    text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
    text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
    rect(0.62, 0.03, 0.98, 0.18)
    dev.off()
  }
}
auc_df <- arrange(auc_df, desc(AUC)) %>% filter(AUC > 0.05)
write.csv(auc_df, "AUC_table_lncRNA.csv", row.names = F)   ### change filename to save AUC values.
write.csv(auc_df, "AUC_table_lncRNA.txt", row.names = F)   ### change filename to save AUC values.

############################################################################################
rna <- PAAD.FPKM.UQ.protein.ductal
rna <- (rna-min(rna))/(max(rna)-min(rna))
ind_Exp <- as.data.frame(rna)
#ind_Exp <- ind_Exp[common_names,]
ind_Exp$Index <- rownames(ind_Exp)
ind_Exp$TargetID <- rownames(ind_Exp)
#rownames(ind_Exp) <- ind_Exp$TargetID
ind_Exp <- ind_Exp %>% filter(!is.na(Index)) %>% droplevels()
#
dat <- ind_Exp %>% select(Index, TargetID, contains("TCGA")) %>% melt(id = c("Index", "TargetID")) %>% mutate(group = ifelse(grepl("-01", variable), 1, 0)) %>% split(ind_Exp$TargetID)


se_auc <- function(auc, n1, n2) { #n1 is positive, n2 is negative
  q1 <- auc / (2 - auc)
  q2 <- (2 * auc ^ 2) / (1 + auc)
  sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc ^ 2) + (n2 - 1) * (q2 - auc ^ 2)) / (n1 * n2))
}

ci_auc <- function(auc, n1, n2, level = 0.95) {
  ci <- auc + c(-1, 1) * qnorm((1 + level) / 2) * se_auc(auc, n1, n2)
  ci[1] <- ifelse(ci[1] < 0, 0, ci[1])
  ci[2] <- ifelse(ci[2] > 1, 1, ci[2])
  ci }
save_roc_plot <- function(roc_obj) {}

############ R code for AUC ###########
aucs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "auc")@y.values)
rocs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "tpr", "fpr"))
# Create directory (if it doesn't exist)
if (!dir.exists("Protein_ROC_Plot")) dir.create("Protein_ROC_Plot")   ### change directory name here
# Loop through ROCs to make plots.
ns <- unname(table(dat[[1]]$group))
auc_df <- data.frame(name = character(length(dat)), 
                     AUC = numeric(length(dat)), 
                     CI_lower = numeric(length(dat)), 
                     CI_upper = numeric(length(dat)), stringsAsFactors = F)
for (i in 1:length(rocs)) {
  roc <- rocs[[i]]
  auc <- aucs[[i]]
  name <- names(aucs[i])
  attributes(roc)$alpha.values[[1]][1] <- 1
  # Flip results if AUC < 0.5
  ## Why he is flipping AUC values
  if (auc < 0.5) {
    auc <- 1 - unlist(auc)
    temp <- roc@x.values
    roc@x.values <- roc@y.values
    roc@y.values <- temp
  }
  ci <- ci_auc(unlist(auc), ns[1], ns[2])
  ci_chr <- paste0(round(ci, 2), collapse = ", ") 
  auc_df[i, "name"] <- name
  auc_df[i, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
  if (auc > 0.5) {
    file_name = paste0("Protein_ROC_Plot/roc_", name, ".png")    ### change directory name here (the one you created above.
    png(file_name)
    #plot(roc, main = name, col = "red")
    plot(roc, main = name, colorize=T)  
    abline(0, 1, lty = 2)
    ## We can change position by changing 0.65, 0.20 etc.
    ## These range are based on ROC plot X and Y axis value [1,1] range
    ## AUC start from X aix=0.8, Y axis = 0.13, similarly CI start from position 0.8 on X and 0.08 on Y
    ## rect(0.62, 0.03, 0.98, 0.18); It start on X=0.62 end X=0.98 (go from 0.66 to 0.98 on X)
    ## On Y axis start Y= 0.03 go to 0.18
    text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
    text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
    rect(0.62, 0.03, 0.98, 0.18)
    dev.off()
  }
}
auc_df <- arrange(auc_df, desc(AUC)) %>% filter(AUC > 0.05)
write.csv(auc_df, "AUC_table_Protein.csv", row.names = F)   ### change filename to save AUC values.
write.csv(auc_df, "AUC_table_Protein.txt", row.names = F)   ### change filename to save AUC values.


###########################################################
rna <- PAAD.miRNA.RPM.ductal
rna <- (rna-min(rna))/(max(rna)-min(rna))
ind_Exp <- as.data.frame(rna)
#ind_Exp <- ind_Exp[common_names,]
ind_Exp$Index <- rownames(ind_Exp)
ind_Exp$TargetID <- rownames(ind_Exp)
#rownames(ind_Exp) <- ind_Exp$TargetID
ind_Exp <- ind_Exp %>% filter(!is.na(Index)) %>% droplevels()
#
dat <- ind_Exp %>% select(Index, TargetID, contains("TCGA")) %>% melt(id = c("Index", "TargetID")) %>% mutate(group = ifelse(grepl("-01", variable), 1, 0)) %>% split(ind_Exp$TargetID)


se_auc <- function(auc, n1, n2) { #n1 is positive, n2 is negative
  q1 <- auc / (2 - auc)
  q2 <- (2 * auc ^ 2) / (1 + auc)
  sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc ^ 2) + (n2 - 1) * (q2 - auc ^ 2)) / (n1 * n2))
}

ci_auc <- function(auc, n1, n2, level = 0.95) {
  ci <- auc + c(-1, 1) * qnorm((1 + level) / 2) * se_auc(auc, n1, n2)
  ci[1] <- ifelse(ci[1] < 0, 0, ci[1])
  ci[2] <- ifelse(ci[2] > 1, 1, ci[2])
  ci }
save_roc_plot <- function(roc_obj) {}

############ R code for AUC ###########
aucs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "auc")@y.values)
rocs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "tpr", "fpr"))
# Create directory (if it doesn't exist)
if (!dir.exists("miRNA_ROC_Plot")) dir.create("miRNA_ROC_Plot")   ### change directory name here
# Loop through ROCs to make plots.
ns <- unname(table(dat[[1]]$group))
auc_df <- data.frame(name = character(length(dat)), 
                     AUC = numeric(length(dat)), 
                     CI_lower = numeric(length(dat)), 
                     CI_upper = numeric(length(dat)), stringsAsFactors = F)
for (i in 1:length(rocs)) {
  roc <- rocs[[i]]
  auc <- aucs[[i]]
  name <- names(aucs[i])
  attributes(roc)$alpha.values[[1]][1] <- 1
  # Flip results if AUC < 0.5
  ## Why he is flipping AUC values
  if (auc < 0.5) {
    auc <- 1 - unlist(auc)
    temp <- roc@x.values
    roc@x.values <- roc@y.values
    roc@y.values <- temp
  }
  ci <- ci_auc(unlist(auc), ns[1], ns[2])
  ci_chr <- paste0(round(ci, 2), collapse = ", ") 
  auc_df[i, "name"] <- name
  auc_df[i, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
  if (auc > 0.5) {
    file_name = paste0("miRNA_ROC_Plot/roc_", name, ".png")    ### change directory name here (the one you created above.
    png(file_name)
    #plot(roc, main = name, col = "red")
    plot(roc, main = name, colorize=T)  
    abline(0, 1, lty = 2)
    ## We can change position by changing 0.65, 0.20 etc.
    ## These range are based on ROC plot X and Y axis value [1,1] range
    ## AUC start from X aix=0.8, Y axis = 0.13, similarly CI start from position 0.8 on X and 0.08 on Y
    ## rect(0.62, 0.03, 0.98, 0.18); It start on X=0.62 end X=0.98 (go from 0.66 to 0.98 on X)
    ## On Y axis start Y= 0.03 go to 0.18
    text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
    text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
    rect(0.62, 0.03, 0.98, 0.18)
    dev.off()
  }
}
auc_df <- arrange(auc_df, desc(AUC)) %>% filter(AUC > 0.05)
write.csv(auc_df, "AUC_table_miRNA.csv", row.names = F)   ### change filename to save AUC values.
write.csv(auc_df, "AUC_table_miRNA.txt", row.names = F)   ### change filename to save AUC values.

###########################
save.image("ROC_Plot_RNASeq.RData")
