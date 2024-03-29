library(dplyr)
library(reshape2)
library(ROCR)
setwd("D:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
load("PAAD_Methylation_Oct-10.RData")
rm(list=setdiff(ls(), c("m1")))

TCGA_PAAD <- read.csv("study_view_clinical_data.txt", header = TRUE, sep = "\t")
TCGA_PAAD_select <- subset(TCGA_PAAD, select = c(Patient.ID, Sample.ID, Neoplasm.Histologic.Type.Name, Person.Neoplasm.Status, Patient.s.Vital.Status, Disease.Free.Status, Tumor.Other.Histologic.Subtype, Overall.Survival..Months.))
TCGA_PAAD_select.ductal <- TCGA_PAAD_select[grep(pattern = "Pancreas-Adenocarcinoma Ductal", TCGA_PAAD_select$Neoplasm.Histologic.Type.Name),]

ind_Exp <- as.data.frame(m1)
ind_Exp$Index <- rownames(ind_Exp)
ind_Exp$TargetID <- rownames(ind_Exp)
ind_Exp <- ind_Exp %>% filter(!is.na(Index)) %>% droplevels()
rownames(ind_Exp) <- ind_Exp$TargetID
ind_Exp <- ind_Exp[rownames(diffMeth),]
dat <- ind_Exp %>% select(Index, TargetID, contains("TCGA")) %>% melt(id = c("Index", "TargetID")) %>% mutate(group = ifelse(grepl("-01A", variable), 1, 0)) %>% split(ind_Exp$TargetID)
#### Now code for R function ########### 
#It will write two functio se_auc ci_auc
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
if (!dir.exists("DiffMeth_ROC_Plot")) dir.create("DiffMeth_ROC_Plot")   ### change directory name here
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
    file_name = paste0("DiffMeth_ROC_Plot/roc_", name, ".png")    ### change directory name here (the one you created above.
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
write.csv(auc_df, "AUC_table_DiffMeth_ROC_Plot.csv", row.names = F)   ### change filename to save AUC values.
write.csv(auc_df, "AUC_table_DiffMeth_ROC_Plot.txt", row.names = F)   ### change filename to save AUC values.
save.image("DiffMeth_ROC_Plot.RData")

#save.image("Survival_AUC_Final.RData")
