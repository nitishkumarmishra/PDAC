####################################################
##### Created By Nitish Mishra (March 7th 2019)
##### This is R code for multivariate analysis by using tidyverse approach.
##### For this analysis I used R package "survivalAnalysis" 
##### End of this code I also added multiple univariate analysis.
####################################################

setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
load("PAAD_Methylation_Oct-10.RData")
load("Survival_For_All_RANASeq.RData")
rm(list=setdiff(ls(), c("PAAD.Clinical", "BMIQ.Meth", "PAAD.FPKM.UQ.protein.ductal","PAAD.FPKM.UQ.lncRNA.ductal", "PAAD.miRNA.RPM.ductal")))
Clinical <- read.table("Clinical_data_for_manuscript.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = TRUE)

BMIQ.Meth.Tumor <- BMIQ.Meth[,grep("-01A", colnames(BMIQ.Meth))]
miRNA <- PAAD.miRNA.RPM.ductal[c("hsa-mir-196b"),]
Protein <- PAAD.FPKM.UQ.protein.ductal[c("ASPM", "B3GNT3", "BMF", "CD300LB", "CD68", "CENPF", "DEPDC1B", "DMBT1", "DTL", "ERCC6L", "FAM111B", "HIST1H2BC", "HIST1H2BJ", "HIST1H3H", "KIF4A", "NEK2", "RASSF4"),]
#c("ASPM", "B3GNT3", "BMF", "CD300LB", "CD68", "CENPF", "DEPDC1B", "DMBT1", "DTL", "ERCC6L", "FAM111B", "HIST1H2BC", "HIST1H2BJ", "HIST1H3H", "KIF4A", "NEK2", "RASSF4", "hsa-mir-196b")
Meth <- BMIQ.Meth.Tumor[c("cg02587316", "cg18630667", "cg05020604", "cg03234186"),]


colnames(Meth) <- substr(colnames(Meth), 1, 12)
Protein$Gene_Symbol <- NULL
colnames(Protein) <- substr(colnames(Protein), 1, 12)
colnames(miRNA) <- substr(colnames(miRNA), 1, 12)
rownames(Clinical) <- Clinical$Patient.ID


Meth <- t(Meth); miRNA <- t(miRNA); Protein <- t(Protein)
Protein <- as.data.frame(Protein); miRNA <- as.data.frame(miRNA); Meth <- as.data.frame(Meth)
Meth$Patient.ID <- rownames(Meth)
Protein$Patient.ID <- rownames(Protein)
miRNA$Patient.ID <- rownames(miRNA)

####################################################
##R library for tidyverse based survival analysis ##
library(plyr)
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)

Join <- join_all(list(Clinical, Protein, miRNA, Meth), by='Patient.ID', type = 'full')
Survival_data <- Join[-c(1:7,154),] ## Remove samples which have missing miRNA and Protein data

####################################################
######### Multivariate survival analysis ###########

Survival_data <- Survival_data %>% 
  mutate(Status = if_else(Patients.Vital.Status=='Alive', 1, 2))

colnames(Survival_data) <- gsub("hsa-mir-196b", "Mir196b", colnames(Survival_data))
covariate_names <- c(Diagnosis.Age="Age at Diagnosis",
                     Sex="Sex",
                     ASPM="ASPM",                         
                     B3GNT3="B3GNT3",                       
                     BMF="BMF",                         
                     CD300LB="CD300LB",                     
                     CD68="CD68",                         
                     CENPF="CENPF",                        
                     DEPDC1B="DEPDC1B",                      
                     DMBT1="DMBT1",                        
                     DTL="DTL",                          
                     ERCC6L="ERCC6L",                      
                     FAM111B="FAM111B",                      
                     HIST1H2BC="HIST1H2BC",                    
                     HIST1H2BJ="HIST1H2BJ",                    
                     HIST1H3H="HIST1H3H",                     
                     KIF4A="KIF4A",                       
                     NEK2="NEK2",                         
                     RASSF4="RASSF4",                       
                     Mir196b="Mir 196b",                 
                     cg02587316="cg02587316",                   
                     cg18630667="cg18630667",                   
                     cg05020604="cg05020604",                   
                     cg03234186="cg03234186"   
                     )

Survival_data %>%
  mutate(sex=rename_factor(Sex, `1` = "Male", `2` = "Female")) %>%
  analyse_multivariate(vars(Overall.Survival.Months, Status),
                       covariates = vars(Sex, ASPM,   B3GNT3, BMF, CD300LB,                     
                                         CD68, CENPF, DEPDC1B, DMBT1, DTL,                          
                                         ERCC6L, FAM111B, HIST1H2BC, HIST1H2BJ,                    
                                         HIST1H3H, KIF4A, NEK2, RASSF4, Mir196b,                 
                                         cg02587316, cg18630667, cg05020604, cg03234186),
                       covariate_name_dict = covariate_names) ->
  results

print(results)
forest_plot(results)

####################################################
### Multivariate analysis of only DNA methylation ##
####################################################
Survival_data %>%
  mutate(sex=rename_factor(Sex, `1` = "Male", `2` = "Female")) %>%
  analyse_multivariate(vars(Overall.Survival.Months, Status),
                       covariates = vars(cg02587316, cg18630667, cg05020604, cg03234186),
                       covariate_name_dict = covariate_names) ->
  results
print(results)
forest_plot(results)

####################################################
########### Multiple univariate analysis ###########
####################################################
df <- Survival_data
map(vars(Sex, ASPM,   B3GNT3, BMF, CD300LB,                     
         CD68, CENPF, DEPDC1B, DMBT1, DTL,                          
         ERCC6L, FAM111B, HIST1H2BC, HIST1H2BJ,                    
         HIST1H3H, KIF4A, NEK2, RASSF4, Mir196b,                 
         cg02587316, cg18630667, cg05020604, cg03234186), function(by)
{
  analyse_multivariate(df,
                       vars(Overall.Survival.Months, Status),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="OS"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 10))

####################################################
######## Save the workspace after analysis #########
save.image("Multivariate_Survival_Analysis.RData")
####################################################