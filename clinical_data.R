load("PAAD.RData")

TCGA_PAAD <- read.csv("study_view_clinical_data.txt", header = TRUE, sep = "\t")
TCGA_PAAD_select <- subset(TCGA_PAAD, select = c(Patient.ID, Sample.ID, Neoplasm.Histologic.Type.Name, Person.Neoplasm.Status, Patient.s.Vital.Status, Disease.Free.Status,Diagnosis.Age,Overall.Survival..Months.,Sex,Ethnicity.Category,Race.Category ))
TCGA_PAAD_select.ductal <- TCGA_PAAD_select[grep(pattern = "Pancreas-Adenocarcinoma Ductal", TCGA_PAAD_select$Neoplasm.Histologic.Type.Name),]


geneExp <- PAAD.HTSeq.Counts # RNA
miRNAexp <- PAAD.miRNA.Counts # miRNA
meth_data <- PAAD.Meth




TCGA_PAAD_select.ductal$GeneExpression <- ifelse(match(TCGA_PAAD_select.ductal$Patient.ID,(substr(colnames(geneExp),1,12))),"Yes","No")

TCGA_PAAD_select.ductal$miRNA <- ifelse(match(TCGA_PAAD_select.ductal$Patient.ID,(substr(colnames(miRNAexp),1,12))),"Yes","No")


TCGA_PAAD_select.ductal$Methylation <- ifelse(match(TCGA_PAAD_select.ductal$Patient.ID,(substr(colnames(meth_data),1,12))),"Yes","No")

write.csv(TCGA_PAAD_select.ductal, file ="Clinical_data_for_manuscript.txt", row.names = FALSE)
