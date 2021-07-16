library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")

diffMeth <- read.csv("DiffMeth_No_chrXYM_P.val_0.01_UCSC.txt", header = TRUE, sep = ",", row.names = 2)
Total.CpG <- read.csv("DiffMeth_Total_CpG_No_chrXYM.txt",header = TRUE, sep = ",", row.names = 1)

sigCpGs <- rownames(diffMeth)
gst <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(Total.CpG), collection="GO")
gst.KEGG <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(Total.CpG), collection="KEGG")


diffMeth.Hyper <- diffMeth[diffMeth$MeanDiff >0 ,]
diffMeth.Hypo <- diffMeth[diffMeth$MeanDiff <0 ,]

gst.hyper <- gometh(sig.cpg=rownames(diffMeth.Hyper), all.cpg=rownames(Total.CpG), collection="GO")
gst.KEGG.hyper <- gometh(sig.cpg=rownames(diffMeth.Hyper), all.cpg=rownames(Total.CpG), collection="KEGG")

gst.hyper <- gometh(sig.cpg=rownames(diffMeth.Hyper), all.cpg=rownames(Total.CpG), collection="GO")
gst.KEGG.hyper <- gometh(sig.cpg=rownames(diffMeth.Hyper), all.cpg=rownames(Total.CpG), collection="KEGG")


gst.hypo <- gometh(sig.cpg=rownames(diffMeth.Hypo), all.cpg=rownames(Total.CpG), collection="GO")
gst.KEGG.hypo <- gometh(sig.cpg=rownames(diffMeth.Hypo), all.cpg=rownames(Total.CpG), collection="KEGG")

spearman_Negative <- read.csv("spearman_Negative_Protein.txt", header = TRUE, sep = ",")
spearman_Negative.CpG <- unique(spearman_Negative$Probe)

dmCpG.neg <- intersect(spearman_Negative.CpG, rownames(diffMeth.Hyper))
gst.GO.dmCpG.neg <- gometh(sig.cpg=dmCpG.neg, all.cpg=rownames(Total.CpG), collection="GO")
gst.KEGG.dmCpG.neg <- gometh(sig.cpg=dmCpG.neg, all.cpg=rownames(Total.CpG), collection="KEGG")

save.image("missMethyl_KEGG.RData")
