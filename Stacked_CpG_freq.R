library(RColorBrewer)

diffMeth <- read.csv("DiffMeth_No_chrXYM_P.val_0.01.txt", header = TRUE, sep = ",", row.names = 1)
Hyper <- diffMeth[diffMeth$MeanDiff >0, ]
Hypo <- diffMeth[diffMeth$MeanDiff < 0, ]
dmCpGs.Hyper <- table(Hyper$CpG_chrm)
dmCpGs.Hypo <- table(Hypo$CpG_chrm)
CpG.all <- cbind(dmCpGs.total, dmCpGs.Hyper, dmCpGs.Hypo)

Chromosome_Size <- read.csv2("Chromosome_Size.txt", header = TRUE, sep = "\t")
Chromosome_Size$Chr <- paste0("chr", Chromosome_Size$Chr)
Chromosome_Size$Mb <- Chromosome_Size$Size/1000000

rownames(Chromosome_Size) <- Chromosome_Size$Chr
Chromosome_Size <- cbind(Chromosome_Size, CpG.all)
Chromosome_Size$Name <- NULL
Chromosome_Size$Code <- NULL

Chromosome_Size$CpG.per.Mb <- round(Chromosome_Size$dmCpGs.total/Chromosome_Size$Mb, 3)
Chromosome_Size$Hyper.per.Mb <- round(Chromosome_Size$dmCpGs.Hyper/Chromosome_Size$Mb, 3)
Chromosome_Size$Hypo.per.Mb <- round(Chromosome_Size$dmCpGs.Hypo/Chromosome_Size$Mb, 3)
Chromosome_Size <- Chromosome_Size[order(Chromosome_Size$CpG.per.Mb, decreasing = TRUE),]

stacked.plot <- Chromosome_Size[,c("Hyper.per.Mb", "Hypo.per.Mb")]
rownames(stacked.plot) <- gsub("chr", "", rownames(stacked.plot))

colors <- c("red4", "blue4")
pdf("stacked_Meth_perMb.pdf", width = 8, height = 8)
par(mar=c(8,8,4,3))
#barplot(t(as.matrix(chr)), col = sequential, xlab = "Chromosome", ylab = "Methylation/Mb")
barplot(t(as.matrix(stacked.plot)), horiz = TRUE,  las=2, col = colors, xlab = "Methylation/Mb", ylab = "Chromosome", axes = TRUE)
#legend("topright", legend = c("Hypermethylation", "Hypomethylation"),fill = sequential[2:1], title = "Direction")
legend("topright", legend = c("Hypermethylation","Hypomethylation"), fill = colors, title = "Direction")
abline(0,0)
dev.off()

write.csv(Chromosome_Size, "Stacked_CpG_data.txt")
save.image("Stacked_CpG_freq.RData")
