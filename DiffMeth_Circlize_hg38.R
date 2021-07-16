library(gtrellis);library(circlize); library(ComplexHeatmap); library(naturalsort)

setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
bed <- read.csv("DiffMeth_No_chrXYM_P.val_0.01.txt", header = TRUE, row.names = 1)
bed <- subset(bed, select = c(CpG_chrm, CpG_beg, CpG_end, log2FC))
names(bed)[1] <- c("chr"); names(bed)[2] <- c("start"); names(bed)[3] <- c("end")
bed$direction <- ifelse(bed$log2FC >0, "Hyper", "Hypo")
bed$log2FC <- NULL
bed <- bed[naturalorder(bed$chr),]### Order dataframe with chromosome name

DMR_hyper = bed[bed$direction == "Hyper", ]
DMR_hypo = bed[bed$direction == "Hypo", ]
bed_list = list(DMR_hyper, DMR_hypo)

distance_CpG <- rainfallTransform(bed)
density_CpG <- genomicDensity(bed, window.size = 1000000)

density_CpG.ordered <- density_CpG[order(density_CpG$pct, decreasing = TRUE),]

Hyper <- bed[grep("Hyper", bed$direction),]
Hypo <- bed[grep("Hypo", bed$direction),]
density_Hyper <- genomicDensity(Hyper, window.size = 1000000)
density_Hypo <- genomicDensity(Hypo, window.size = 1000000)
density_Hyper.ordered <- density_Hyper[order(density_Hyper$pct, decreasing = TRUE),]
density_Hypo.ordered <- density_Hypo[order(density_Hypo$pct, decreasing = TRUE),]

pdf("Circlize Rainfall Plot logFC2.pdf", width = 8, height = 8)
circos.clear() # this line needed before circos.par
circos.par("start.degree" = 90) ## This line to make chr1 on top and center for starting plot
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22), species = "hg38")
circos.genomicRainfall(bed_list, pch = 16, cex = 0.5, track.height = 0.2, col = c("red4", "blue4"))
circos.genomicDensity(DMR_hyper, col = c("red4"), track.height = 0.2, window.size = 5000000)
circos.genomicDensity(DMR_hypo, col = c("blue4"), track.height = 0.2, window.size = 5000000)
legend("center", bty="n", cex=1, legend = c("Hypermethylated", "Hypomethylated"), col = c("red4","blue4"), text.col = c("red4","blue4"), pch = 16)
dev.off()

###########################################
save.image("DiffMeth_Circlize_hg38.RData")