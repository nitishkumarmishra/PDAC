## This is the R code to generate circos plot of DNA methylation and genomic density
library(circlize)
load(paste0(system.file(package = "circlize"), "/extdata/DMR.RData"))
par(mar =c(1, 1, 1, 1))
#circos.initializeWithIdeogram(plotType =c("axis", "labels"))
circos.initializeWithIdeogram()
bed_list =list(DMR_hyper, DMR_hypo)
circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col =c("#FF000080", "#0000FF80"))
circos.genomicDensity(DMR_hyper, col =c("#FF000080"), track.height = 0.1)
circos.genomicDensity(DMR_hypo, col =c("#0000FF80"), track.height = 0.1)