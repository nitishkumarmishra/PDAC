## DMRcate analysis worked on R 4.3.3 (not on R 5.1)
library(DMRcate)
library(limma)
setwd("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal")
load("PAAD_Methylation_Oct-10.RData")
rm(list=setdiff(ls(), c("PAAD.Clinical", "BMIQ.Meth")))

myMs <- logit2(BMIQ.Meth)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)


cancerID <- grep("01A", substr(colnames(myMs.noSNPs), 14, 16))
normalID <- grep("11A", substr(colnames(myMs.noSNPs), 14, 16))

c1 <- length(cancerID)
c2 <- length(normalID)

Tumor.BMIQ <- myMs.noSNPs[,cancerID]
Normal.BMIQ <- myMs.noSNPs[, normalID]

groups <- as.factor(c(rep("Tumor",c1),rep("Normal",c2)))


design<-model.matrix(~0+groups)
colnames(design)=levels(groups)
contrast.matrix <- makeContrasts(Tumor-Normal, levels = design)


myannotation <- cpg.annotate("array",myMs.noSNPs, analysis.type = "differential", design = design, contrasts = TRUE, what="M", fdr = 0.01, cont.matrix = contrast.matrix, coef = "Tumor - Normal")
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, p.adjust.method = "BH", pcutoff = "fdr", betacutoff=0.2)
#### Annotate overlapping promoter regions (+/- 2000 bp from TSS)
results.ranges <- extractRanges(dmrcoutput, genome = "hg38")
groups <- c(Tumor="magenta", Normal="forestgreen")
type <- as.factor(c(rep("Tumor", c1), rep("Normal", c2)))
cols <- groups[as.character(type)]
samps <- c(1:153, 153+(1:9))


DMR.plot(ranges=results.ranges, dmr=1, CpGs=meth, phen.col=cols, genome="hg38", arraytype = "450K", what="Beta", samps=samps, showSampleNames = FALSE, col.axis="Black", col.title="Black",cex.sampleNames = 0.6, col.sampleNames="Black", 
         separator = 1.75, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE, cex = 0.6, fontcolor="Black", fontcolor.group="Black", fontsize.group=8, fontcolor.feature = 1, cex.feature = 0.6, min.width = 3, min.distance = 5, cex.title=0.5, rotation.title=360)

## Plot is not working. To many samples and to many CpGs
par(mfrow=c(1,1))


#DMR.plot(ranges=results.ranges, dmr=797, CpGs=BMIQ.Meth, phen.col=cols, genome="hg38", arraytype = "450K", what="Beta", samps=samps, showSampleNames = TRUE, cex.sampleNames = 0.8, separator = 1, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE)
DMR.plot(ranges=results.ranges, dmr=797, CpGs=BMIQ.Meth, phen.col=cols, genome="hg38", arraytype = "450K", what="Beta", samps=samps, showSampleNames = FALSE, col.axis="Black", col.title="Black",cex.sampleNames = 0.6, col.sampleNames="Black", 
         separator = 1.75, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE, cex = 0.6, fontcolor="Black", fontcolor.group="Black", fontsize.group=8, fontcolor.feature = 1, cex.feature = 0.6, min.width = 3, min.distance = 5, cex.title=0.5, rotation.title=360)


rm(c1, c2, numNAs, myNorm, Meth, tx.hg19, tx.hg38, tx.mm10, probe.features, probeInfoALL.lv)


df1 <- data.frame(seqnames=seqnames(results.ranges),
                  starts=start(results.ranges),
                  ends=end(results.ranges),
                  names=c(rep(".", length(results.ranges))),
                  strand=strand(results.ranges)
)
df <- mcols(results.ranges)
df2 <- cbind(df1, df)
write.csv(df2, file = "DMRcateDMR.txt")
write.csv(df2, file = "DMRcateDMR.csv")

save.image("DMRcate_analysis")

###################################
## If I use all tumor and sample DMR.plot will not work.
## This part to make plot. I tried upto 100 cancer and 9 normal, its workin. But it's quality is not good, so I used only 75 tumor samples.
colnames(BMIQ.Meth) <- substr(colnames(BMIQ.Meth), 1, 15)
random <- sample(1:153, 150, replace = TRUE)
random <- unique(random)
random <- random[1:75]
meth <- BMIQ.Meth[,c(random, normalID)]

type <- as.factor(c(rep("Tumor", 75), rep("Normal", 9)))
cols <- groups[as.character(type)]
samps <- c(1:75, 75+(1:9))

if (!dir.exists("DMRcate_Plot")) dir.create("DMRcate_Plot")
for (i in 1:100) {
  
  DMR.plot(ranges=results.ranges, dmr=i, CpGs=meth, phen.col=cols, genome="hg38", arraytype = "450K", what="Beta", samps=samps, showSampleNames = FALSE, col.axis="Black", col.title="Black",cex.sampleNames = 0.6, col.sampleNames="Black", 
           separator = 1.75, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE, cex = 0.6, fontcolor="Black", fontcolor.group="Black", fontsize.group=8, fontcolor.feature = 1, cex.feature = 0.6, min.width = 3, min.distance = 5, cex.title=0.5, rotation.title=360)
  
  file_name = paste0("DMRcate_Plot", "/DMR_", i, ".pdf")
  dev.print(pdf, file_name, width = 9, height = 9)
}
