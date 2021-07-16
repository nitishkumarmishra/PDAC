library(ggpubr)
### R code for whisker plot for gene regions and genomic regions
load("F:/OneDrive - University of Nebraska Medical Center/Pancreatic Endocrine/TCGA Pancreatic Ductal/PAAD_Methylation_Oct-10.RData")
rm(list=setdiff(ls(), c("m1", "BMIQ.Meth")))
hm450 <- read.csv("../../Human450K-Extra/HuiShen/September 2018/hm450.hg38.manifest.tsv.gz", header = TRUE, sep = "\t")
rownames(hm450) <- hm450$probeID
hm450.subset <- subset(x = hm450, select = c("CpG_chrm", "MASK_general"))
hm450.subset <- hm450.subset[!as.character(hm450.subset$CpG_chrm) %in% c("chrX", "chrY", "chrM"),]
hm450.subset <- hm450.subset[hm450.subset$MASK_general==FALSE,]

BMIQ.mask <- BMIQ.Meth[rownames(BMIQ.Meth)%in%rownames(hm450.subset),]


hm450.orig <- read.csv("../../Human450K-Extra/HuiShen/September 2018/hm450.hg19.manifest.original.tsv.gz", header = TRUE, sep = "\t")
rownames(hm450.orig) <- hm450.orig$Name

hm450.orig.subset <- subset(x = hm450.orig,select = c("Relation_to_UCSC_CpG_Island", "UCSC_RefGene_Group", "CHR"))


annotation <- apply(hm450.orig.subset, 2, function(x) gsub("^$|^ $", NA, x))
tmp1 <- sapply(strsplit(annotation[,2], ";"), "[[", 1)

hm450.orig.subset$UCSC_RefGene_Group <- tmp1
BMIQ.mask <- merge(BMIQ.mask, hm450.orig.subset, by="row.names")
rownames(BMIQ.mask) <- BMIQ.mask$Row.names
BMIQ.mask$Row.names <- NULL

rm(list=setdiff(ls(), c("BMIQ.mask")))

BMIQ.mask.copy <- BMIQ.mask

test<-read.csv("DiffMeth_No_chrXYM_P.val_0.01_UCSC.txt",header = TRUE)
BMIQ.mask <- BMIQ.mask[rownames(BMIQ.mask) %in% test$Row.names,]

BMIQ.mask.TSS200 <- BMIQ.mask[grep("TSS200", BMIQ.mask$UCSC_RefGene_Group),]
BMIQ.mask.TSS1500 <- BMIQ.mask[grep("TSS1500", BMIQ.mask$UCSC_RefGene_Group),]
BMIQ.mask.5UTR <- BMIQ.mask[grep("5'UTR", BMIQ.mask$UCSC_RefGene_Group),]
BMIQ.mask.1stExon <- BMIQ.mask[grep("1stExon", BMIQ.mask$UCSC_RefGene_Group),]
BMIQ.mask.Body <- BMIQ.mask[grep("Body", BMIQ.mask$UCSC_RefGene_Group),]
BMIQ.mask.3UTR <- BMIQ.mask[grep("3'UTR", BMIQ.mask$UCSC_RefGene_Group),]

BMIQ.mask.Island <- BMIQ.mask[grep("Island", BMIQ.mask$Relation_to_UCSC_CpG_Island),]
BMIQ.mask.S_Shore <- BMIQ.mask[grep("S_Shore", BMIQ.mask$Relation_to_UCSC_CpG_Island),]
BMIQ.mask.N_Shore <- BMIQ.mask[grep("N_Shore", BMIQ.mask$Relation_to_UCSC_CpG_Island),]
BMIQ.mask.S_Shelf <- BMIQ.mask[grep("S_Shelf", BMIQ.mask$Relation_to_UCSC_CpG_Island),]
BMIQ.mask.N_Shelf <- BMIQ.mask[grep("N_Shelf", BMIQ.mask$Relation_to_UCSC_CpG_Island),]


#Tumor.BMIQ <- read.csv("Tumor_BMIQ.txt", header = TRUE)

#as.vector(t(Tumor.BMIQ))
################### This is for TSS200 ########################
################### Similar for other region ##################

for(Region_type in c("TSS200", "TSS1500","5UTR","1stExon","Body","3UTR",
                     "Island","S_Shore","N_Shore","S_Shelf","N_Shelf")){


#Tumor <- BMIQ.mask.TSS200[,grep("01A", substr(colnames(BMIQ.mask.TSS200), 14, 16))]
#Normal <- BMIQ.mask.TSS200[,grep("11A", substr(colnames(BMIQ.mask.TSS200), 14, 16))]

Tumor <- eval(parse(text = (paste0("BMIQ.mask.",Region_type))))[,grep("01A", substr(colnames(eval(parse(text = (paste0("BMIQ.mask.",Region_type))))), 14, 16))]
Normal <- eval(parse(text = (paste0("BMIQ.mask.",Region_type))))[,grep("11A", substr(colnames(eval(parse(text = (paste0("BMIQ.mask.",Region_type))))), 14, 16))]



tmp <- as.vector(t(Tumor))
tmp1 <- as.vector(t(Normal))

tt <- rep("Tumor", length(tmp))
tt1 <- rep("Normal", length(tmp1))
tumor <- cbind(tmp, tt)
normal <- cbind(tmp1, tt1)
df <- as.data.frame(rbind(tumor, normal))
df$tmp <- as.numeric(as.character(df$tmp))

#compare_means(tmp ~ tt,  data = df, ref.group = ".all.", method = "wilcox.test")
ggboxplot(df, x = "tt", y = "tmp", xlab = "Status",ylab = "Beta value", color = "tt", palette = "aaas",legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(df$tmp), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "wilcox.test", label.y = max(df$tmp)+0.2)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method  = "t.test", 
                     ref.group = ".all.") # Pairwise comparison against all


  ggsave(paste0(Region_type,"Boxplot.pdf"), width = 8, height = 8)
  

}
