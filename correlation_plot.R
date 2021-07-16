library(data.table)
library(cowplot)
library(ggplot2)
library(ggpubr)

exp_meth_protein_coding <- read.csv("exp_meth_protein_coding_11_6_2018.txt", header = TRUE, check.names = FALSE)

##columns : 6 to 129 exp_meth_protein_coding[sample_gene,6:129] data 124 Cancer)
##columms :130:253 exp_meth_protein_coding[sample_gene,130:253] data 124 cancer)


# exp <- exp_meth_protein_coding[,6:129]
# exp<-exp[ , order(names(exp))]
# exp <- log(exp)
# 
# meth <- exp_meth_protein_coding[,130:253]
# meth<-meth[ , order(names(meth))]


df <- data.frame(matrix(ncol = 3, nrow = 584))  ## 137*3  # 146 Tumor *4 Probes = 584
colnames(df) <- c("exp","meth","Probe")
datalist =list()
count = 1

for (gene_name in c("cg26696870","cg05292954","cg04956949","cg01181415")){
  
  Probe = which(gene_name == exp_meth_protein_coding$ProbeID)
  

  
  datalist[[count]] <- exp_meth_protein_coding[Probe,]
  
  count = count + 1
  
  
} 
  
dataset <- rbind(datalist[[1]],datalist[[2]],datalist[[3]],datalist[[4]])

df$exp <- c(dataset[1,6:151],dataset[2,6:151],dataset[3,6:151],dataset[4,6:151])
df$meth <- c(dataset[1,152:297],dataset[2,152:297],dataset[3,152:297],dataset[4,152:297])
df$Probe <- as.factor(c(rep("cg26696870",146),rep("cg20127859",146),rep("cg04956949",146),rep("cg01181415",146)))


df$exp <- unlist(df$exp)
df$meth <- unlist(df$meth)

 # Main plot
  pmain <- ggplot(df, aes(x = exp, y = meth, color = Probe)) +
    geom_point()+ geom_smooth(method='lm',se=FALSE) +
    ggpubr::color_palette("aaas") +
    stat_cor(aes(x = exp, y = meth, color = Probe), method = "spearman",label.x.npc = 0.5, 
             label.y.npc = 1, hjust = -0.5) 
     
  
  # Marginal densities along x axis
  xdens <- axis_canvas(pmain, axis = "x")+
    geom_density(data = df, aes(x = exp, fill = Probe),
                 alpha = 0.7, size = 0.2)+
    ggpubr::fill_palette("aaas")
  # Marginal densities along y axis
  # Need to set coord_flip = TRUE, if you plan to use coord_flip()
  ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
    geom_density(data = df, aes(x = meth, fill = Probe),
                 alpha = 0.7, size = 0.2)+
    coord_flip()+
    ggpubr::fill_palette("aaas")
  p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.4, "null"), position = "top")
  p2<- insert_yaxis_grob(p1, ydens, grid::unit(.4, "null"), position = "right")
  ggdraw(p2)
  

  ggsave(filename = "Correlation_plot_11_13_2018.pdf", width = 10, height = 8, dpi = 800)
  
  
  
  