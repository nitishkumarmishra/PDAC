############## Positive correlations ##############
logp<-(-log10(diffmeth.merge.na.01.hyper$adj.P.Val))
disttss<-as.numeric(diffmeth.merge.na.01.hyper$dist)
pick<-disttss<=1500 & disttss>=-1500
disttssf<-disttss[pick]
logpf<-logp[pick]
############## Negative correlations###############
nlogp<-(-log10(diffmeth.merge.na.01.hypo$adj.P.Val))
ndisttss<-as.numeric(diffmeth.merge.na.01.hypo$dist)
pick<-ndisttss<=1500 & ndisttss>=-1500
ndisttssf<-ndisttss[pick] 
nlogpf<-nlogp[pick]

pdf("Plot_50Kb.pdf",width=10,height=8)
plot(ndisttssf,nlogpf,col="blue4",pch=20,cex=.5,
     ylim=c(min(c(logpf,nlogpf)),max(c(logpf,nlogpf))),
     xlab="Distance between dm-CpGs and TSS (bp)",
     ylab="-log10(BH adjusted P-value)",
     xaxt="n"
)

axis(1,at=seq(-1.5e3,1.5e3,length.out=7),labels=c("-1500","-1000", "-500","0", "500","1000","1500"))
legend("topleft", c("Hypermethylation","Hypomethylation"), text.col=c("red4","blue4"),bty="n", col=c("red4","blue4"), lwd=2, lty=c(0,0),pch=c(20,20))
points(disttssf,logpf,col="red4",pch=20,cex=.5)
dev.off()
