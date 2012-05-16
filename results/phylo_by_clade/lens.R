pdf("K27_peplens.pdf")
lH3K27me3_bordergenes <- read.table("H3K27me3_bordergenes_lengths.dat",header=T,sep="\t");
summary(lH3K27me3_bordergenes)
lH3K27me3_genes <- read.table("H3K27me3_genes_lengths.dat",header=T,sep="\t");
summary(lH3K27me3_genes)
lnonH3K27me3_genes <- read.table("nonH3K27me3_genes_lengths.dat",header=T,sep="\t");
summary(lnonH3K27me3_genes)
boxplot(subset(lH3K27me3_bordergenes$LEN,lH3K27me3_bordergenes$LEN < 4000),subset(lH3K27me3_genes$LEN,lH3K27me3_genes$LEN < 4000),subset(lnonH3K27me3_genes$LEN,lnonH3K27me3_genes$LEN < 4000),xlab="K27 Status", main="Protein size for K27 status",ylab="Protein length",range=2,names=c("lH3K27me3_bordergenes","lH3K27me3_genes","lnonH3K27me3_genes"),col=rainbow(11, start=0, end=1));