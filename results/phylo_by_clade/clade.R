pdf("K27_byclades.pdf");
cH3K27me3_bordergenes <- read.table("H3K27me3_bordergenes_clade.dat",header=F,sep="\t");
pct <- round(cH3K27me3_bordergenes$V2/sum(cH3K27me3_bordergenes$V2)*100);
pct <- round(cH3K27me3_bordergenes$V2/sum(cH3K27me3_bordergenes$V2)*100);
lbls <- paste(cH3K27me3_bordergenes$V1, pct);
lbls <- paste(lbls,"%",sep="");
par(mar=c(10,0,2,2))
pie(cH3K27me3_bordergenes$V2,labels=lbls,col=rainbow(length(lbls)+1),radius=0.8,main="K27 H3K27me3_bordergenes Conservation");
par(mar=c(10,4,4,2))
barplot(cH3K27me3_bordergenes$V2,names.arg=lbls,las=2,ylab="count",main="K27 H3K27me3_bordergenes % (Conservation Counts)");
cH3K27me3_genes <- read.table("H3K27me3_genes_clade.dat",header=F,sep="\t");
pct <- round(cH3K27me3_genes$V2/sum(cH3K27me3_genes$V2)*100);
pct <- round(cH3K27me3_genes$V2/sum(cH3K27me3_genes$V2)*100);
lbls <- paste(cH3K27me3_genes$V1, pct);
lbls <- paste(lbls,"%",sep="");
par(mar=c(10,0,2,2))
pie(cH3K27me3_genes$V2,labels=lbls,col=rainbow(length(lbls)+1),radius=0.8,main="K27 H3K27me3_genes Conservation");
par(mar=c(10,4,4,2))
barplot(cH3K27me3_genes$V2,names.arg=lbls,las=2,ylab="count",main="K27 H3K27me3_genes % (Conservation Counts)");
cnonH3K27me3_genes <- read.table("nonH3K27me3_genes_clade.dat",header=F,sep="\t");
pct <- round(cnonH3K27me3_genes$V2/sum(cnonH3K27me3_genes$V2)*100);
pct <- round(cnonH3K27me3_genes$V2/sum(cnonH3K27me3_genes$V2)*100);
lbls <- paste(cnonH3K27me3_genes$V1, pct);
lbls <- paste(lbls,"%",sep="");
par(mar=c(10,0,2,2))
pie(cnonH3K27me3_genes$V2,labels=lbls,col=rainbow(length(lbls)+1),radius=0.8,main="K27 nonH3K27me3_genes Conservation");
par(mar=c(10,4,4,2))
barplot(cnonH3K27me3_genes$V2,names.arg=lbls,las=2,ylab="count",main="K27 nonH3K27me3_genes % (Conservation Counts)");
