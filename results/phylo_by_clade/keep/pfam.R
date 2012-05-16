pdf("K27_bypfam.pdf");
pH3K27me3_bordergenes <- read.table("H3K27me3_bordergenes_pfam.dat",header=F,sep="\t");
pct <- round(pH3K27me3_bordergenes$V2/sum(pH3K27me3_bordergenes$V2)*100);
pct <- round(pH3K27me3_bordergenes$V2/sum(pH3K27me3_bordergenes$V2)*100);
lbls <- paste(pH3K27me3_bordergenes$V1, pct);
lbls <- paste(lbls,"%",sep="");
par(mar=c(10,4,4,2))
barplot(pH3K27me3_bordergenes$V2,names.arg=lbls,las=2,ylab="count",main="K27 H3K27me3_bordergenes% (Pfam counts)");
par(mar=c(10,0,2,2))
pie(pH3K27me3_bordergenes$V2,labels=lbls,col=rainbow(length(lbls)),main="K27 H3K27me3_bordergenes% Pfam");
pH3K27me3_genes <- read.table("H3K27me3_genes_pfam.dat",header=F,sep="\t");
pct <- round(pH3K27me3_genes$V2/sum(pH3K27me3_genes$V2)*100);
pct <- round(pH3K27me3_genes$V2/sum(pH3K27me3_genes$V2)*100);
lbls <- paste(pH3K27me3_genes$V1, pct);
lbls <- paste(lbls,"%",sep="");
par(mar=c(10,4,4,2))
barplot(pH3K27me3_genes$V2,names.arg=lbls,las=2,ylab="count",main="K27 H3K27me3_genes% (Pfam counts)");
par(mar=c(10,0,2,2))
pie(pH3K27me3_genes$V2,labels=lbls,col=rainbow(length(lbls)),main="K27 H3K27me3_genes% Pfam");
pnonH3K27me3_genes <- read.table("nonH3K27me3_genes_pfam.dat",header=F,sep="\t");
pct <- round(pnonH3K27me3_genes$V2/sum(pnonH3K27me3_genes$V2)*100);
pct <- round(pnonH3K27me3_genes$V2/sum(pnonH3K27me3_genes$V2)*100);
lbls <- paste(pnonH3K27me3_genes$V1, pct);
lbls <- paste(lbls,"%",sep="");
par(mar=c(10,4,4,2))
barplot(pnonH3K27me3_genes$V2,names.arg=lbls,las=2,ylab="count",main="K27 nonH3K27me3_genes% (Pfam counts)");
par(mar=c(10,0,2,2))
pie(pnonH3K27me3_genes$V2,labels=lbls,col=rainbow(length(lbls)),main="K27 nonH3K27me3_genes% Pfam");
