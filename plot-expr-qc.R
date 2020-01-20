################################################################################
#
# Functions:
#		1. Plot the corrlation among all samples
#		2. generate the MADScore plot
#		3. boxplots of gene expressions
#		4. top/highly expressed genes
#
# Usage
#	Rscript plot-expr-qc.R tpm-expression-file output_basename Figure-size[optional]
#
# For example:
#	Rscript plot-expr-qc.R salmon-tx-tpm.txt expr-gene 8x7
#
#  Output
#			expr-gene-corr.txt
#			expr-gene-corr.*.png
#			expr-gene-MADScore.txt
#			expr-gene-MADScore.*.png
#			expr-gene-boxplot.*.png
#			expr-gene-top30-across-samples.txt
#			expr-gene-top30-across-samples.png
#			expr-gene-top10x.details.txt
#			expr-gene-top10x.txt
#			expr-gene-top10x.*.png
#
# Shanrong Zhao
#
# December 6, 2019
#
#################################################################################

library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

len=length(args)
if ( len < 2) {
	print ("Rscript plot-expr-count.R tpm-expression-file output_basename Figure-size[optional]")
	print ("Figure_Size:  $widthx$height (e.g  8x7). The unit is inch.")
	print ("For example: ")
	print ("	Rscript plot-expr-count.R salmon-tx-tpm.txt expr-gene-corr 8x7")
	q(save = "no")
}


if ( ! file.exists(args[1])) {
	print ( paste("TPM Expression File", args[1], "does not exist!", sep=" "))
	q(save = "no")
}

# smart figure size
size <- "8x7"
if (len >2 ) {
	size <- tolower(args[3])
}
wh <- as.numeric(unlist(strsplit(size, "x")))


# input data
data = read.table(args[1], header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=F)
expr <- as.matrix(data[,-c(1:6)])
samples <- colnames(expr)
no <- length( samples )

#get genes expression
expr.gene <- aggregate(expr, by=list(gene_ID=data$gene_ID, gene_name=data$gene_name), FUN=sum)
rownames(expr.gene) <- expr.gene[,1]
expr <- as.matrix(expr.gene[,-c(1,2)])
expr.annot <- expr.gene[, 1:2]


################################################################
#
# Top 30 genes, rank by the average of gene expressions
#
################################################################

n=30
total_exprs <- sum(expr)
total_feature_exprs <- rowSums(expr)
oo <- order(total_feature_exprs, decreasing = TRUE)

top30_pctage <- 100 * sum(total_feature_exprs[oo[1:n]]) / total_exprs

top30 <- (100 * t(expr[oo[1:n],]) / colSums(expr))
top30.genes <- expr.annot[colnames(top30),"gene_name"]
colnames(top30) <- top30.genes

# 
top30.out <- cbind( data.frame(Gene=top30.genes, round(t(top30),digits=2)))
write.table( top30.out, paste0(args[2],"-top30-across-samples.txt"), sep="\t", row.names = F, quote = F)

top30.mean <- colMeans(top30)
top30.mean <- data.frame(value=top30.mean,
	Gene= factor(names(top30.mean),levels= rev(top30.genes)) )

## Melt dataframe so it is conducive to ggplot
top30_long <- reshape2::melt(top30)
top30_long$Var2 <- factor(top30_long$Var2, levels= rev(top30.genes))

g.30 <- ggplot(top30_long, aes_string(y = "Var2", x = "value") ) +
	geom_point(alpha = 0.6, shape = 124) +
	ggtitle(paste0("Top ", n, " genes account for ",
				   format(top30_pctage, digits = 3), "% of total transcripts")) +
	ylab("Gene") + xlab(paste0("% of total TPM")) +
	theme( title = element_text(size=11))
	
g.30 <- g.30 + geom_point(
	aes_string(x = "value", y = "Gene"),
	data = top30.mean, colour = "red", shape = 21 )

outfile <- paste0(args[2],"-top30-across-samples.7x5.png")
ggsave(plot=g.30, filename=outfile, width=7, height=5)



################################################################
#
# Top 1 or 10 genes, rank by gene expression in individual samples
#
################################################################

total_exprs <- colSums(expr)
samples <- colnames(expr)
genes <- rownames(expr)

top10.tpm <- apply(expr, 2, function(x) { 
	oo <- order(x, decreasing = TRUE)
	x[oo[1:10]]
	} )
top10.pct <- t(top10.tpm/total_exprs * 100)
top10.gene <- t(apply(expr, 2, function(x) { 
	oo <- order(x, decreasing = TRUE)
	expr.gene[genes[oo[1:10]],"gene_name"]
	} ))

top10 <- cbind(data.frame(Sample=samples),data.frame(top10.gene),data.frame(round(top10.pct,digits=1)))
colnames(top10) <-c("Sample",paste0("Gene_",1:10),paste0("Pct_",1:10))
write.table( top10, paste0(args[2],"-top10x.details.txt"), sep="\t", row.names = F, quote = F)

mt.idx <- grep("^MT-|^Mt-", expr.annot$gene_name)
Pct_MT <- 100*colSums(expr[mt.idx,])/total_exprs

top10.summary <- cbind(data.frame(Sample=samples, Pct_MT=round(Pct_MT,digits=1)), 
	Pct_Top1=round(top10.pct[,1],digits=1), 
	Pct_Top10=round(rowSums(top10.pct), digits=1))
write.table( top10.summary, paste0(args[2],"-top10x.txt"), sep="\t", row.names = F, quote = F)

#
#plot top1 and top10
#
colnames(top10.summary) <- c("Sample","Mitochondria", "Top_1","Top_10")
data <- melt(top10.summary, id.vars=colnames(top10.summary)[1],
	measure.vars=colnames(top10.summary)[-1],
    variable.name="Top",
    value.name="Percent"
    )
data$Sample <- factor(data$Sample,  levels=as.character(unique(data$Sample)))

g.top10<- ggplot(data,aes(x=Sample, y=Percent))+
	geom_bar(stat="identity", width=0.7)+
	xlab("") + ylab("Percentage") + 
	theme(
		axis.text=element_text(size=9, face="bold"), 
		axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5),
		axis.title=element_text(size=9,face="bold"),
		strip.text=element_text(size=10, face="bold")) + 
	facet_grid(Top ~ .)


if (no/8 > 9) {
	ggsave(plot=g.top10+ scale_x_discrete(breaks=NULL), filename=paste0(args[2],"-top10x.9x7.png"), width=9, height=7)
} else {
	ggsave(plot=g.top10, filename=paste0(args[2],"-top10x.9x7.png"), width=9, height=7)
}
				
ggsave(plot=g.top10, filename=paste0(args[2],"-top10x.png"), width=wh[1], height=7,limitsize=FALSE)


#
#plot top1 and top10 without Mitochondria
#
top10.summary <- cbind(data.frame(Sample=samples), 
	Pct_Top1=round(top10.pct[,1],digits=1), 
	Pct_Top10=round(rowSums(top10.pct), digits=1))

colnames(top10.summary) <- c("Sample", "Top_1","Top_10")
data <- melt(top10.summary, id.vars=colnames(top10.summary)[1],
	measure.vars=colnames(top10.summary)[-1],
    variable.name="Top",
    value.name="Percent"
    )
data$Sample <- factor(data$Sample,  levels=as.character(unique(data$Sample)))

g2.top10<- ggplot(data,aes(x=Sample, y=Percent))+
	geom_bar(stat="identity", width=0.7)+
	xlab("") + ylab("Percentage") + 
	theme(
		axis.text=element_text(size=9, face="bold"), 
		axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5),
		axis.title=element_text(size=9,face="bold"),
		strip.text=element_text(size=10, face="bold")) + 
	facet_grid(Top ~ .)

if (no/8 > 9) {
	ggsave(plot=g2.top10+ scale_x_discrete(breaks=NULL), filename=paste0(args[2],"-top10x-noMT.9x7.png"), width=9, height=7)
} else {
	ggsave(plot=g2.top10, filename=paste0(args[2],"-top10x-noMT.9x7.png"), width=9, height=7)
}
				
ggsave(plot=g2.top10, filename=paste0(args[2],"-top10x-noMT.png"), width=wh[1], height=7,limitsize=FALSE)


################################################################
#
# correlation
#
################################################################

#filter out genes if 80% of them have tpm < 0.2
#in another word, a gene is excluded if it does not express in more than 80% samples
keep <- rowSums(expr > 0.2) > ncol(expr)*0.2
expr <- expr[keep,]
expr.log2 <- log2(expr +1 )
expr.annot <- expr.annot[keep,]

rm(list=c("expr.gene","data"))
gc()

corr <- cor(expr.log2)

corr.out <- cbind( data.frame(Sample=rownames(corr)), round(corr,digits=3))
write.table( corr.out, paste0(args[2],"-corr.txt"), sep="\t", row.names = F, quote = F)

n <- nrow(corr)
corr.data = melt(corr)

p <- ggplot(corr.data, aes(x=Var1, y=Var2, fill= value))

g.pair <- p + geom_tile() + 
	scale_fill_gradient2("Correlation",low="blue", high="red", 
		midpoint=median(corr.data$value)) + 
	theme_bw() + 
	labs( x="", y="") +
	scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
	theme( axis.ticks = element_blank(),
			axis.text.y = element_text(size=8),
			axis.text.x = element_text(size=8, hjust=1, vjust=0.5, angle = 90) )

outfile <- paste0(args[2],"-corr.9x7.png")
# disable x-axis and y-axis labelling or not
if (n/8 > 9) {
	ggsave(plot=g.pair+scale_x_discrete(breaks=NULL) + scale_y_discrete(breaks=NULL), filename=outfile, width=9, height=7)
} else {
	ggsave(plot=g.pair, filename=outfile, width=9, height=7)
}

outfile <- paste0(args[2],"-corr.png")
ggsave(plot=g.pair, filename=outfile, width=wh[1], height=wh[2], limitsize=FALSE)


################################################################
#
# MADscore
#
################################################################
MADscore <- function (M) {
	n <- nrow(M)
	diag(M) <- 0
	
	#the difference of averages
	sum.all <- sum(M)
	sum.each <- apply(M,1,sum)
	
	avg.each <- sum.each /(n-1)
	avg.rest <- (sum.all - 2 * sum.each)/((n-1)*(n-2))
	
	avg.delta <- avg.each - avg.rest

	#convert to z-score
	med <- median( avg.delta)
	MAD <- median (abs(avg.delta-med))
	MADScore <- (avg.delta-med)/ (MAD * 1.4826) 
	fail <- rep("N", n)
	fail[ MADScore < -5] <- "Y"
	
	data.frame(Avg_Corr=avg.each, Corr_diff=avg.delta, MADScore=MADScore, isOutlier=fail)
}

score=MADscore(corr)
score[,1:3] <- round(score[,1:3], digits=2)
qc <- cbind( data.frame(Sample=rownames(score)), score)

outfile <- paste0(args[2],"-MADScore.txt")
write.table( qc, outfile, sep="\t", row.names = F, quote = F)

#
# plot MADScore
#
qc$Sample <- factor(qc$Sample,  levels=as.character(unique(qc$Sample)))

g<- ggplot(qc, aes(x = Sample, y = MADScore)) + 
	geom_point() + xlab("") + ylab("MADScore") +  
	scale_color_discrete(guide=FALSE) +
	geom_point(data=qc[qc$MADScore < -3,], aes(x = Sample, y = MADScore), color="pink") +
	geom_point(data=qc[qc$MADScore < -5,], aes(x = Sample, y = MADScore), color="red") +
	geom_hline( yintercept = -5, lty="dashed") + 
	theme(
		axis.text=element_text(size=9, face="bold"), 
		axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
		axis.title=element_text(size=10,face="bold"))

ggsave(plot=g, filename=paste0(args[2],"-MADScore.9x7.png"), width=9, height=7)
if (no/8 > 9) {
	ggsave(plot=g+scale_x_discrete(breaks=NULL), filename=paste0(args[2],"-MADScore.9x7.png"), width=9, height=7)
} else {
	ggsave(plot=g, filename=paste0(args[2],"-MADScore.9x7.png"), width=9, height=7)
}
ggsave(plot=g, filename=paste0(args[2],"-MADScore.png"), width=wh[1], height=7,limitsize=FALSE)



################################################################
#
# Boxplot of gene expression
#
################################################################

png(filename=paste0(args[2],"-boxplot.9x7.png"), width = 9, height = 7, units = "in", res=300)
boxplot(expr.log2, ylab="log2TPM", las=2, cex.axis=0.7, ces.lab=0.9)
dev.off()

png(filename=paste0(args[2],"-boxplot.png"), width = wh[1], height = 7, units = "in", res=300)
boxplot(expr.log2, ylab="log2TPM", las=2, cex.axis=0.7, ces.lab=0.9)
dev.off()



################################################################
#
# PCA, output the top 25 components
#
################################################################

pca  <- prcomp(t(expr.log2))
pca.x <- pca$x

pca.no <- min(10, ncol(pca.x))
pca.out <- cbind(data.frame(Sample=rownames(pca.x)), round(pca.x[,1:pca.no],digits=2))

#PCA coordinate
outfile <- paste0(args[2],"-pca.x.txt")
write.table( pca.out, outfile, sep="\t", row.names = F, quote = F)

eig <- ( pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
pcs <- paste0("PC",1:length(pca$sdev))
eig.summary <- data.frame(PC=pcs, eigvalue = eig, variance = variance, cumvariance = cumvar)
eig.summary[,2:4] <- round(eig.summary[,2:4], digits=2)

#PCA variance
outfile <- paste0(args[2],"-pca.var.txt")
write.table( eig.summary[1:pca.no,], outfile, sep="\t", row.names = F, quote = F)

#contrib <- function(pca.x, comp.sdev, n.ind){
#  100*(1/n.ind)*pca.x^2/comp.sdev^2
#}
#ind.contrib <- t(apply(pca.x,1, contrib, pca$sdev, nrow(pca.x)))
#head(ind.contrib[, 1:4])

#PCA plot
xlabs <- paste0("PC1 (", round(eig.summary[1,"variance"],digits=2), "%)")
ylabs <- paste0("PC2 (", round(eig.summary[2,"variance"],digits=2), "%)")

png(filename=paste0(args[2],"-pca.png"), width = 8, height = 7, units = "in", res=300)
plot(pca.x[,1], pca.x[,2], pch = 19,  xlab=xlabs, ylab=ylabs)
abline(h=0, v=0, lty = 2)
if (nrow(pca.x) < 50) {
	text(pca.x[,1], pca.x[,2], labels=rownames(pca.x), cex=0.5, pos = 3)
}
dev.off()

