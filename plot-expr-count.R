################################################################################
#
# Plot the number of expressed genes and transcripts at given cut-offs
#
# Usage
#	Rscript plot-expr-count.R tpm-expression-file output_basename Figure-size [optional]
#
# For example:
#	Rscript plot-expr-count.R salmon-tx-tpm.txt expr-count 8x7
#
#  Output
#			expr-count-gene.png
#			expr-count-tx.png
#
# Shanrong Zhao
#
# November 6, 2018
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
	print ("	Rscript plot-expr-count.R salmon-tx-tpm.txt expr-count 8x7")
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

#
# transcripts
#
base_name <- paste0(args[2], "-tx")

dat <- data.frame(Sample=samples, TPM=rep(0,no), count=apply(expr>0, 2, sum))
for (cutoff in c(1, 10)) {
	temp <- data.frame(Sample=samples, TPM=rep(cutoff,no), count=apply(expr>cutoff,2, sum))
	dat <- rbind(dat, temp)
}

# export table for interactive plot
dat_wide0 <- dcast(dat, Sample ~ TPM, value.var="count")
samplesMatrix <- data.frame(Sample=samples)
dat_wide <- merge (samplesMatrix, dat_wide0, all.x=TRUE, sort=FALSE, by="Sample")
write.table(format(dat_wide, digits=4), file = paste0(base_name,".txt"), sep="\t", quote=FALSE, row.names=FALSE )

dat$Sample <- factor(dat$Sample,  levels=as.character(unique(dat$Sample)))
g<- ggplot(dat, aes(x = Sample, y = count, colour = factor(TPM))) + 
	geom_point() + xlab("") + ylab("Number of expressed transcripts") +  
	scale_color_discrete(name  ="TPM") + 
	guides(colour = guide_legend(override.aes = list(size=3), 
		title.theme = element_text(size=10,angle=0,face="bold"), 
		label.theme = element_text(size=9,angle=0,face="bold"))) +  
	theme(
		axis.text=element_text(size=9, face="bold"), 
		axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
		axis.title=element_text(size=10,face="bold"))

if (no/8 > 9) {
	ggsave(plot=g+ scale_x_discrete(breaks=NULL), filename=paste0(base_name,".9x7.png"), width=9, height=7)
} else {
	ggsave(plot=g, filename=paste0(base_name,".9x7.png"), width=9, height=7)
}
ggsave(plot=g, filename=paste0(base_name,".png"), width=wh[1], height=7,limitsize=FALSE)



#
# gene
#
base_name <- paste0(args[2], "-gene")

expr.gene <- aggregate(expr, by=list(gene_ID=data$gene_ID), FUN=sum)
expr <- as.matrix(expr.gene[,-1])

dat <- data.frame(Sample=samples, TPM=rep(0,no), count=apply(expr>0, 2, sum))
for (cutoff in c(1, 10)) {
	temp <- data.frame(Sample=samples, TPM=rep(cutoff,no), count=apply(expr>cutoff,2, sum))
	dat <- rbind(dat, temp)
}

# export table for interactive plot
dat_wide0 <- dcast(dat, Sample ~ TPM, value.var="count")
samplesMatrix <- data.frame(Sample=samples)
dat_wide <- merge (samplesMatrix, dat_wide0, all.x=TRUE, sort=FALSE, by="Sample")

write.table(format(dat_wide, digits=4), file = paste0(base_name,".txt"), sep="\t", quote=FALSE, row.names=FALSE )

# 
dat$Sample <- factor(dat$Sample,  levels=as.character(unique(dat$Sample)))
g<- ggplot(dat, aes(x = Sample, y = count, colour = factor(TPM))) + 
	geom_point() + xlab("") + ylab("Number of expressed genes") +  
	scale_color_discrete(name  ="TPM") + 
	guides(colour = guide_legend(override.aes = list(size=3), 
		title.theme = element_text(size=10,angle=0,face="bold"), 
		label.theme = element_text(size=9,angle=0,face="bold"))) +  
	theme(
		axis.text=element_text(size=9, face="bold"), 
		axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
		axis.title=element_text(size=10,face="bold"))

if (no/8 > 9) {
	ggsave(plot=g+ scale_x_discrete(breaks=NULL), filename=paste0(base_name,".9x7.png"), width=9, height=7)
} else {
	ggsave(plot=g, filename=paste0(base_name,".9x7.png"), width=9, height=7)
}
ggsave(plot=g, filename=paste0(base_name,".png"), width=wh[1], height=7,limitsize=FALSE)

