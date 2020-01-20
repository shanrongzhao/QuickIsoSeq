################################################################################
#
# Plot summaries for RNA-Seq read mapping, counting or read distribution
#
# Usage
#	Rscript plot-lib-size.R star-mapping-summary.txt output_basename Figure_Size[optional, default 8x6]
#
# Shanrong Zhao
#
# November 16, 2018
#
#################################################################################

library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

arglen=length(args)
if ( arglen < 2) {
	print ("Rscript plot-lib-size.R star-mapping-summary.txt output_basename Figure_size[optional]")
	print ("Figure_Size:  $widthx$height (e.g  9x7). The unit is inch.")
	q(save = "no")
}

# input data
data = read.table(args[1], header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=F)
data.reads <- data[, 1:2]
colnames(data.reads) <- c("Sample","Total_reads")
data.reads$Sample <- factor(data.reads$Sample,  levels=data.reads$Sample)
no = nrow(data.reads)

# smart figure size
w <- floor (no/8) + 1
if ( w < 4 ) {
	w <- 4
}
size <- paste0(w,"x6")
if (arglen >2 ) {
	size <- tolower(args[3])
}
wh <- as.numeric(unlist(strsplit(size, "x")))

#
# plot
g.reads <- ggplot(data.reads, aes(x=Sample, y=Total_reads) )+
	geom_bar(stat="identity", width=0.6, fill="gray")+ 
	xlab(NULL) + ylab("Total number of sequenced reads")  +
	theme(axis.text=element_text(size=8,colour="black"), 
		axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
		axis.title=element_text(size=10,face="bold")) 
#g.reads

if (no/8 > 9 ) {
	ggsave(plot=g.reads + scale_x_discrete(breaks=NULL), filename=paste0(args[2],".9x7.png"), width=9, height=7)
} else {
	ggsave(plot=g.reads, filename=paste0(args[2],".9x7.png"), width=9, height=7)
}
				
ggsave(plot=g.reads, filename=paste0(args[2],".png"), width=wh[1], height=wh[2],limitsize=FALSE)
