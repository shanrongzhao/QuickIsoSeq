################################################################################
#
# Plot summaries for RNA-Seq read mapping, counting or read distribution
#
# Usage
#	Rscript plot-rnaseq-metrics.R Datafile output_basename Figure_Size[optional, default 8X6]
#
# Shanrong Zhao
#
# November 14, 2018
#
#################################################################################

library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

arglen=length(args)
if ( arglen < 2) {
	print ("Rscript plot-rnaseq-metrics.R Datafile Output_basename Figure_size[optional] Feature_Name[optional]")
	print ("Figure_Size:  $widthx$height (e.g  9x7). The unit is inch.")
	print ("Feature name:  Mapping, Counting or Reads_split and etc.")
	q(save = "no")
}

# input data
data = read.table(args[1], header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=F)
data = data[, -2]
no=nrow(data)

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

feature_name <- "Feature"
if (arglen >3 ) {
	feature_name <- args[4]
}

# 
# 
datam <- melt(data, id.vars=colnames(data)[1],
	measure.vars=colnames(data)[-1],
    variable.name="Feature",
    value.name="Percent"
    )
	
datam$Sample <- factor(datam$Sample,  levels=as.character(unique(datam$Sample)))

g<- ggplot(datam,aes(x=Sample, y=Percent))+
	geom_bar(stat="identity",aes(fill=Feature))+
	xlab("") + ylab("Percentage of Reads") + labs( fill = feature_name) + 
	theme(
		axis.text=element_text(size=9, face="bold"), 
		axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5),
		axis.title=element_text(size=9,face="bold")) + 
	guides(fill = guide_legend(
				title.theme = element_text(size=10,angle=0,,face="bold"), 
				label.theme = element_text(size=9,angle=0,face="bold"), 
				keywidth = 0.8, keyheight = 0.8
				)) 

#ggsave(plot=g, filename="Test.9x7.png"), width=9, height=7)
ggsave(plot=g, filename=paste0(args[2],".9x7.png"), width=9, height=7)

if (no/8 > 9) {
	ggsave(plot=g+scale_x_discrete(breaks=NULL), filename=paste0(args[2],".9x7.png"), width=9, height=7)
} else {
	ggsave(plot=g, filename=paste0(args[2],".9x7.png"), width=9, height=7)
}
				
ggsave(plot=g, filename=paste0(args[2],".png"), width=wh[1], height=wh[2],limitsize=FALSE)