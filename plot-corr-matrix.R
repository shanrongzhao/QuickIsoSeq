################################################################################
#
# Plot a correlation matrix
#
# Usage
#	Rscript plot-corr-matrix.R matrix.txt output_basename Figure-size[optional]
#
# Output
#	${output_basename}.corr.png 
#
# Shanrong Zhao
#
# November 6, 2018
#
#################################################################################

library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
arglen=length(args)
if ( arglen < 2) {
	print ("Rscript plot-corr-matrix.R tpm-matrix.txt $output_basename Figure-size[optional, default 9x7]")
	print ("Figure_Size:  $widthx$height (e.g  9x7). The unit is inch.")
	print ("For example: ")
	print ("	Rscript plot-expr-count.R salmon-tx-tpm.txt expr-count 8x7")
	q(save = "no")
}

# input data
corr = as.matrix( read.table(args[1], header=TRUE, sep="\t", row.names=1))
no <- nrow(corr)

# smart figure size
w <- floor (no/8) + 1
if ( w < 4 ) {
	w <- 4
}
h <- w-1
size <- paste0(w,"x",h)
#size <- "9x7"
if (arglen >2 ) {
	size <- tolower(args[3])
}
wh <- as.numeric(unlist(strsplit(size, "x")))

#
# plot
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

outfile <- paste(args[2],"9x7.png", sep=".")

# disable x-axis labelling or not
if ( no/8 > 9) {
	g.pair2 <- g.pair + scale_x_discrete(breaks=NULL) + scale_y_discrete(breaks=NULL)
	ggsave(plot=g.pair2, filename=outfile, width=9, height=7)
} else {
	ggsave(plot=g.pair, filename=outfile, width=9, height=7)
}

outfile <- paste(args[2],"png", sep=".")
ggsave(plot=g.pair, filename=outfile, width=wh[1], height=wh[2], limitsize=FALSE)
