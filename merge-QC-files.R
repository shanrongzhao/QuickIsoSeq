################################################################################
#
# Usage
#	Rscript merge-QC-files.R file1 file2 ......
#
#
# The output:  
#			RNASeq-merged-metrics.txt
#
# Shanrong Zhao
#
# November 8, 2018
#
################################################################################

args <- commandArgs(trailingOnly = TRUE)

len=length(args)
if ( len < 1) {
	print ("Rscript merge-QC-files.R file1 file2 ......")
	print ("The output merged file: RNASeq-merged-metrics.txt")
	q(save = "no")
}


#read all sample id
merged <- read.table( args[1], sep="\t",header=TRUE, stringsAsFactors=F, check.names = F)

if ( len > 1 ) {
	for (i in 2:len ) {
		f2 <- read.table(args[i], sep="\t", header=TRUE, row.names=1, stringsAsFactors=F, check.names = F)
		merged <- cbind(merged,f2)
	}
}

colnames(merged)[1] <- "Sample"
write.table(merged, "RNASeq-merged-metrics.txt", sep="\t",row.names = F, quote = F)
