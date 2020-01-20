################################################################################
#
# Usage
#	Rscript extractBy.R annot_file key_file new_annot_file
#
#
# Shanrong Zhao
#
# November 4, 2018
#
################################################################################

args <- commandArgs(trailingOnly = TRUE)

len=length(args)
if ( len < 3) {
	print ("Rscript extractBy.R annot_file key_file new_annot_file")
	print ("key_file Ccontains a single column without header")
	q(save = "no")
}

#read annotation
annot <- read.table( args[1], sep="\t",header=TRUE, stringsAsFactors=F, check.names = F)
rownames(annot) <- annot[,1]

#read keys
keys <- read.table( args[2], sep="\t",header=FALSE, stringsAsFactors=F, check.names = F)[,1]

#write extratced rows in the order of key
extracted <- annot[keys,]
write.table(extracted, args[3], sep="\t",row.names = F, quote = F)
