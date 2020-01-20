################################################################################
#
# REORDER snp matrix by subject ID in sample.annotation.txt so that all samples 
# are in the right order
#
# Usage
#	Rscript reorder-snp-corr.R SNP.corr.txt sample.annotation.txt <sort_method>[optional]
#
#
# Shanrong Zhao
#
# November 16, 2018
#
#################################################################################

args <- commandArgs(trailingOnly = TRUE)
len=length(args)
if ( len < 2) {
	print ("Rscript reorder-snp-corr.R <correlation file> <annotation file> <sort_method<[optional]")
	print ("Annotation_file: sample_id:1st column; subject_id:2nd column")
	print ("sort_method: byorder: by the order of unique names; byname: by the subject_id name")
	q(save = "no")
}

corr = read.table(args[1], header=TRUE, sep="\t", check.names = F)
rownames(corr) <- corr[,1]

annot = read.table(args[2], header=TRUE, sep="\t", check.names = F)
rownames(annot) <- annot[,1]

sort_method <- "byorder"
if (len == 3) {
	sort_method <- tolower(args[3])
}
if ( sort_method == "byorder") { 
	annot[,2] = factor( annot[,2], levels=unique(annot[,2]) )
} 

annot <- annot[order(annot[,2]),]
samples <- rownames(annot)
corr <- corr[samples, c(colnames(corr)[1],samples)]

write.table(corr,args[1], sep="\t",row.names = F, quote = F)

