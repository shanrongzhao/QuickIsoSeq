################################################################################
#
# merge counts (and annotation) from salmon_aln 
#
# Usage
#	      Rscript get-salmon_aln-tpm.R allIDs.txt transcript_annotation_file
#
# The output:  
#			salmon_aln-tx-tpm.txt
#			salmon_aln-gene-tpm.txt
#			salmon_aln-tx-counts.txt
#			salmon_aln-gene-counts.txt
#
# Input
#
#Name	Length	EffectiveLength	TPM	NumReads
#ENST00000456328.2	1657	1493.62	0.0380841	2.98629
#ENST00000450305.2	632	468.69	0.026642	0.655543
#ENST00000488147.1	1351	1187.62	15.3098	954.542
#ENST00000619216.1	68	8.89449	0	0
#
# Chi Zhang, Shanrong Zhao
#
# November 4, 2018
#
################################################################################

args <- commandArgs(trailingOnly = TRUE)
len<-length(args)
if ( len < 2) {
	print ("Please input <id file> and <transcript annotation file> \n")
	q(save = "no")
}

manifest.file <- args[1]
manifest <- read.table(manifest.file, header=F, check.names = F)
samples <- manifest[,1]

annot_file <- args[2]
if ( ! file.exists(annot_file)) {
        print ( paste("Annotation file", annot_file, " does not exist!", sep=" "))
        q(save = "no")
}
annote <- read.table(annot_file, header=T, stringsAsFactors=F, check.names=F, sep="\t")


tpm <- list()
count <- list()
i <- 0
for (id in samples) {
	filename <- paste(id, "/salmon_aln/quant.sf", sep="")
	
	data <- read.table(filename, header=TRUE, stringsAsFactors=F, check.names = F)

	tpm[[id]] <- data[,4]
	count[[id]] <- data[,5]
	if ( i < 1 ) {
		transcript_ID <- data[,1]
	}
	i <- i+1
}


#
# tpm
#
tpm <- round(as.data.frame( do.call("cbind", tpm)), digits=3)
tpm <- cbind(data.frame(transcript_ID = transcript_ID), tpm)
tpm <- merge( annote, tpm, by = "transcript_ID" )

tpm <- tpm[with(tpm, order(gene_ID,transcript_ID)),]
write.table(tpm, "salmon_aln-tx-tpm.txt", sep="\t", row.names =F, quote = F)

tpm.gene <- aggregate(tpm[,-c(1:6)], 
	by=list(gene_ID=tpm$gene_ID, gene_name=tpm$gene_name), FUN=sum)
write.table(tpm.gene, "salmon_aln-gene-tpm.txt", sep="\t", row.names =F, quote = F)

#
# count
#
count <- round(as.data.frame( do.call("cbind", count)), digits=3)
count <- data.frame(transcript_ID = transcript_ID, count)
count <- merge( annote, count, by = "transcript_ID" )

count <- count[with(count, order(gene_ID,transcript_ID)),]
write.table(count, "salmon_aln-tx-counts.txt", sep="\t", row.names =F, quote = F)

count.gene <- aggregate(count[,-c(1:6)], 
	by=list(gene_ID=count$gene_ID, gene_name=count$gene_name), FUN=sum)
write.table(count.gene, "salmon_aln-gene-counts.txt", sep="\t", row.names =F, quote = F)
