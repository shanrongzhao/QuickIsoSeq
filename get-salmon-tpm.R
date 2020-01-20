################################################################################
#
# merge counts (and annotation) from salmon 
#
# Usage
#	      Rscript get-salmon-tpm.R allIDs.txt transcript_annotation_file
#
# The output:  
#			salmon-tx-tpm.txt
#			salmon-gene-tpm.txt
#			salmon-tx-counts.txt
#			salmon-gene-counts.txt
#
# Input
#Name	Length	EffectiveLength	TPM	NumReads
#ENST00000456328.2	1657	1493.04	0.228311	18.3346
#ENST00000450305.2	632	468.476	0.076341	1.92361
#ENST00000488147.1	1351	1187.04	14.1208	901.568
#
# Chi Zhang, Shanrong Zhao
#
# November 4, 2018
#
################################################################################

args <- commandArgs(trailingOnly = TRUE)
len<-length(args)
if ( len < 2) {
	print ("Please input <id file> and <transcript annotation file> ")
	q(save = "no")
}


manifest.file <- args[1]
manifest <- read.table(manifest.file, header=F)
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
	filename <- paste(id, "/salmon/quant.sf", sep="")
	
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
tpm <- round(as.data.frame( do.call("cbind", tpm)), digits=2)
tpm <- cbind(data.frame(transcript_ID = transcript_ID), tpm)
tpm <- merge( annote, tpm, by = "transcript_ID" )

tpm <- tpm[with(tpm, order(gene_ID,transcript_ID)),]
write.table(tpm, "salmon-tx-tpm.txt", sep="\t", row.names =F, quote = F)

tpm.gene <- aggregate(tpm[,-c(1:6)], 
	by=list(gene_ID=tpm$gene_ID, gene_name=tpm$gene_name), FUN=sum)
write.table(tpm.gene, "salmon-gene-tpm.txt", sep="\t", row.names =F, quote = F)

#
# count
#
count <- round(as.data.frame( do.call("cbind", count)), digits=1)
count <- data.frame(transcript_ID = transcript_ID, count)
count <- merge( annote, count, by = "transcript_ID" )

count <- count[with(count, order(gene_ID,transcript_ID)),]
write.table(count, "salmon-tx-counts.txt", sep="\t", row.names =F, quote = F)

count.gene <- aggregate(count[,-c(1:6)], 
	by=list(gene_ID=count$gene_ID, gene_name=count$gene_name), FUN=sum)
write.table(count.gene, "salmon-gene-counts.txt", sep="\t", row.names =F, quote = F)



