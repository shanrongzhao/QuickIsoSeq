################################################################################
#
# merge counts (and annotation) from kallisto 
#
# Usage
#	Rscript get-kallisto-tpm.R allIDs.txt transcript_annotation_file 
#
# The output:  
#			kallisto-tx-tpm.txt
#			kallisto-gene-tpm.txt
#			kallisto-tx-counts.txt
#			kallisto-gene-counts.txt
#
# Input
#target_id	length	eff_length	est_counts	tpm
#ENST00000456328.2	1657	1495.73	31.5647	0.37312
#ENST00000450305.2	632	470.766	1.80423	0.0677622
#ENST00000488147.1	1351	1189.73	866.902	12.8831
#
# Shanrong Zhao
#
# November 4, 2018
#
################################################################################
args <- commandArgs(trailingOnly = TRUE)
len<-length(args)
if ( len < 2) {
	print ("Please input <id file> and <transcript annotation file>")
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

#genes
tpm <- list()
count <- list()
i <- 0
for (id in samples) {
	filename <- paste(id, "/kallisto/abundance.tsv", sep="")
	
	data <- read.table(filename, header=TRUE, stringsAsFactors=F, check.names = F)

	tpm[[id]] <- data[,5]
	count[[id]] <- data[,4]
	
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
write.table(tpm, "kallisto-tx-tpm.txt", sep="\t", row.names =F, quote = F)

tpm.gene <- aggregate(tpm[,-c(1:6)], 
	by=list(gene_ID=tpm$gene_ID, gene_name=tpm$gene_name), FUN=sum)
write.table(tpm.gene, "kallisto-gene-tpm.txt", sep="\t", row.names =F, quote = F)

#
# count
#
count <- round(as.data.frame( do.call("cbind", count)), digits=1)
count <- data.frame(transcript_ID = transcript_ID, count)
count <- merge( annote, count, by = "transcript_ID" )

count <- count[with(count, order(gene_ID,transcript_ID)),]
write.table(count, "kallisto-tx-counts.txt", sep="\t", row.names =F, quote = F)

count.gene <- aggregate(count[,-c(1:6)], 
	by=list(gene_ID=count$gene_ID, gene_name=count$gene_name), FUN=sum)
write.table(count.gene, "kallisto-gene-counts.txt", sep="\t", row.names =F, quote = F)

