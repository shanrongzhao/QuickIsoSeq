################################################################################
#
# merge counts (and annotation) from RSEM 
#
# Usage
#	Rscript get-rsem-tpm.R allIDs.txt transcript_annotation_file
#
# The output:  
#			rsem-tx-tpm.txt
#			rsem-gene-tpm.txt
#			rsem-tx-counts.txt
#			rsem-gene-counts.txt
#
# Input
#
#transcript_id   gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct
#ENST00000373020.4       ENSG00000000003.10      2206    2042.02 2561.79 26.07   19.18   97.93
#ENST00000494424.1       ENSG00000000003.10      820     656.02  8.45    0.27    0.20    1.01
#ENST00000496771.1       ENSG00000000003.10      1025    861.02  11.76   0.28    0.21    1.07
#
# Shanrong Zhao
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

#isoforms
tpm <- list()
count <- list()
i <- 0
for (id in samples) {
	filename = paste(id, "/rsem/", id, ".rsem.isoforms.results", sep="")
	
	data <- read.table(filename, header=TRUE, stringsAsFactors=F, check.names = F)
	#data <- data[with(data, order(Geneid)),]

	tpm[[id]] <- data[,6]
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
write.table(tpm, "rsem-tx-tpm.txt", sep="\t", row.names =F, quote = F)

tpm.gene <- aggregate(tpm[,-c(1:6)], 
	by=list(gene_ID=tpm$gene_ID, gene_name=tpm$gene_name), FUN=sum)
write.table(tpm.gene, "rsem-gene-tpm.txt", sep="\t", row.names =F, quote = F)

#
# count
#
count <- round(as.data.frame( do.call("cbind", count)), digits=1)
count <- data.frame(transcript_ID = transcript_ID, count)
count <- merge( annote, count, by = "transcript_ID" )

count <- count[with(count, order(gene_ID,transcript_ID)),]
write.table(count, "rsem-tx-counts.txt", sep="\t", row.names =F, quote = F)

count.gene <- aggregate(count[,-c(1:6)], 
	by=list(gene_ID=count$gene_ID, gene_name=count$gene_name), FUN=sum)
write.table(count.gene, "rsem-gene-counts.txt", sep="\t", row.names =F, quote = F)
