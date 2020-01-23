################################################################################
#
# get the summary counts and tpm table from different tools 
#
# Usage
#		merge-tpm-table.sh <ID file> <run configure file> <Output folder>[optional]
#
# Output:
#			*tpm.txt
#			*counts.txt
#			*rpkm.txt
#
# Shanrong Zhao
#
# December 12th, 2018
#
################################################################################

if [[ $# -ne 2 && $# -ne 3 ]]; then
    echo "$0 <ID file> <run configure file> <Output folder>[optional]"
	exit
fi

if [ -z "$QuickIsoSeq" ]; then
	echo "QuickIsoSeq environment variable has not been set"
	exit
fi

FLIST=$1
if [ ! -f $FLIST ]
then
	echo "$1 does not exist"
	exit
fi

CONFIG=$2
if [ -f $CONFIG ]
then
	source $CONFIG
else
	echo "The configuration file $CONFIG does not exist"
	exit
fi

if [ -z "$TRANSCRIPT_ANNOTATION" ]; then
	echo "TRANSCRIPT_ANNOTATION not set in the configure file"
	exit
fi

if [ ! -f $TRANSCRIPT_ANNOTATION ]
then
	echo "$TRANSCRIPT_ANNOTATION does not exist"
	exit
fi

#
# Merge featureCounts results
#
if [ -z "$GENE_ANNOTATION" ]; then
	echo "GENE_ANNOTATION not set in the configure file"
	exit
fi

if [ ! -f $GENE_ANNOTATION ]
then
	echo "$GENE_ANNOTATION does not exist"
	exit
fi

#featureCount
Rscript $QuickIsoSeq/get-fc-gene.R $1 $GENE_ANNOTATION


if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "KALLISTO" ]]; then
	Rscript $QuickIsoSeq/get-kallisto-tpm.R $1 $TRANSCRIPT_ANNOTATION
fi

if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "RSEM" ]]; then
	Rscript $QuickIsoSeq/get-rsem-tpm.R $1 $TRANSCRIPT_ANNOTATION
fi

if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "SALMON_ALN" ]]; then
	Rscript $QuickIsoSeq/get-salmon_aln-tpm.R $1 $TRANSCRIPT_ANNOTATION
fi

if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "SALMON_ALN" ]]; then
	Rscript $QuickIsoSeq/get-salmon-tpm.R $1 $TRANSCRIPT_ANNOTATION
fi


#
# move tables to output folder
#
RESULT_FOLDER=Results
if [[ $# -eq 3 ]]; then
    RESULT_FOLDER=$3

	mkdir -p $RESULT_FOLDER
	mv *tpm.txt $RESULT_FOLDER
	mv *counts.txt $RESULT_FOLDER
	mv *rpkm.txt $RESULT_FOLDER
fi

