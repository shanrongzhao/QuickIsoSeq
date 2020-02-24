#######################################################################
#
# Usage
#       merge-isoseq.sh allIDs.txt run.config output_directory[optional, default to Results]
#
# Output:
#			star-mapping-summary*
#			fc-counting-summary*
#			RNASeq-lib-size
#			RNASeq-snp-corr* 
#			*tpm.txt
#			*counts.txt
#			*rpkm.txt
#			expr-count-*
#			expr-gene-*
#
#
# Author: Shanrong Zhao
#
# December 23, 2019
#
# Note: run this script under the project root directory
#
#######################################################################

if [[ $# -ne 2 && $# -ne 3 ]]; then
    echo "$0 <ID file> <run configure file> <Output folder>[optional]"
	echo "Default output folder:  Results"
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

#
# Extract annotation from $SAMPLE_ANNOTATION
#
if [ -n "$SAMPLE_ANNOTATION" ]; then
	dos2unix $SAMPLE_ANNOTATION
	Rscript $QuickIsoSeq/extractBy.R $SAMPLE_ANNOTATION $FLIST sample.annot.tmp
	SAMPLE_ANNOTATION=sample.annot.tmp
fi


#
# calculate most optimal width/heigth for plots
#
no_sample=`cat $FLIST | wc -l | awk '{print $1}'`
width=$(expr $no_sample / 8 + 1 )
if [[ $width -lt 4 ]]; then
	width=4
fi
height=$(expr $width - 1)


# 
# Part A: Merging mapping, counting, and SNP
#
echo ""
echo "Part A: Merge mapping, counting, and SNP"

get-star-summary.pl $FLIST  star-mapping-summary.txt
get-fc-summary.pl $FLIST  fc-counting-summary.txt

get-snp-corr.pl $FLIST RNASeq-snp-corr.txt
if [ -n "$SAMPLE_ANNOTATION" ]; then
	Rscript $QuickIsoSeq/reorder-snp-corr.R RNASeq-snp-corr.txt "$SAMPLE_ANNOTATION"
fi

Rscript $QuickIsoSeq/plot-rnaseq-metrics.R star-mapping-summary.txt star-mapping-summary ${width}x6 Mapping

Rscript $QuickIsoSeq/plot-lib-size.R star-mapping-summary.txt RNASeq-lib-size ${width}x6

Rscript $QuickIsoSeq/plot-rnaseq-metrics.R fc-counting-summary.txt fc-counting-summary ${width}x6 Counting

Rscript $QuickIsoSeq/plot-corr-matrix.R RNASeq-snp-corr.txt RNASeq-snp-corr ${width}x${height}

if [[ $SAMPLE_SPECIES = "human" ]]; then
	get-reads-split.pl $FLIST reads-split-summary.txt
	Rscript $QuickIsoSeq/plot-rnaseq-metrics.R reads-split-summary.txt reads-split-summary ${width}x6 Reads_split
fi


# 
# Part B: Merging quantification results
#
#$QuickIsoSeq/merge-tpm-table.sh $FLIST $CONFIG
echo ""
echo "Part B: Merging quantification results"

if [ -z "$TRANSCRIPT_ANNOTATION" ]; then
	echo "TRANSCRIPT_ANNOTATION not set in the configure file"
	exit
fi
if [ ! -f $TRANSCRIPT_ANNOTATION ]
then
	echo "$TRANSCRIPT_ANNOTATION does not exist"
	exit
fi

if [ -z "$GENE_ANNOTATION" ]; then
	echo "GENE_ANNOTATION not set in the configure file"
	exit
fi
if [ ! -f $GENE_ANNOTATION ]
then
	echo "$GENE_ANNOTATION does not exist"
	exit
fi

# merge featureCounts results
Rscript $QuickIsoSeq/get-fc-gene.R $FLIST $GENE_ANNOTATION

if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "KALLISTO" ]]; then
	echo "Merging counts table from KALLISTO"
	Rscript $QuickIsoSeq/get-kallisto-tpm.R $FLIST $TRANSCRIPT_ANNOTATION
fi

if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "RSEM" ]]; then
	echo "Merging counts table from RSEM"
	Rscript $QuickIsoSeq/get-rsem-tpm.R $FLIST $TRANSCRIPT_ANNOTATION
fi

if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "SALMON_ALN" ]]; then
	echo "Merging counts table from SALMON_ALN"
	Rscript $QuickIsoSeq/get-salmon_aln-tpm.R $FLIST $TRANSCRIPT_ANNOTATION
fi

if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "SALMON" ]]; then
	echo "$Merging counts table from SALMON"
	Rscript $QuickIsoSeq/get-salmon-tpm.R $FLIST $TRANSCRIPT_ANNOTATION
fi


#
# Part C: correlation-based QC
#
echo ""
echo "Part C: correlation-based sample QC"

#
# only 1 algorithm is selected in the order
# salmon > salmon_aln > > kallisto > RSEM

TX_TPM_FILE=rsem-tx-tpm.txt 
if [[ $ISOFORM_ALGORITHM = "KALLISTO" ]]; then
	TX_TPM_FILE=kallisto-tx-tpm.txt
fi

if [[ $ISOFORM_ALGORITHM = "SALMON_ALN" ]]; then
	TX_TPM_FILE=salmon_aln-tx-tpm.txt
fi

if [[ $ISOFORM_ALGORITHM = "SALMON" || $ISOFORM_ALGORITHM = "ALL" ]]; then
	TX_TPM_FILE=salmon-tx-tpm.txt
fi

Rscript $QuickIsoSeq/plot-expr-count.R $TX_TPM_FILE expr-count ${width}x${height}

Rscript $QuickIsoSeq/plot-expr-qc.R  $TX_TPM_FILE expr-gene  ${width}x${height}


#
# Part D: merge QC metrics
#
echo ""
echo "Part D: merge QC metrics and move results to the output folder"

FILES="star-mapping-summary.txt fc-counting-summary.txt expr-gene-MADScore.txt expr-gene-top10x.txt"
if [[ $SAMPLE_SPECIES = "human" ]]; then
	FILES="reads-split-summary.txt star-mapping-summary.txt fc-counting-summary.txt expr-gene-MADScore.txt expr-gene-top10x.txt"
fi

if [ -n "$SAMPLE_ANNOTATION" ]; then
	Rscript $QuickIsoSeq/merge-QC-files.R "$SAMPLE_ANNOTATION" $FILES
else
	Rscript $QuickIsoSeq/merge-QC-files.R $FILES
fi


RESULT_FOLDER=Results
if [[ $# -eq 3 ]]; then
    RESULT_FOLDER=$3
fi
SUMMARY_FOLDER=$RESULT_FOLDER/Summary
rm -rf $RESULT_FOLDER

echo "The summary results are under $SUMMARY_FOLDER"


mkdir -p $RESULT_FOLDER
mkdir -p $SUMMARY_FOLDER

if [ -n "$SAMPLE_ANNOTATION" ]; then
	mv "$SAMPLE_ANNOTATION" $SUMMARY_FOLDER/sample.annotation.txt
fi

cp $CONFIG $SUMMARY_FOLDER

mv *-summary* $SUMMARY_FOLDER

mv RNASeq-lib-size* $SUMMARY_FOLDER
mv RNASeq-snp-corr* $SUMMARY_FOLDER
mv RNASeq-merged-metrics.txt $SUMMARY_FOLDER

mv *tpm.txt $SUMMARY_FOLDER
mv *counts.txt $SUMMARY_FOLDER
mv *rpkm.txt $SUMMARY_FOLDER

mv expr-count* $SUMMARY_FOLDER
mv expr-gene*  $SUMMARY_FOLDER


#
# Part E: multiqc
#

command_exists () {
    type "$1" &> /dev/null ;
}
if command_exists multiqc; then
    echo ""
	echo "Part E: run multiqc and the multiqc report multiqc_report.html is under $RESULT_FOLDER"
	multiqc -f -m star -m featureCounts -m bowtie1 -m fastqc .
	mv multiqc_report.html $RESULT_FOLDER
	mv multiqc_data $RESULT_FOLDER
else
	echo "multiqc (a python-based package) is recommeded but not installed"
fi

