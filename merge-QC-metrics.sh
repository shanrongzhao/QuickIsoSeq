#######################################################################
#
# Usage
#       merge-QC-metrics.sh allIDs.txt output_directory[optional, default to Results]
#
# Output:
#			star-mapping-summary*
#			fc-counting-summary*
#			RNASeq-snp-corr* 
#
#
# Author: Shanrong Zhao
#
# October 12, 2018
#
# Note: run this script under the project root directory
#
#######################################################################

if [[ $# -ne 1 && $# -ne 2 ]]; then
    echo "$0 <ID file> output_directory[optional, default to Results]"
	exit
fi


if [ -z "$QuickIsoSeq" ]; then
	echo "QuickIsoSeq environment variable has not been set"
	exit
fi

if [ ! -f $1 ]; then
	echo "$1 does not exist"
	exit
fi

FLIST=$1

# 
# Part A: Merging mapping, counting, and SNP
#
get-star-summary.pl $FLIST  star-mapping-summary.txt
get-fc-summary.pl $FLIST  fc-counting-summary.txt
get-snp-corr.pl $FLIST RNASeq-snp-corr.txt

#
no_sample=`cat $FLIST | wc -l | awk '{print $1}'`
width=$(expr $no_sample / 8 + 1 )
if [[ $width -lt 5 ]]; then
	width=5
fi
height=$(expr $width - 1)

Rscript $QuickIsoSeq/plot-rnaseq-metrics.R star-mapping-summary.txt star-mapping-summary ${width}x6

Rscript $QuickIsoSeq/plot-rnaseq-metrics.R fc-counting-summary.txt fc-counting-summary ${width}x6

Rscript $QuickIsoSeq/plot-corr-matrix.R RNASeq-snp-corr.txt RNASeq-snp-corr ${width}x${height}


#
RESULT_FOLDER=Results
if [[ $# -eq 2 ]]; then
    RESULT_FOLDER=$2
	mkdir -p $RESULT_FOLDER

	mv star-mapping-summary* $RESULT_FOLDER
	mv fc-counting-summary* $RESULT_FOLDER
	mv RNASeq-snp-corr* $RESULT_FOLDER
fi
