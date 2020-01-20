#!/bin/bash
#
# map RNA-seq to rRNA
#
# check-RNA.sh allIDs.txt <FASTQ_DIR> fastq.gz $Indexes/human/bowtie_rRNA/rRNA rRNA
#  eg: 
#  export PATH=/lustre/workspace/projects/ECD/zhaos25/QuickIsoSeq/Tools/bowtie-1.2.2/:$PATH
#
#  For ribosomal RNA:
#  check-RNA.sh allIDs.txt /lustre/workspace/projects/ECD/zhaos25/GTEx fastq.gz  
#	/lustre/workspace/projects/ECD/zhaos25/QuickIsoSeq/Indexes/human/bowtie_rRNA/rRNA rRNA
#
#  For hemoglobin:
#  check-RNA.sh allIDs.txt /lustre/workspace/projects/ECD/zhaos25/GTEx fastq.gz  
#	/lustre/workspace/projects/ECD/zhaos25/QuickIsoSeq/Indexes/human/bowtie_hgRNA/hgRNA hgRNA
#
# Shanrong Zhao
#
# December 24, 2018
#

if [[ $# -ne 5 ]]; then
    echo "check-RNA.sh <ID file> <FASTQ_DIR> <SUFFIX> <Bowtie index> <category>"
	echo ""
	echo "category: rRNA or hgRNA"
	echo ""
	echo "For exmaple:"
	echo "check-RNA.sh allIDs.txt /lustre/workspace/projects/ECD/zhaos25/GTEx fastq.gz /lustre/workspace/projects/ECD/zhaos25/QuickIsoSeq/Indexes/human/bowtie_rRNA/rRNA"
	exit
fi

FLIST=$1
FASTQ_DIR=$2
#fastq suffix: fq, fq.gz, fastq, fastq.gz
SUFFIX=$3
BWT_INDEX=$4
CATEGORY=$5
echo $CATEGORY
exit

LOGDIR=/hpc/grid/scratch/$USER/log
mkdir -p  $LOGDIR


for f in `cat $FLIST`
do
	mkdir -p $f
	cd $f

	fastq_file=$FASTQ_DIR/${f}.${SUFFIX}
	if [ ! -f $fastq_file ]; then
		#use Read 1 for pair ended RNA-seq
		fastq_file=$FASTQ_DIR/${f}_1.${SUFFIX} 
	fi
	
	#https://github.com/BenLangmead/bowtie/issues/31
	# .gz is supported by bowtie from v1.2.1. This is a fantastic news.
	#unfortunately, this feature is not documented by bowtie user manual.
	
	bsub -app medium -n 8 -J "$f.${CATEGORY}" -o "$LOGDIR/$f.${CATEGORY}.log" -e "$LOGDIR/$f.${CATEGORY}.err" \
		"bowtie -p 8 -v 2 $BWT_INDEX $fastq_file 1>/dev/null 2>${f}.bwt.${CATEGORY}.log"
		
	cd ..
done


