############################################################################
##  1. map reads using STAR
##  2. SNP call from aligned BAM files
##  3. quantify transcript expression level
##  
##  Usage:
##  	run-isoseq.sh allIDs.txt run.config
##
##  Shanrong Zhao (shanrong.zhao@pfizer.com)
##
##  December 26, 2019
##
############################################################################


if [[ $# -ne 2 ]]; then
    echo "$0 <ID file> <run configure file>"
	exit
fi

FLIST=$1
CONFIG=$2

if [ ! -f $FLIST ]
then
	echo "$1 does not exists"
	exit
fi

if [ -f $CONFIG ]
then
	source $CONFIG
	FASTQ_DIR=$(readlink -f $FASTQ_DIR)
else
	echo "Please provide a configuration file"
	exit
fi


##############################################################################
#
# STAR Mapping
#
##############################################################################

# fq.gz or fastq.gz or fastq  or fq
READCMD="--readFilesCommand zcat"
if [[ $FASTQ_SUFFIX = "fastq" || $FASTQ_SUFFIX = "fq" || $FASTQ_SUFFIX = "fasta" || $FASTQ_SUFFIX = "fa" ]]; then
     READCMD=""
fi

#
#To handle very deep sequencing samples
#

#By default: --limitBAMsortRAM=0, will use genome size as RAM limit for BAM sorting
# SRR821518.star.err: EXITING because of fatal ERROR: not enough memory for BAM sorting:
# SOLUTION: re-run STAR with at least --limitBAMsortRAM 32847009527
# Below: --limitBAMsortRAM was set to 50GB. In case there is still a problem, please 
# set SEQUENCE_DEPTH=regular in yout run.config and re-run your analysis

limitBAMsortRAM="--limitBAMsortRAM 53687091200"
CORE=8
MEMORY_USAGE=""
M_PARAMETER=""
if [[ $SEQUENCE_DEPTH = "deep" ]]; then
	CORE=12
    MEMORY_USAGE="rusage [mem=131072]"
	M_PARAMETER="-M 131072"
	limitBAMsortRAM="--limitBAMsortRAM 107374182400"
fi

for f in `cat $FLIST`
do
   mkdir -p $f
   cd $f
   
   # handle pair or single RNA-seq
	FASTQS="$FASTQ_DIR/${f}_1.$FASTQ_SUFFIX $FASTQ_DIR/${f}_2.$FASTQ_SUFFIX"
	if [[ $SEQUENCE_TYPE = "single" ]]; then
		FASTQS="$FASTQ_DIR/${f}.$FASTQ_SUFFIX"
	fi

	#transcript-based BAM files are usually very large. Generate them only when they are needed, such as by RSEM
	QUANT_MODE=""
	RENAME_TX_BAM=""
	if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "RSEM" || $ISOFORM_ALGORITHM = "SALMON_ALN" ]]; then
		QUANT_MODE="--quantMode TranscriptomeSAM"
		RENAME_TX_BAM="mv Aligned.toTranscriptome.out.bam $f.tx.bam;"
	fi
	bsub -app large -J "$f.map-star" -n $CORE -R "span[hosts=1] $MEMORY_USAGE" $M_PARAMETER -o "$LOGDIR/$f.star.log" -e "$LOGDIR/$f.star.err" \
        "STAR --genomeDir $STAR_INDEX --readFilesIn $FASTQS $READCMD --runThreadN $CORE $limitBAMsortRAM --alignEndsType EndToEnd $QUANT_MODE --outSAMtype BAM SortedByCoordinate $STAR_PARAMETER; 
         mv Aligned.sortedByCoord.out.bam $f.sort.bam; $RENAME_TX_BAM 
		 samtools index $f.sort.bam;
         rm -rf _STARtmp"
   cd ..
done


###############################################################################
#
# featureCounts counting
#
###############################################################################
PAIRS="-p -B"
if [[ $SEQUENCE_TYPE = "single" ]]; then
        PAIRS=""
fi

for f in `cat $FLIST`
do
	cd $f
	bsub -app medium -n 6 -J "$f.featureCounts" -w "ended('$f.map-star')"  -o "$LOGDIR/$f.featureCounts.log" -e "$LOGDIR/$f.featureCounts.err" \
       		"featureCounts $PAIRS -T 6 -F GTF -a $GTF_FILE -t exon -g gene_id -s $STRAND -C $FEATURECOUNTS_OVERLAP -o ${f}.featureCounts.txt $f.sort.bam"
	cd ..
done


#############################################################################
#
# SNP CALL by VarScan
#
#############################################################################
for f in `cat $FLIST`
do
    cd $f
     bsub -J "$f.snp-call" -w "ended('$f.map-star')"  -M 12000 -R "rusage[mem=12000]" -o "$LOGDIR/$f.snp-call.log" -e "$LOGDIR/$f.snp-call.err" \
     "samtools mpileup -r $CHR_REGION -f $GENOME_FASTA $f.sort.bam | awk '\$4 != 0' | java -jar $VARSCAN_JAR pileup2snp - $VARSCAN_PARAMETER > $f.snp.txt"
    cd ..
done


##############################################################################
#
# Transcript quantification by one or more algorithms
#
#############################################################################

#
#COUNTING using RSEM
#
if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "RSEM" ]]; then

	RSEM_IS_PAIR="--paired-end"
	if [[ $SEQUENCE_TYPE = "single" ]]; then
		RSEM_IS_PAIR=" "
	fi
	
	#From RSEM 1.3.0, introduced `--strandedness <none|forward|reverse>` option, 
	#`--strand-specific` and `--forward-prob` are deprecated (but still supported). 
	#	
	# STRANDNESS="none"      #non-stranded
	# if [ $STRAND -eq 1 ]; then
		# STRANDNESS="forward"
	# fi
	# if [ $STRAND -eq 2 ]; then
		# STRANDNESS="reverse"
	# fi
	
	PROB=0.5      #non-stranded
	if [ $STRAND -eq 1 ]; then
		PROB=1
	fi
	if [ $STRAND -eq 2 ]; then
		PROB=0
	fi

	for f in `cat $FLIST`
	do
		cd $f
		#RSEM is typically slow, and the default job queue (express) does not always work.
		# The short job queue is recommeded and used
		bsub -q short -app medium -n 8 -J "$f.rsem" -w "ended('$f.map-star')"  -o "$LOGDIR/$f.rsem.log" -e "$LOGDIR/$f.rsem.err" \
		"rsem-calculate-expression -p 8 $RSEM_PARAMETER $RSEM_IS_PAIR  --forward-prob=$PROB  $f.tx.bam $RSEM_INDEX $f.rsem;
		mkdir -p rsem; mv -f -t rsem $f.rsem.* "
		cd ..
	done
fi


#
#COUNTING using salmon_aln
#
if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "SALMON_ALN" ]]; then
	SALMON_LIBTYPE="IU "      #non-stranded
	if [[ $SEQUENCE_TYPE = "single" ]]; then
		SALMON_LIBTYPE="U"
	fi
	if [ $STRAND -eq 1 ]; then
		SALMON_LIBTYPE="ISF"
	fi
	if [ $STRAND -eq 2 ]; then
		SALMON_LIBTYPE="ISR"
	fi
	for f in `cat $FLIST`
	do
	cd $f
		bsub -app medium -n 8 -J "$f.salmon_aln" -w "ended('$f.map-star')" -o "$LOGDIR/$f.salmon_aln.log" -e "$LOGDIR/$f.salmon_aln.err" \
		"salmon quant -t $TRANSCRIPT_FASTA -l $SALMON_LIBTYPE -a $f.tx.bam -p 8 -o salmon_aln $SALMON_ALN_PARAMETER"
	cd ..
	done
fi


#
#COUNTING using salmon
#
if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "SALMON" ]]; then
	SALMON_LIBTYPE="IU "      #non-stranded
	if [[ $SEQUENCE_TYPE = "single" ]]; then
		SALMON_LIBTYPE="U"
	fi
	if [ $STRAND -eq 1 ]; then
		SALMON_LIBTYPE="ISF"
	fi
	if [ $STRAND -eq 2 ]; then
		SALMON_LIBTYPE="ISR"
	fi
	for f in `cat $FLIST`
	do
	cd $f
		FASTQS="-1 $FASTQ_DIR/${f}_1.$FASTQ_SUFFIX -2 $FASTQ_DIR/${f}_2.$FASTQ_SUFFIX"
		if [[ $SEQUENCE_TYPE = "single" ]]; then
			FASTQS="-r $FASTQ_DIR/${f}.$FASTQ_SUFFIX"
		fi
		bsub -app medium -n 8 -J "$f.salmon" -o "$LOGDIR/$f.salmon.log" -e "$LOGDIR/$f.salmon.err" \
		"salmon quant -i $SALMON_INDEX -l $SALMON_LIBTYPE $FASTQS -p 8 -o salmon $SALMON_PARAMETER"
	cd ..
	done
fi


#
#COUNTING using kallisto
#
if [[ $ISOFORM_ALGORITHM = "ALL" || $ISOFORM_ALGORITHM = "KALLISTO" ]]; then
	KALLISTO_STRAND=" "      #non-stranded
	if [ $STRAND -eq 1 ]; then
		KALLISTO_STRAND="--fr-stranded"
	fi
	if [ $STRAND -eq 2 ]; then
		KALLISTO_STRAND="--rf-stranded"
	fi
	
	for f in `cat $FLIST`
	do
	cd $f
		FASTQS="$FASTQ_DIR/${f}_1.$FASTQ_SUFFIX $FASTQ_DIR/${f}_2.$FASTQ_SUFFIX"
		if [[ $SEQUENCE_TYPE = "single" ]]; then
			FASTQS="--single -l 190 -s 30 $FASTQ_DIR/${f}.$FASTQ_SUFFIX"
		fi
		bsub -app medium -n 8 -J "$f.kallisto" -o "$LOGDIR/$f.kallisto.log" -e "$LOGDIR/$f.kallisto.err" \
		"kallisto quant -i $KALLISTO_INDEX $KALLISTO_STRAND $KALLISTO_PARAMETER -t 8 -o kallisto $FASTQS; 
		 kallisto h5dump -o kallisto ./kallisto/abundance.h5"
	cd ..
	done
fi


#
# split sequences into rRNA, hgRNA and other
#
# This is available for human RNA samples, and especially useful for blood RNA-seq
#
if [[ $SAMPLE_SPECIES = "human" ]]; then

	for f in `cat $FLIST`
	do
		cd $f
		
		fastq_file=$FASTQ_DIR/${f}.${FASTQ_SUFFIX}
		if [ ! -f $fastq_file ]; then
			#use Read 1 for pair ended RNA-seq
			fastq_file=$FASTQ_DIR/${f}_1.${FASTQ_SUFFIX} 
		fi
		
		# https://github.com/BenLangmead/bowtie/issues/31
		# .gz is supported by bowtie from v1.2.1. This is a fantastic news.
		# unfortunately, this feature is not documented by bowtie user manual.
		
		bsub -app medium -n 8 -J "$f.rRNA" -o "$LOGDIR/$f.rRNA.log" -e "$LOGDIR/$f.rRNA.err" \
			"bowtie -p 8 -v 2 $rRNA_BWT_INDEX $fastq_file 1>/dev/null 2>${f}.bwt.rRNA.log"
			
		bsub -app medium -n 8 -J "$f.hgRNA" -o "$LOGDIR/$f.hgRNA.log" -e "$LOGDIR/$f.hgRNA.err" \
			"bowtie -p 8 -v 2 $hgRNA_BWT_INDEX $fastq_file 1>/dev/null 2>${f}.bwt.hgRNA.log"
		
		cd ..
	done
fi
