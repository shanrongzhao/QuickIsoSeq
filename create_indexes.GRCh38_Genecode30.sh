##################################################################
#
# Download all required files from Gencode and prepare indexes 
# required for QuickIsoSeq pipeline
#
# Shanrong Zhao
#
# December 24, 2019
#
##################################################################

#
#export QuickIsoSeq=/lustre/workspace/projects/ECD/zhaos25/QuickIsoSeq_Pub

if [ -z "$QuickIsoSeq" ]; then
	echo "QuickIsoSeq environment variable has not been set"
	exit
fi

#
# Set the right paths to all executables.
# $QuickIsoSeq/tools_path.sh was generated when you run install-tools.sh
#
source $QuickIsoSeq/Tools/tools_path.txt

#
# It is recommended to name index_root as refGenone_Annotation
#
SPECIES=human

GENCODE_VERSION=30

INDEX_ROOT=$QuickIsoSeq/Indexes/GRCh38_Genecode30

mkdir -p $INDEX_ROOT
cd $INDEX_ROOT
rm -rf *

#
# you need to download genome(in fasta) and annotation (in GTF) 
# Either in advance. If so, please set GENOME_FASTA and GTF_FILE properly
# Or download them on the fly
#

# GENOME_FASTA=$QuickIsoSeq/Indexes/human/fasta/GRCh38.primary.genome.fa
# GTF_FILE=$QuickIsoSeq/Indexes/human/annotation/gencode.v30.gtf
# cp $GENOME_FASTA genome.fa
# cp $GTF_FILE gene.gtf

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa genome.fa

GTF_FTP=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf.gz
wget $GTF_FTP
gunzip gencode.v${GENCODE_VERSION}.annotation.gtf.gz 
mv gencode.v${GENCODE_VERSION}.annotation.gtf gene.gtf


#
# Extract transcript sequences from reference genome/tranascriptome
# Extract annotations from annotations
#
samtools faidx  genome.fa
gffread -w transcript.fa  -g genome.fa gene.gtf
gtf2annot.pl  gene.gtf > gene.annot
gtf2annotTx.pl gene.gtf > transcript.annot


#
# create indexes
#
STAR_INDEX=${INDEX_ROOT}/STAR_100
STAR_INDEX_150=${INDEX_ROOT}/STAR_150
RSEM_INDEX=${INDEX_ROOT}/rsem/rsem
SALMON_INDEX=${INDEX_ROOT}/salmon
KALLISTO_INDEX=${INDEX_ROOT}/kallisto/index

mkdir -p $STAR_INDEX
mkdir -p $STAR_INDEX_150
mkdir -p $RSEM_INDEX
mkdir -p $SALMON_INDEX
mkdir -p ${INDEX_ROOT}/kallisto


LOGDIR=/hpc/grid/scratch/${USER}/log
mkdir -p $LOGDIR

#STAR INdex
bsub -app large -n 8 -R "span[hosts=1]" -J "star_index" -o "$LOGDIR/star_index.log" -e "$LOGDIR/star_index.err" "STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $STAR_INDEX --genomeFastaFiles genome.fa --sjdbGTFfile gene.gtf --sjdbOverhang 99; rm -rf _STARtmp Log.out"
bsub -app large -n 8 -R "span[hosts=1]" -J "star_index2" -o "$LOGDIR/star_index2.log" -e "$LOGDIR/star_index2.err" "STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $STAR_INDEX_150 --genomeFastaFiles genome.fa --sjdbGTFfile gene.gtf --sjdbOverhang 149; rm -rf _STARtmp Log.out"

#SALMON Index
bsub -app large -R "span[hosts=1]" -n 8 -J "salmon_index" -o "$LOGDIR/salmon_index.log" -e "$LOGDIR/salmon_index.err" "salmon index -t transcript.fa -p 8 -i $SALMON_INDEX"

#KALLISTO Index
bsub -app large -R "span[hosts=1]" -J "kal_index" -o "$LOGDIR/kal_index.log" -e "$LOGDIR/kal_index.err" "kallisto index --index=$KALLISTO_INDEX transcript.fa" 

#RSEM index
bsub -app medium -R "span[hosts=1]" -J "RSEM-index" -o "$LOGDIR/RSEM-index.log" -e "$LOGDIR/RSEM-index.err" "rsem-prepare-reference --gtf gene.gtf genome.fa $RSEM_INDEX"


#
# hgRNA and rRNA check is performed only for human samples
#
if [ $SPECIES = "human" ]; then
	#
	# hgRNA
	#
	$QuickIsoSeq/extract-hgRNA.v2.pl transcript.fa transcript.annot hgRNAs.fa

	#
	# rRNA
	# https://software.broadinstitute.org/cancer/cga/rnaseqc_download
	#
	wget https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/rnaseqc/rRNA.tar.gz
	tar -xzvf rRNA.tar.gz
	rm -rf  rRNA.tar.gz human_all_rRNA.fasta.*
	mv human_all_rRNA.fasta rRNA.fa
	touch rRNA.fa

	BOWTIE_rRNA_INDEX=${INDEX_ROOT}/bowtie_rRNA
	BOWTIE_hgRNA_INDEX=${INDEX_ROOT}/bowtie_hgRNA
	mkdir -p $BOWTIE_rRNA_INDEX
	mkdir -p $BOWTIE_hgRNA_INDEX

	#
	# rRNA index for bowtie
	#
	bsub -J "Bowtie-index-rRNA" -o "$LOGDIR/rRNA-index.log" -e "$LOGDIR/rRNA-index.err" \
		"bowtie-build rRNA.fa $BOWTIE_rRNA_INDEX/rRNA"

	#
	# hgRNA index for bowtie
	#
	bsub -J "Bowtie-index-hgRNA" -o "$LOGDIR/hgRNA-index.log" -e "$LOGDIR/hgRNA-index.err" \
		"bowtie-build hgRNAs.fa $BOWTIE_hgRNA_INDEX/hgRNA"
fi

#
# Add the following output to run.config
#
# please go to STAR index fold and find out the coordinates 
# for the chromosome where MHC genes are located
# 
# For human: chr6;  for mouse:  chr17
#
echo "
Note: The default CHR_REGION parameter is correct only for huamn GRCh38
1. Please go to $STAR_INDEX 
  folder and find out the coordinates for the chromosome
  where MHC genes are located and then update CHR_REGION
  in both run.config 
  and $INDEX_ROOT/indexes_path.txt
2. Please add the following output to run.config

"

echo "
SAMPLE_SPECIES=$SPECIES

INDEX_ROOT=$INDEX_ROOT

GENOME_FASTA=\$INDEX_ROOT/genome.fa
GTF_FILE=\$INDEX_ROOT/gene.gtf
TRANSCRIPT_FASTA=\$INDEX_ROOT/transcript.fa
TRANSCRIPT_ANNOTATION=\$INDEX_ROOT/transcript.annot
GENE_ANNOTATION=\$INDEX_ROOT/gene.annot
CHR_REGION=chr6:1-170805979

STAR_INDEX=\${INDEX_ROOT}/STAR_100
RSEM_INDEX=\${INDEX_ROOT}/rsem/rsem
SALMON_INDEX=\${INDEX_ROOT}/salmon
KALLISTO_INDEX=\${INDEX_ROOT}/kallisto/index"

#
# Record indexes paths to $QuickIsoSeq/indexes_path.txt
# You can insert those indexes paths to yuor run.config
#
echo "
SAMPLE_SPECIES=$SPECIES

INDEX_ROOT=$INDEX_ROOT

GENOME_FASTA=\$INDEX_ROOT/genome.fa
GTF_FILE=\$INDEX_ROOT/gene.gtf
TRANSCRIPT_FASTA=\$INDEX_ROOT/transcript.fa
TRANSCRIPT_ANNOTATION=\$INDEX_ROOT/transcript.annot
GENE_ANNOTATION=\$INDEX_ROOT/gene.annot
CHR_REGION=chr6:1-170805979

STAR_INDEX=\${INDEX_ROOT}/STAR_100
RSEM_INDEX=\${INDEX_ROOT}/rsem/rsem
SALMON_INDEX=\${INDEX_ROOT}/salmon
KALLISTO_INDEX=\${INDEX_ROOT}/kallisto/index
"  > $INDEX_ROOT/indexes_path.txt


if [ $SPECIES = "human" ]; then
echo "
rRNA_BWT_INDEX=\${INDEX_ROOT}/bowtie_rRNA/rRNA
hgRNA_BWT_INDEX=\${INDEX_ROOT}/bowtie_hgRNA/hgRNA
"
echo "
rRNA_BWT_INDEX=\${INDEX_ROOT}/bowtie_rRNA/rRNA
hgRNA_BWT_INDEX=\${INDEX_ROOT}/bowtie_hgRNA/hgRNA
" >> $INDEX_ROOT/indexes_path.txt
fi
