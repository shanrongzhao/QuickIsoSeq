#
#Below are command lines to install tools needed for QuickIsoSeq/Tools
#You can customize it to fit your environment. 
#
#Shanrong Zhao
#
#January 14, 2020
#

#export QuickIsoSeq=/lustre/workspace/projects/ECD/zhaos25/QuickIsoSeq_Pub

if [ -z "$QuickIsoSeq" ]; then
	echo "QuickIsoSeq environment variable has not been set"
	exit
fi

APPLICATION_ROOT=$QuickIsoSeq/Tools

mkdir -p $APPLICATION_ROOT
cd $APPLICATION_ROOT

rm -rf *

#
# STAR installation
#
# Please go to https://github.com/alexdobin/STAR and find out the 
# corresponding STAR version and change STAR_VERSION
STAR_VERSION=2.7.3a
wget https://github.com/alexdobin/STAR/archive/master.zip
unzip master.zip
rm -rf  master.zip
mv STAR-master STAR_$STAR_VERSION

#
# http://subread.sourceforge.net/
#
# Subread/FeatureCount installation. Please change SUBREAD_VERSION 
# if yuo download a different version
# 
#
SUBREAD_VERSION=subread-2.0.0
wget https://sourceforge.net/projects/subread/files/$SUBREAD_VERSION/${SUBREAD_VERSION}-Linux-x86_64.tar.gz
tar -xzf ${SUBREAD_VERSION}-Linux-x86_64.tar.gz
rm -rf ${SUBREAD_VERSION}-Linux-x86_64.tar.gz
mv ${SUBREAD_VERSION}-Linux-x86_64 ${SUBREAD_VERSION}

#
# https://github.com/dkoboldt/varscan
#
# VARSCAN
# Note: v2.4.0 is the latest version, and no newer version any more
# 
wget https://github.com/dkoboldt/varscan/raw/master/VarScan.v2.4.0.jar

#
# http://www.htslib.org/download/
#
# samtools
# 
SAMTOOLS_VERSION=1.9
wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
bzip2 -d samtools-${SAMTOOLS_VERSION}.tar.bz2
tar -xvf  samtools-${SAMTOOLS_VERSION}.tar
rm -rf samtools-${SAMTOOLS_VERSION}.tar

mv samtools-${SAMTOOLS_VERSION} samtools-${SAMTOOLS_VERSION}.source
mkdir samtools-${SAMTOOLS_VERSION}

cd samtools-${SAMTOOLS_VERSION}.source
./configure --prefix=$APPLICATION_ROOT/samtools-${SAMTOOLS_VERSION}
make
make install

cd ..

#
# https://deweylab.github.io/RSEM/
# https://github.com/deweylab/RSEM/releases
#
# RSEM
# 
RSEM_VERSIOM=v1.3.1
RSEM_PKG=RSEM-1.3.1

wget https://github.com/deweylab/RSEM/archive/${RSEM_VERSIOM}.zip
unzip ${RSEM_VERSIOM}.zip
rm -f ${RSEM_VERSIOM}.zip

cd $RSEM_PKG
make 
cd ..

#
# https://github.com/COMBINE-lab/salmon/releases
#
# SALMON
# 
SALMON_VERSION=v1.1.0
SALMON_PKG=salmon-1.1.0
wget https://github.com/COMBINE-lab/salmon/releases/download/$SALMON_VERSION/${SALMON_PKG}_linux_x86_64.tar.gz
tar -xzvf  ${SALMON_PKG}_linux_x86_64.tar.gz
rm -rf ${SALMON_PKG}_linux_x86_64.tar.gz
mv ${SALMON_PKG}_linux_x86_64 ${SALMON_PKG}
mv salmon-latest_linux_x86_64 ${SALMON_PKG}

#
# http://pachterlab.github.io/kallisto/download
#
# Kallisto
# 
KALLISTO_VERSION=v0.46.1
wget https://github.com/pachterlab/kallisto/releases/download/${KALLISTO_VERSION}/kallisto_linux-${KALLISTO_VERSION}.tar.gz
tar -xzvf kallisto_linux-${KALLISTO_VERSION}.tar.gz
rm -f kallisto_linux-${KALLISTO_VERSION}.tar.gz


#
# https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.3/
# https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.3/bowtie-1.2.3-linux-x86_64.zip/download
#
# Bowtie  (not Botie2)
# 
BOWTIE_VERSION=1.2.3
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/$BOWTIE_VERSION/bowtie-$BOWTIE_VERSION-linux-x86_64.zip
unzip bowtie-$BOWTIE_VERSION-linux-x86_64.zip
rm -f bowtie-$BOWTIE_VERSION-linux-x86_64.zip


#
# http://ccb.jhu.edu/software/stringtie/dl/
#
# gffread
# 
wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.11.4.Linux_x86_64.tar.gz
tar -xzvf gffread-0.11.4.Linux_x86_64.tar.gz
rm -f gffread-0.11.4.Linux_x86_64.tar.gz

#
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
#
# fastgc
# 
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
rm  fastqc_v0.11.9.zip
chmod +x FastQC/fastqc


#
# Add the following output to run.config
#
echo "
Please add the following exccutables and PATH to run.config file
	"

echo "
APPLICATION_ROOT=\$QuickIsoSeq/Tools

STAR=\$APPLICATION_ROOT/STAR_$STAR_VERSION/bin/Linux_x86_64_static
FEATURECOUNTS=\$APPLICATION_ROOT/${SUBREAD_VERSION}/bin
VARSCAN_JAR=\$APPLICATION_ROOT/VarScan.v2.4.0.jar
SAMTOOLS=\$APPLICATION_ROOT/samtools-${SAMTOOLS_VERSION}/bin
RSEM=\$APPLICATION_ROOT/$RSEM_PKG
SALMON=\$APPLICATION_ROOT/${SALMON_PKG}/bin
KALLISTO=\$APPLICATION_ROOT/kallisto
GFFREAD=\$APPLICATION_ROOT/gffread-0.11.4.Linux_x86_64
BOWTIE=\$APPLICATION_ROOT/bowtie-$BOWTIE_VERSION-linux-x86_64
FASTQC=\$APPLICATION_ROOT/FastQC
export PATH=\$SAMTOOLS:\$STAR:\$FEATURECOUNTS:\$RSEM:\$SALMON:\$KALLISTO:\$GFFREAD:\$BOWTIE:\$FASTQC:\$PATH
"

#
# Put exccutables and PATH to $QuickIsoSeq/tools_path.sh
#
echo "
APPLICATION_ROOT=\$QuickIsoSeq/Tools

STAR=\$APPLICATION_ROOT/STAR_$STAR_VERSION/bin/Linux_x86_64_static
FEATURECOUNTS=\$APPLICATION_ROOT/${SUBREAD_VERSION}/bin
VARSCAN_JAR=\$APPLICATION_ROOT/VarScan.v2.4.0.jar
SAMTOOLS=\$APPLICATION_ROOT/samtools-${SAMTOOLS_VERSION}/bin
RSEM=\$APPLICATION_ROOT/$RSEM_PKG
SALMON=\$APPLICATION_ROOT/${SALMON_PKG}/bin
KALLISTO=\$APPLICATION_ROOT/kallisto
GFFREAD=\$APPLICATION_ROOT/gffread-0.11.4.Linux_x86_64
BOWTIE=\$APPLICATION_ROOT/bowtie-$BOWTIE_VERSION-linux-x86_64
FASTQC=\$APPLICATION_ROOT/FastQC

export PATH=\$SAMTOOLS:\$STAR:\$FEATURECOUNTS:\$RSEM:\$SALMON:\$KALLISTO:\$GFFREAD:\$BOWTIE:\$FASTQC:\$PATH
" > $QuickIsoSeq/Tools/tools_path.txt

#
# Check whether multiqc installed
#
command_exists () {
    type "$1" &> /dev/null ;
}
if ! command_exists multiqc; then
    echo "multiqc (a python-based package) is recommeded but not installed"
	echo "Install MultiQC:  pip install multiqc"
fi
