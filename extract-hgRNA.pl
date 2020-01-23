#!/usr/bin/perl
##
# Extract hemoglobin sequences form gencode transcript fasta files
#
# Shanrong Zhao
#
# December 24, 2018
#

use strict;

if (@ARGV < 2) {
        print "Usage: $0 <transcript FASTA file>[Input]  <hgRNA FASTA file>[Output]\n";
        exit;
}

open(FASTA, "$ARGV[0]") || die "Can't open $ARGV[0].\n";

open RNA, ">$ARGV[1]";

#
# hemoglobin RNA sequences from Ensembl (r91) were utilized, specifically 
# ten protein-coding hemoglobin subunit genes and two pseudogenes
#
my @hgGenes = qw( HBA1 HBA2 HBB HBBP1 HBD HBE1 HBG1 HBG2 HBM HBQ1 HBZ HBZP1);
my %hgGenes_hash = map { $_ => 1 } @hgGenes;

my $hgRNA = "no";
foreach my $line (<FASTA>){
	
	if ($line =~ /^>/){
		chomp $line;
		my @temp = split(/\|/, $line);
		
		$hgRNA = "no";
		##ENST00000606857.1|ENSG00000268020.3|OTTHUMG00000185779.1|OTTHUMT00000471235.1|OR4G4P-201|OR4G4P|840|unprocessed_pseudogene|
		if( exists( $hgGenes_hash{$temp[5]} ) )  {
			$hgRNA = "yes";
			print RNA "$temp[0]\t$temp[4]\t$temp[5]\t$temp[6]\n";
		}
	}
	else{
		print RNA $line if ( $hgRNA eq "yes" );
	}
}

#
# A total of 37 transcript in Gencode Release version 29
#
# extract-hgRNA.pl gencode.v29.transcripts.fa hgRNAs.fasta
#
# grep ">" hgRNAs.fasta | sort -k2,2
# ID					Transcipt_Name	Gene_Name	Trancript_length
#ENST00000320868.9      HBA1-201        HBA1    577
#ENST00000397797.1      HBA1-202        HBA1    504
#ENST00000472694.1      HBA1-203        HBA1    674
#ENST00000487791.1      HBA1-204        HBA1    410
#ENST00000251595.11     HBA2-201        HBA2    576
#ENST00000397806.1      HBA2-202        HBA2    513
#ENST00000482565.1      HBA2-203        HBA2    640
#ENST00000484216.1      HBA2-204        HBA2    403
#ENST00000335295.4      HBB-201 HBB     628
#ENST00000380315.2      HBB-202 HBB     502
#ENST00000475226.1      HBB-203 HBB     319
#ENST00000485743.1      HBB-204 HBB     680
#ENST00000633227.1      HBB-205 HBB     609
#ENST00000647020.1      HBB-206 HBB     754
#ENST00000433329.1      HBBP1-201       HBBP1   439
#ENST00000454892.2      HBBP1-202       HBBP1   1450
#ENST00000292901.7      HBD-201 HBD     486
#ENST00000380299.3      HBD-202 HBD     785
#ENST00000417377.1      HBD-203 HBD     305
#ENST00000429817.1      HBD-204 HBD     566
#ENST00000643122.1      HBD-205 HBD     806
#ENST00000650601.1      HBD-206 HBD     620
#ENST00000292896.3      HBE1-201        HBE1    918
#ENST00000380237.5      HBE1-202        HBE1    913
#ENST00000396895.2      HBE1-203        HBE1    616
#ENST00000330597.4      HBG1-201        HBG1    779
#ENST00000632727.1      HBG1-202        HBG1    369
#ENST00000648735.1      HBG1-203        HBG1    1464
#ENST00000336906.5      HBG2-201        HBG2    614
#ENST00000380252.6      HBG2-202        HBG2    694
#ENST00000444587.1      HBG2-203        HBG2    307
#ENST00000356815.3      HBM-201 HBM     506
#ENST00000472539.5      HBM-202 HBM     592
#ENST00000496585.1      HBM-203 HBM     592
#ENST00000199708.2      HBQ1-201        HBQ1    536
#ENST00000252951.2      HBZ-201 HBZ     755
#ENST00000354915.3      HBZP1-201       HBZP1   429
