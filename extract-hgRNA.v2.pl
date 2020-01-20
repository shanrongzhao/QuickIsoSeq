#!/usr/bin/perl
##
# Extract hemoglobin sequences from gencode transcript fasta files
# with reference transcript annotation file
#
# Shanrong Zhao
#
# December 24, 2019
#

use strict;

if (@ARGV < 3) {
        print "Usage: $0 <transcript FASTA file>[Input] <transcript annotation file>[Input] <hgRNA FASTA file>[Output]\n";
        exit;
}

open(FASTA, "$ARGV[0]") || die "Can't open $ARGV[0].\n";
open(ANNOT, "$ARGV[1]") || die "Can't open $ARGV[1].\n";
open RNA, ">$ARGV[2]";

#
# hemoglobin RNA sequences from Ensembl (r91) were utilized, specifically 
# ten protein-coding hemoglobin subunit genes and two pseudogenes
# track transcript IDs corresponding to hemoglobin
#
my @hgGenes = qw( HBA1 HBA2 HBB HBBP1 HBD HBE1 HBG1 HBG2 HBM HBQ1 HBZ HBZP1);
my %hgGenes_hash = map { $_ => 1 } @hgGenes;

my @hgIDs = ();
foreach my $line (<ANNOT>){
	chomp $line;
	my @temp = split(/\t/, $line);
	my $ID = ">" . $temp[0];
	push @hgIDs, $ID if exists $hgGenes_hash{$temp[3]};
	#print "$ID\n" if exists $hgGenes_hash{$temp[3]};
}

my %hgGenes_hash = map { $_ => 1 } @hgIDs;

my $hgRNA = "no";
foreach my $line (<FASTA>){
	
	if ($line =~ /^>/){
		chomp $line;
		my @temp = split(/ /, $line);
		
		$hgRNA = "no";
		if( exists( $hgGenes_hash{$temp[0]} ) )  {
			$hgRNA = "yes";
			print RNA "$temp[0]\n";
		}
	}
	else{
		print RNA $line if ( $hgRNA eq "yes" );
	}
}

