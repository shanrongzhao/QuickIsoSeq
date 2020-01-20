#!/usr/bin/perl
#
# Split the total number of sequenced reads into
#		1. rRNA reads  
#		2. hgRNA reads (hemoglobin, check the effiency of globin clear)
#		3. other
#
# Usage
#	get-reads-split.pl allIDs.txt reads-split-summary.txt[optional]
#
#
# Shanrong Zhao
#
# December 24, 2018
#
use strict;

if (@ARGV < 1) {
	#category: rRNA, or hgRNA
	print "Usage: $0 <id file> <summary file>[optional, default to reads-split-summary.txt]";
	exit;
}

my $Summary_File="reads-split-summary.txt";
if (@ARGV > 1) {
	$Summary_File=$ARGV[1];
}

sub getSummary {
	my $bwt_file = shift;
	open (SUMMARY, "<$bwt_file") || die $@;
	# reads processed: 7147623
	# reads with at least one reported alignment: 88415 (1.24%)
	# reads that failed to align: 7059208 (98.76%)
	# Reported 88415 alignments to 1 output stream(s)

	my ($total_reads, $RNA_reads, $rate );
	
	while (my $aline = <SUMMARY>) {
	#	next unless $aline =~ /reads/;
		chomp($aline);

		$aline =~ s/^\s+|\s+$//g;
		my @aitems = split(/\:/, $aline);

		my $key = $aitems[0];
		$key =~ s/\s+$//g;
		my $val= $aitems[1];
		$val =~ s/\s+//g;

		$total_reads=$val if ($key eq "# reads processed");
		$RNA_reads=$val if ($key eq "# reads with at least one reported alignment");
	}
	
	close(SUMMARY);
	$RNA_reads =~ s/\(/ /;
	$RNA_reads =~ s/%\)//;
	my @RNAs = split(/ /, $RNA_reads);

	my @results = ($total_reads, $RNAs[1]);
	return @results;
}


open (OUTPUT, ">$Summary_File") || die $@;
print OUTPUT "Sample\tTotal_seq\trRNA\thgRNA\tOther\n";

open(FILE, $ARGV[0]) || die $@;
while(my $line = <FILE>) {
	chomp($line);
	my @items = split(/\t/, $line);
	my $id = $items[0];

	my @rRNAs = getSummary("./${id}/${id}.bwt.rRNA.log");
	my @hgRNAs = getSummary("./${id}/${id}.bwt.hgRNA.log");
	my $other = 100- $rRNAs[1] - $hgRNAs[1];
	print OUTPUT "$id\t$rRNAs[0]\t$rRNAs[1]\t$hgRNAs[1]\t$other\n";
}
close(FILE);
close(OUTPUT);
