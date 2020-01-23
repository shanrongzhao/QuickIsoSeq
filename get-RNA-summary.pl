#!/usr/bin/perl
#
# Merge the RNA summary 
#
# Usage
#	get-RNA-summary.pl allIDs.txt Category  ${Category}-summary.txt
#
# The %RNA was calculated bu check-RNA.sh, the output look like
#
# -bash-4.2$ cat SRR607214/SRR607214.bwt.hgRNA.log
# reads processed: 39769361
# reads with at least one reported alignment: 15939527 (40.08%)
# reads that failed to align: 23829834 (59.92%)
#
# Shanrong Zhao
#
# December 24, 2018
#
use strict;


if (@ARGV < 2) {
	#category: rRNA, or hgRNA
	print "Usage: $0 <id file> <category> <summary file>[optional, default to Category-summary.txt]";
	exit;
}

my $CATEGORY = $ARGV[1];
my $Summary_File="${CATEGORY}-summary.txt";
if (@ARGV > 2) {
	$Summary_File=$ARGV[2];
}

sub getSummary {
	my $id = shift;
	open (SUMMARY, "<./${id}/${id}.bwt.${CATEGORY}.log") || die $@;
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

	return "$id\t$total_reads\t$RNAs[0]\t$RNAs[1]";
}


open (OUTPUT, ">$Summary_File") || die $@;
print OUTPUT "Sample\tTotal_reads\t${CATEGORY}_Reads\t${CATEGORY}_Rate\n";

open(FILE, $ARGV[0]) || die $@;
while(my $line = <FILE>) {
	chomp($line);
	my @items = split(/\t/, $line);
	my $id = $items[0];

	print OUTPUT getSummary($id), "\n";
}
close(FILE);
close(OUTPUT);
