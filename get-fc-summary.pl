#!/usr/bin/perl
#
# Merge the counting summay from the *featureCounts.txt.summary generated from featureCounts
#
# Usage
#	perl get-fc-summary.pl allIDs.txt <summary file>[optional, default to fc-counting-summary.txt]
#
# Shanrong Zhao
#
# October 8, 2018
#
# you need to run this script at the PARENT directory for featureCounts runs

use strict;

if (@ARGV < 1) {
	print "Usage: $0 <id file> <summary file>[optional, default to fc-counting-summary.txt]\n";
	exit;
}

my $Summary_File="fc-counting-summary.txt";
if (@ARGV > 1) {
	$Summary_File=$ARGV[1];
}


################################################################
# get the counting summary
#
# cat US-1495821.featureCounts.txt.summary
# Status  US-1495821.sort.bam
# Assigned        45707059
# Unassigned_Ambiguity    1390102
# Unassigned_MultiMapping 6394295
# Unassigned_NoFeatures   6524502
# Unassigned_Unmapped     0
# Unassigned_MappingQuality       0
# Unassigned_FragementLength      0
# Unassigned_Chimera      0
# Unassigned_Secondary    0
# Unassigned_Nonjunction  0
# Unassigned_Duplicate    0

sub getCountingSummary {
	my $id = shift;
	open (SUMMARY, "<$id/$id.featureCounts.txt.summary") || die $@;

	my ($uniq_reads, $gene_reads, $ambiguity_reads, $no_feature_reads,$gene_rate, $ambiguity_rate, $no_feature_rate);
	while (my $aline = <SUMMARY>) {
		chomp($aline);

		$aline =~ s/^\s+|\s+$//g;
		my @aitems = split(/\t/, $aline);

		my $key = $aitems[0];
		my $val = $aitems[1];
		
		$gene_reads = $val if ($key eq "Assigned");
		$ambiguity_reads = $val if ($key eq "Unassigned_Ambiguity");
		$no_feature_reads = $val if ($key eq "Unassigned_NoFeatures");
	}
	
	close(SUMMARY);
	
	$uniq_reads = $gene_reads + $ambiguity_reads + $no_feature_reads;
	$gene_rate = sprintf("%5.2f",$gene_reads *100/$uniq_reads);
	$ambiguity_rate = sprintf("%5.2f",$ambiguity_reads *100/$uniq_reads);
	$no_feature_rate = sprintf("%5.2f",$no_feature_reads *100/$uniq_reads);
	
	return "$id\t$uniq_reads\t$gene_rate\t$ambiguity_rate\t$no_feature_rate";
}


open (OUTPUT, ">$Summary_File") || die $@;
print OUTPUT "Sample\tUniq_Mapped_Reads\tGene_Rate\tAmbiguity_Rate\tNo_Feature_Rate\n";

open(FILE, $ARGV[0]) || die $@;
while(my $line = <FILE>) {
	chomp($line);
	my @items = split(/\t/, $line);
	my $id = $items[0];

	print OUTPUT getCountingSummary($id), "\n";
}

close(FILE);
close(OUTPUT);