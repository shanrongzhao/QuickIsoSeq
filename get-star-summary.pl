#!/usr/bin/perl
#
# Merge the mapping summary from the Log.final.out generated from STAR
#
# Usage
#	perl get-star-summary.pl allIDs.txt <summary file>[optional, default to star-mapping-summary.txt]
#
# Shanrong Zhao
#
# October 3th, 2018
#
# you have to run this script at the PARENT directory for STAR runs
#

use strict;

if (@ARGV < 1) {
	print "Usage: $0 <id file> <summary file>[optional, default to star-mapping-summary.txt]\n";
	exit;
}

my $Summary_File="star-mapping-summary.txt";
if (@ARGV > 1) {
	$Summary_File=$ARGV[1];
}

#fields to be extracted
	
# $ awk '{print "#" $0}' Log.final.out
#                                 Started job on |      Nov 21 22:22:21
#                             Started mapping on |      Nov 21 22:24:39
#                                    Finished on |      Nov 21 22:48:46
#       Mapping speed, Million of reads per hour |      103.87
#
#                          Number of input reads |      41751259
#                      Average input read length |      150
#                                    UNIQUE READS:
#                   Uniquely mapped reads number |      36897631
#                        Uniquely mapped reads % |      88.37%
#                          Average mapped length |      149.88
#                       Number of splices: Total |      17118933
#            Number of splices: Annotated (sjdb) |      16946761
#                       Number of splices: GT/AG |      16968171
#                       Number of splices: GC/AG |      115307
#                       Number of splices: AT/AC |      13661
#               Number of splices: Non-canonical |      21794
#                      Mismatch rate per base, % |      0.81%
#                         Deletion rate per base |      0.01%
#                        Deletion average length |      1.56
#                        Insertion rate per base |      0.01%
#                       Insertion average length |      1.32
#                             MULTI-MAPPING READS:
#        Number of reads mapped to multiple loci |      1923250
#             % of reads mapped to multiple loci |      4.61%
#        Number of reads mapped to too many loci |      20882
#             % of reads mapped to too many loci |      0.05%
#                                  UNMAPPED READS:
#       % of reads unmapped: too many mismatches |      0.00%
#                 % of reads unmapped: too short |      6.95%
#                     % of reads unmapped: other |      0.02%


sub getSummary {
	my $id = shift;
	open (SUMMARY, "<$id/Log.final.out") || die $@;

	my ($total_reads, $uniq_reads, $read_length, $map_length, $multi_rate, $uniq_rate, $unmap_rate);
	
	while (my $aline = <SUMMARY>) {
	#	next unless $aline =~ /reads/;
		chomp($aline);

		$aline =~ s/^\s+|\s+$//g;
		my @aitems = split(/\|/, $aline);

		my $key = $aitems[0];
		$key =~ s/\s+$//g;
		my $val= $aitems[1];
		$val =~ s/\s+//g;

		$total_reads=$val if ($key eq "Number of input reads");
		$uniq_reads=$val if ($key eq "Uniquely mapped reads number");
		$read_length=$val if ($key eq "Average input read length");
		$map_length=$val if ($key eq "Average mapped length");
		$uniq_rate=$val if ($key eq "Uniquely mapped reads %");
		$multi_rate=$val if ($key eq "% of reads mapped to multiple loci");
	}
	
	close(SUMMARY);
	
	$uniq_rate =~ s|(\d+\.\d+)%|$1|eg;
	$multi_rate =~ s|(\d+\.\d+)%|$1|eg;
	$unmap_rate = 100 - $uniq_rate - $multi_rate;
	$unmap_rate = sprintf("%5.2f", $unmap_rate);
	return "$id\t$total_reads\t$uniq_rate\t$multi_rate\t$unmap_rate";
}


open (OUTPUT, ">$Summary_File") || die $@;
print OUTPUT "Sample\tTotal_reads\tUniq_Rate\tMulti_Rate\tUnmap_Rate\n";

open(FILE, $ARGV[0]) || die $@;
while(my $line = <FILE>) {
	chomp($line);
	my @items = split(/\t/, $line);
	my $id = $items[0];

	print OUTPUT getSummary($id), "\n";
}

close(FILE);
close(OUTPUT);
