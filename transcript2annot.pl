#!/usr/bin/perl
##
# polish Gencode transcript fasta files
#
# Shanrong Zhao
#
# September 29, 2018
#

use strict;

if (@ARGV < 3) {
        print "Usage: $0 <transcript FASTA file>[Input]  <polished transcript FASTA file>[Output] <transcript FASTA annotation file>[Output] \n";
        exit;
}

open(FASTA, "$ARGV[0]") || die "Can't open $ARGV[0].\n";

open POLISH, ">$ARGV[1]";
open ANNOT, ">$ARGV[2]";

print ANNOT "transcript_ID\tgene_ID\ttranscript_name\tgene_name\ttranscript_length\ttranscript_type\n";
foreach my $line (<FASTA>){
	if ($line =~ /^>/){
		chomp $line;
		my @temp = split(/\|/, $line);
		print POLISH $temp[0], "\n";
		$temp[0] =~ s/^>//;
		print ANNOT "$temp[0]\t$temp[1]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\n";
	}
	else{
		print POLISH $line;
	}
}
