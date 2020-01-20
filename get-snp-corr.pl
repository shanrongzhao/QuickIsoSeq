#!/usr/bin/perl
#
# SNP concordance
#
# Usage
#	perl get-snp-corr.pl allIDs.txt <summary file>[optional, default to RNASeq.SNP.corr.txt] <Method> [default to concordant]
#
#
# Shanrong Zhao, Baohong Zhang
#
# October 12, 2018
#
use strict;

if (@ARGV < 1) {
	print "Usage: $0 <id file> <correlation file>[optional, default to RNASeq.snp.corr.txt] <Method> [optional,default to concordant]\n";
	print "Method: concordant or total\n";
	exit;
}

my $SNP_Corr_File="RNASeq.snp.corr.txt";
if (@ARGV > 1) {
	$SNP_Corr_File=$ARGV[1];
}

my $Method="concordant";
if (@ARGV > 2) {
	$Method=lc($ARGV[2]);
}


# /hpc/grid/ngsws/molmed/zhaos25/SNP_Coor$ head US-1562856/US-1562856.snp.txt
# Chrom   Position        Ref     Cons    Reads1  Reads2  VarFreq Strands1        Strands2        Qual1   Qual2   Pvalue  MapQual1        MapQual2     Reads1Plus       Reads1Minus     Reads2Plus      Reads2Minus     VarAllele
# chr6    292447  T       Y       44      5       10.2%   2       1       55      30      0.98    1       1       22      22      0       5       C
# chr6    292512  G       K       102     13      11.3%   2       2       66      64      0.98    1       1       72      30      10      3       T
# chr6    304644  A       R       183     4       2.14%   2       2       66      66      0.98    1       1       105     78      2       2       G
# chr6    311938  C       Y       134     4       2.9%    2       2       68      69      0.98    1       1       63      71      2       2       T

#
#get ID list and read SNPs
#
my @ids =();
my %SNPS = ();
my %SNPS2 = ();  #no requirements on coverage and allele abundance

open(FILE, $ARGV[0]) || die $@;

while(my $line = <FILE>) {
	chomp($line);
	my @items = split(/\t/, $line);
	push @ids, $items[0];
	my $id = $items[0];
	
	my %refSNP = ();
	my %refSNP2 = ();
	
	open(FILE2, "$id/$id.snp.txt") || die $@;
	<FLIE2>;
	
	while ($line = <FILE2>) {
        $line =~ s/\s+$//g;
        my @items = split(/\t/, $line);
        $items[6] =~ s/%//g;
		
		$refSNP2{$items[1]} = $items[3];
		
		# at least 10 reads for each allele and at least 25% alternative allele
        next if ($items[4] < 10 || $items[5] < 10 || $items[6] < 25); 

        $refSNP{$items[1]} = $items[3];
	}
	
	close(FILE2);
	
	$SNPS{$id} = \%refSNP;
	$SNPS2{$id} = \%refSNP2;
}

close(FILE);


#
# Calculate correlation
#
my $no = @ids;
my %scores = ();

for (my $i=0; $i < $no; $i++) {
	my $id1 = $ids[$i];
	my %SNP1 = %{$SNPS{$id1}};
	my $total_refSNP = scalar(keys %SNP1);

	
	for (my $j=0; $j < $no; $j++) {
		next if ($i == $j);
		
		my $id2 = $ids[$j];
		#my %SNP2 = %{$SNPS{$id2}};
		my %SNP2 = %{$SNPS2{$id2}};
		
		my $m = 0;
		my $common = 0;
		
		while ( my ($key, $value) = each %SNP2 ) {
			if (defined $SNP1{$key}) {
				$common ++;
				
				if ($SNP1{$key} eq $value) {
					$m++;
				}
			}
		}
		
		my $score = sprintf("%.2f",$m*100/($total_refSNP+0.01));
		if ( $Method eq "concordant") {
			$score = sprintf("%.2f",$m*100/($common + 0.01));   #"concordant" is the default
		} 
		
		$scores{$id1 . "|" . $id2} = $score;
		#$scores{$id2 . "|" . $id1} = $score;
	}
}


#
# Output correlation
#

open (OUTPUT, ">$SNP_Corr_File") || die $@;

my $ids_all = join "\t", @ids;
print OUTPUT "Pair\t$ids_all\n";

foreach my $id1 (@ids) {
	
	#read SNP in sample B
	my $result = "$id1";

	foreach my $id2 (@ids) {
		if ( $id1 eq $id2) {
			$result = $result . "\t" . "100";
			next;
		} 
		
		$result = $result . "\t" . $scores{ $id1 . "|" . $id2};
	}

	print OUTPUT "$result","\n";
}

close(OUTPUT);

