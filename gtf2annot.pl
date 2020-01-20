#!/usr/bin/perl 
##
# Extract annotation from a  GFT file
#
# Shanrong Zhao
#
# October 6, 2018
#

use strict;

if (@ARGV < 1) {
        print "Usage: $0 <GTF file> \n";
        exit;
}

open(GTF, "$ARGV[0]") || die "Can't open $ARGV[0].\n";

my %all_genes;
while(<GTF>){
	next if(/^##/); #ignore header
	chomp;
	  
	my ($chr, $source, $type, $start, $end, $score, 
		$strand, $phase, $attributes) = split("\t");
		
	next unless $type eq "gene";   #only process gene rows
	
	my %attribs = ();
	my @add_attributes = split(";", $attributes);
	
	foreach my $attr ( @add_attributes ) {
		next unless $attr =~ /^\s*(.+)\s\"(.+)\"$/;
		my $c_type  = $1;
		my $c_val = $2;
		if($c_type  && $c_val){
			$attribs{$c_type} = $c_val;
		}
	}

	#correction for rat gtf download from Ensembl, in which gene_biotype is used instead of gene_type
	if (! exists $attribs{'gene_type'}) {
		 $attribs{'gene_type'} =  $attribs{'gene_biotype'}
	}

        if (! exists $attribs{'gene_name'}) {
                 $attribs{'gene_name'} =  $attribs{'gene_id'}
        }

	$all_genes{$attribs{'gene_id'}} = {
		gene_name  => $attribs{'gene_name'},
		gene_type  => $attribs{'gene_type'},
		chr        => $chr,
		start      => $start,
		end        => $end,
		strand 	   => $strand
	};

}

print "gene\tgene_name\tgene_type\tchr\tstart\tend\tstrand\n";
foreach my $gene (sort keys %all_genes) {
	my %r = %{$all_genes{$gene}};
	print "$gene\t$r{'gene_name'}\t$r{'gene_type'}\t$r{'chr'}\t$r{'start'}\t$r{'end'}\t$r{'strand'}\n";
}
