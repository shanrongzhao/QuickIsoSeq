#!/usr/bin/perl 
##
# Extract isoform annotation from a  GFT file
#
#
# Shanrong Zhao
#
# July 23, 2019
#

use strict;

if (@ARGV < 1) {
        print "Usage: $0 <GTF file> \n";
        exit;
}

open(GTF, "$ARGV[0]") || die "Can't open $ARGV[0].\n";

#gene
#1       ensembl gene    43401   63280   .       -       .       gene_id "ENSMFAG00000044637"; gene_version "1"; gene_name "PGBD2"; gene_source "ensembl"; gene_biotype "protein_coding";
#transcript
#1       ensembl transcript      43401   63280   .       -       .       gene_id "ENSMFAG00000044637"; gene_version "1"; transcript_id "ENSMFAT00000000550"; transcript_version "1"; gene_name "PGBD2"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "PGBD2-201"; transcript_source "ensembl"; transcript_biotype "protein_coding";

my %all_txes;
while(<GTF>){
	next if(/^##/); #ignore header
	chomp;
	  
	my ($chr, $source, $type, $start, $end, $score, 
		$strand, $phase, $attributes) = split("\t");
		
	next unless ($type eq "transcript" || $type eq "exon");   #only process transcript/exon rows
	
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
	# if (! exists $attribs{'gene_type'}) {
		 # $attribs{'gene_type'} =  $attribs{'gene_biotype'}
	# }

    if (! exists $attribs{'gene_name'}) {
        $attribs{'gene_name'} =  $attribs{'gene_id'}
    }
	if (! exists $attribs{'transcript_name'}) {
        $attribs{'transcript_name'} =  $attribs{'transcript_id'}
    }
	if (! exists $attribs{'transcript_biotype'}) {
        $attribs{'transcript_biotype'} =  $attribs{'transcript_type'}
    }

	#transcript
	#1       ensembl transcript      43401   63280   .       -       .       
	#gene_id "ENSMFAG00000044637"; gene_version "1"; 
	#transcript_id "ENSMFAT00000000550"; transcript_version "1"; 
	#gene_name "PGBD2"; gene_source "ensembl"; gene_biotype "protein_coding"; 
	#transcript_name "PGBD2-201"; transcript_source "ensembl"; transcript_biotype "protein_coding";

	if ($type eq "transcript" ) {
		$all_txes{$attribs{'transcript_id'}} = {
			transcript_ID => $attribs{'transcript_id'},
			gene_ID  => $attribs{'gene_id'},
			transcript_name => $attribs{'transcript_name'},
			gene_name =>  $attribs{'gene_name'},
			transcript_type => $attribs{'transcript_biotype'}
		};
	} 

	if ($type eq "exon" ){
		if (! exists $all_txes{$attribs{'transcript_id'}}{'transcript_length'} ) {
			$all_txes{$attribs{'transcript_id'}}{'transcript_length'} = $end - $start + 1;
		} else {
			$all_txes{$attribs{'transcript_id'}}{'transcript_length'} += $end - $start + 1;
		}
	}

}

print "transcript_ID\tgene_ID\ttranscript_name\tgene_name\ttranscript_length\ttranscript_type\n";
foreach my $tx (sort keys %all_txes) {
	my %r = %{$all_txes{$tx}};
	print "$r{'transcript_ID'}\t$r{'gene_ID'}\t$r{'transcript_name'}\t$r{'gene_name'}\t$r{'transcript_length'}\t$r{'transcript_type'}\n";
}
