#!/bin/perl
use strict;
use Getopt::Long;

my $QV = "^";

my $USAGE = "fasta2fastq.pl [options] reads.fa > reads.fq";
my $help;

my $res = GetOptions("help"      => \$help,
                     "qv=n"      => \$QV);
 
if ($help || !$res)
{
  print $USAGE;
  print "\n";
  print "Convert fasta sequences to fastq, assigning fake quality values for the bases\n";
  print "Options\n";
  print "  -qv <q> : Assign this as the fake quality values (default: \'$QV\')\n";
  exit 0;
}

my $header = undef;
my $seq    = undef;

sub printFastq
{
  return if !defined $header;

  print "\@$header\n";
  print "$seq\n";
  print "+\n";
  print $QV x length($seq);
  print "\n";
}

while (<>)
{
  chomp;

  if (/>(\S+)/)
  {
    printFastq();

    $header = $1;
    $seq = "";
  }
  else
  {
    $seq .= $_;
  }
}

printFastq();
