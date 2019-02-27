#!/usr/bin/perl -w
# simple extract the sequences IDs from a multiple seq genbank file
#Jester, 2018/6/14

use strict;
use Bio::SeqIO;
my $usage = "Usage $0 genbank_file\n";
my $filename = shift || die $usage;
my $in = new Bio::SeqIO(-file => $filename, -format => 'genbank');

while( my $seq = $in->next_seq ) {
  print $seq->accession,"\n";
}
