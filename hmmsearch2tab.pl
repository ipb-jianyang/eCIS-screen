#!/usr/bin/perl -w
#use bioperl to parse hmmsearch result to tabular report for importing into Excel.
#by Jester, 2017/12/27

use strict;
use Bio::SearchIO;

@ARGV==1 or die "Usage: $0 <hmmpfam ouput>\n";
-r $ARGV[0] or die "Error reading from $ARGV[0]!\n";

my $in = new Bio::SearchIO(-format => 'hmmer',-file => $ARGV[0]);
my $sw=0;
  while ( my $res = $in->next_result ){
	  $sw++;
	  print join("\t",qw/HMM Score E-value SeqID Desc/),"\n" if($sw==1); #output title line only once
    # get a Bio::Search::Result::HMMERResult object
	my $flag=0;
    my $hmm=$res->query_name;
	next unless $hmm;
#	die "Error: no query HMM name found" unless $hmm;
    while ( my $hit = $res->next_hit ){
		$flag++;
        print join("\t",$hmm,$hit->score,$hit->significance,$hit->name,$hit->description),"\n"
    }
	print "$hmm\t[no hits above thresholds]\n" unless($flag);
  }
print STDERR "Nothing found from your input, maybe an error input format?\n" unless($sw);
