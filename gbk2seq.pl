#!/usr/bin/perl -w
# simple extract the CDS sequences (DNA or protein) from a genbank file, with it's annotation
#Made by Jester, 2004/7/29
#revised at 2004/11/22, jester

use strict;
use Bio::SeqIO;
my $usage = "Usage $0 genbank_file <accession> <FFN(DNA)|FAA(protein)|PG(pseudogene)> [options]\nAvailable options:\n\t-s\tsimple model, i.e. synonym only\n";
my $filename = shift || die $usage;
my $acc= shift || die $usage;
my $type = shift || die $usage;
$type = uc $type;
die $usage unless($type eq 'FFN' or $type eq 'FAA' or $type eq 'PG');
my $simple=0;
if(@ARGV)
{
	$simple=shift;
	die "Error: unknown option '$simple'!\n$usage" unless((uc $simple) eq '-S');
	$simple=1;
}
warn "OK, extract ",$type eq 'FAA' ? 'protein' : ($type eq 'FFN' ? 'gene' : 'pseudogene')," sequences", $simple ? ' in simplified mode': '', " now ...\n";
my $in = new Bio::SeqIO(-file => $filename, -format => 'genbank');

my %geneIDs=(); #check potential duplicate IDs
while( $acc and my $seq = $in->next_seq ) {
	next unless $acc eq $seq->accession;
	$acc='';
#revised for adopting the new policy of NCBI (i.e. no CDS for pseudogenes)
  my @cds = $type eq 'PG' ? grep { $_->primary_tag eq 'gene' and exists $_->{'_gsf_tag_hash'}->{'pseudo'} } $seq->get_SeqFeatures() : grep { $_->primary_tag eq 'CDS' or ($_->primary_tag eq 'gene' and exists $_->{'_gsf_tag_hash'}->{'pseudo'}) } $seq->get_SeqFeatures();
  my $cdsnum=@cds;
  if($cdsnum)
	{
	print STDERR "Total $cdsnum ",$type eq 'PG' ? 'pseudogenes':'CDS'," found in ",$seq->display_id(),".\nStarting extract sequences now ...\n";
	}
	else
	{
		print STDERR "[WARN] No ",$type eq 'PG' ? 'pseudogenes':'CDS'," found in ",$seq->display_id(),". Skipped.\n";
	}
  my $counter=0;
  foreach my $cds ( @cds ) {
	my $pro=exists $cds->{'_gsf_tag_hash'}->{'translation'}->[0]?$cds->{'_gsf_tag_hash'}->{'translation'}->[0]:'';
	next if($type eq 'FAA' and (not $pro));
	next if($type eq 'PG' and (not (exists $cds->{'_gsf_tag_hash'}->{'pseudo'})));
	$pro=~s/(\w{60})/$1\n/g;
    my $featureseq = $cds->spliced_seq;
	my $cdsseq=$featureseq->seq();
	$cdsseq=~s/(\w{60})/$1\n/g;
	print '>';
	if($simple)
	  {
		die "Error: no 'locus_tag' found for CDS at $cds->{'_location'}->{'_start'}\.\.$cds->{'_location'}->{'_end'}!\n" unless(exists $cds->{'_gsf_tag_hash'}->{'locus_tag'}->[0]);
		print $cds->{'_gsf_tag_hash'}->{'locus_tag'}->[0],' ',exists $cds->{'_gsf_tag_hash'}->{'gene'}->[0]?$cds->{'_gsf_tag_hash'}->{'gene'}->[0]:'-';
	  }
	  else
	  {
#find ID for each gene
		my $id='';
		if(exists $cds->{'_gsf_tag_hash'}->{'locus_tag'})
		{ $id=$cds->{'_gsf_tag_hash'}->{'locus_tag'}->[0]; }
		else
		{
			if($cds->has_tag('gene'))
			{
				($id)=$cds->get_tag_values('gene');
				$id='' unless $id=~/^[\w\_]+\d{4,6}$/; #guess it as ID if it looks like an locus_tag
				die "Error: duplicate IDs found '$id'" if exists $geneIDs{$id}; #check for validate
			}
			unless($id)
			{
			 if($cds->has_tag('protein_id'))
			 {($id)=$cds->get_tag_values('protein_id'); }
			}
		}
		die "Fatal error: no ID found" unless $id;
		$geneIDs{$id}=$id;
		print $id,'(',$cds->start,$cds->strand>0?'>':'<',$cds->end,') ', exists $cds->{'_gsf_tag_hash'}->{'product'}->[0]?$cds->{'_gsf_tag_hash'}->{'product'}->[0]:'-';
	  }
	print "\n", $type eq 'FAA' ? $pro : $cdsseq,"\n";
  }
}
print STDERR "Done.\n";
