#!/usr/bin/perl -w
# filter the output from hmmsearch2tab.pl to keep only best hit for each gene (for ScEDS project)
#by Jester, 2018/2/26

use strict;

@ARGV==1 or die "Usage: $0 <hmmsearch2tab.pl ouput>\n";

my $file=shift;
my @result=();
my $format=0;
	open F, $file or die "Error open file $file";
	while(<F>)
	{
		chomp;
		next unless $_;
		my @item=split(/\t/);
		next unless @item>1;
		if($item[0] eq 'HMM') #check format
		{
			$format=$_;
		}
		else
		{
			next if $item[1]=~/no hits/;
			die "Error ID format for '$item[3]'" unless $item[3]=~/^(\S+)\(\d+[\<\>]\d+\)$/;
			my ($id)=($1);
			push @result,[$id,@item];
		}
	}
	close F;
die "Error format of input file" unless $format;

my %ids=();
print $format,"\n" if @result;
foreach my $hit (sort{$b->[2]<=>$a->[2]}@result) #order the results by HMM score decreasely
{
	my $id=shift @{$hit};
	unless(exists $ids{$id}) #if already output the best hit, skip
	{
		print join("\t",@{$hit}),"\n";
		$ids{$id}=1;
	}
}
