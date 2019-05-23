#!/usr/bin/perl -w
# parse the output from hmmsearch2tab.pl for eCIS
#by Jester, 2017/12/28
# revised to use essential & core components for inclusion cutoff, Jester 2018/2/24

use strict;
my $dis_cutoff=30000; #arbitrary setting for AFP system size
my $essential='afp11';
#my @core=qw/afp1_5 afp2_3_4 afp8 afp9 afp15/;
my @core=qw/afp1_5/; #revised to use afp11 and afp1_5 only, Jester 2018/6/9
my $num_genes=1; #lower limit for number of core genes in the AFP system

@ARGV>0 or die "Usage: $0 <hmmsearch2tab.pl ouput1> [output2] ...\n\tDistance cutoff: $dis_cutoff\n\tGene number cutoff: $num_genes\n";

my @afp=qw/afp1_5 afp2_3_4 afp6 afp7 afp8 afp9 afp10 afp11 afp12 afp13 afp14 afp15 afp16/; #predefined
my %result=();

warn "Reading input now ...\n";
foreach my $file (@ARGV)
{
	open F, $file or die "Error open file $file";
	$file=~s/.*\///g; #remove path
	$file=~s/\.hmmout\.txt$//; #use as genome name
	warn "$file ...\n";
	my $format=0;
	my @region=();
	while(<F>)
	{
		chomp;
		next unless $_;
		my @item=split(/\t/);
		next unless @item>1;
		if($item[0] eq 'HMM') #check format
		{
			$format=1;
			$result{$file}=[];
		}
		else
		{
			next if $item[1]=~/no hits/;
			die "Error ID format for '$item[3]'" unless $item[3]=~/^(\S+)\((\d+)[\<\>](\d+)\)$/;
			my ($id,$from,$to)=($1,$2,$3);
			if(@region)
			{
				my @close=(); #keep the closest region id and distance
				foreach my $i (0..$#region) #find the closest region
				{
					my $r=$region[$i];
					if($r->[0]>$to)
					{
						@close=($r->[0]-$to,$i) unless @close and $close[0]<$r->[0]-$to;
					}
					elsif($r->[1]<$from)
					{
						@close=($from-$r->[1],$i) unless @close and $close[0]<$from-$r->[1];
					}
					else #the region cover/overlap the gene
					{
						@close=(0,$i);
						last;
					}
				}
				if($close[0]<=$dis_cutoff) #the current gene can be assigned to exist region
				{
					push @{$result{$file}->[$close[1]]->{$item[0]}},$id;
					@{$region[$close[1]]}=(&min($from,$region[$close[1]]->[0]),&max($to,$region[$close[1]]->[1])); #extend the region if necessary
				}
				else #produce new region
				{
					push @region,[$from,$to];
					push @{$result{$file}->[$#region]->{$item[0]}},$id;
				}
			}
			else
			{
				push @region,[$from,$to];
				push @{$result{$file}->[0]->{$item[0]}},$id;
			}
		}
	}
	close F;
#	warn "\t",scalar @region," regions found.\n"; #for debug
	die "Error format of file $file" unless $format;
}

#output result for each input
warn "Create output now ...\n";
print join("\t",'#Genome',@afp),"\n"; #title line
foreach my $g (sort(keys %result))
{
	foreach my $reg (sort{scalar keys %{$b} <=> scalar keys %{$a}} @{$result{$g}})
	{
		last unless exists $reg->{$essential}; #essential gene must exist!
		my $num_core=0;
		map { $num_core++ if exists $reg->{$_}} @core;
		last if $num_core<$num_genes; #region ignored if no enough core genes found
		print $g;
		foreach my $id (@afp)
		{
			print "\t",exists $reg->{$id} ? join(",",@{$reg->{$id}}) : '-';
		}
		print "\n";
	}
}
warn "Done.\n";

sub max{
my @a=@_;
@a=sort{$b<=>$a}@a;
return $a[0];
}

sub min{
my @a=@_;
@a=sort{$a<=>$b}@a;
return $a[0];
}
