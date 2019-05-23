#!/usr/bin/perl -w
#run hmmsearch of given HMM file against ALL local bacterial genomes
#desigend for eCIS analysis, jester, 2018/1/10
#revised to use newly downloaded GenBank complete genome data, jester, 2018/6/14

use strict;

#settings:
my $hmmsearch='/usr/local/bin/hmmsearch'; #the program from HMMER3 package
my $script_path='/your_path/eCIS-screen'; #location of all eCIS screen scripts

my $gbk2IDs=$script_path.'/gbk2IDs.pl';
my $gbk2seq=$script_path.'/gbk2seq.pl';
my $hs2tab=$script_path.'/hmmsearch2tab.pl';
my $filter=$script_path.'/filter_hmmtab.pl';
my $summary=$script_path.'/parse_hmmtab4eCIS.pl';

-x $gbk2IDs and -x $gbk2seq and -x $hmmsearch and -x $hs2tab and -x $filter and -x $summary or die "error run $gbk2seq or $hs2tab or $filter or $hmmsearch";

@ARGV==3 or die <<EOF;
	<-- Pipeline for screening eCIS loci from bacterial genomes -->
	<--  By Jian YANG, Institute of Pathogen Biology, CAMS&PUMC -->
	<--  Available https://github.com/ipb-jianyang/eCIS-screen  -->

Usage: $0 <query HMM profile> <path to genomes for screen> <output summary file>

[NOTE] 1. The HMM profile file contains multiple HMMs for all components of eCISs.
       2. The path to genomes is expected to contain sub-directories for each genome,
          which should contain the GenBank format file named as *_genomic.gbff.gz.
          An example of the sub-directory within the path is available from:
          ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/196/155/GCA_000196155.1_ASM19615v1
       3. The pipeline will produce intermdiate files (e.g. CDS proteins, HMMER output)
          in the CURRENT working directory (where you run the script), these files can
          be helpful to improve the speed of the pipeline for further runs. But they
          are generally meaningless to users. The suggest way to run the pipeline is:
          # mkdir screen_output
          # cd screen_output
          # /path/HMMsearch_genomesII.pl /path/eCIS.hmm genome_path ../screen_summary.txt
       4. The pipeline may run for a few hours or days depends on the number of genomes
          for screening. The progress message will be printed to the screen by default.
       5. The output is in Tab-delimited text format, view in in Excel for easy read.
EOF

my $query=shift;
-e $query or die "error query hmm profile $query";
my $genomepath=shift;
-d $genomepath or die "error genome path given: $genomepath";
my $output=shift;
die "error: $output already exists" if -e $output;

my @allout=();
my @files=glob "$genomepath/*/*_genomic.gbff.gz"; #the file is gzipped by default
my $total_num=@files;
my $counter=0;
die "Sorry, no genomes found in '$genomepath'" unless @files;
warn 'Total ',scalar @files," genomes found.\n";
warn "Analysis HMM for each genome now ...\n";
foreach my $genome (@files)
{
 my $id=$genome;
 $id=~s!$genomepath!!;
 $id=~s!^/?([\w\.]+)/.*!$1!;
#unzip file first
 my $outgbk="$id.gbk";
 die "Error: $outgbk already exists" if -e $outgbk;
 !system("zcat $genome > $outgbk") or die "unzip file $genome with error";
#check the completeness of the gbk file
 my $tail=`tail -n 1 $outgbk`;
 chomp $tail;
 die "Error: $outgbk is like incomplete" unless $tail eq '//'; #check the last line
#get all replicon IDs from the gbk file
 my @repIDs=`$gbk2IDs $outgbk`;
 my $total_rep=scalar @repIDs;
 die "$gbk2IDs run error" unless $total_rep>0;
 chomp @repIDs; #remove enter key at each item
 warn "\n=> $id: Total $total_rep replicons found.\n";
 ++$counter;
#analysis for each replicon within the gbk file
 my $cc=0;
 foreach my $replicon (@repIDs)
 {
 my $repid=$id.'_'.$replicon;
 my ($outfaa,$outhmm,$outtxt)=("$repid.faa","$repid.hmmout","$repid.hmmout.txt");
 warn '-'x60,"\n[$counter/$total_num - ",++$cc,"/$total_rep] $repid ...\n";
 if(-e $outfaa) #already exists
 { warn "- OK. FASTA file already exists!\n";}
 else
 {!system("$gbk2seq $outgbk $replicon FAA > $outfaa") or die "$gbk2seq run with error";}
 if(-s $outfaa)
 {
  !system("$hmmsearch -E 1e-5 -o $outhmm $query $outfaa") or die "$hmmsearch run with error";
  system("$hs2tab $outhmm | $filter - > $outtxt"); #without check since it can be empty
  if(-s $outtxt) #not empty
  {push @allout, $outtxt;}
  else
  {
   warn "- None HMM found, skipped!\n";
#   unlink $outfaa,$outhmm,$outtxt;
   unlink $outhmm,$outtxt; #revised to keep faa file for further use, 2018/07/26
  }
 }
 else #no CDS, skipped
 {
  unlink $outfaa; 
#  die "error: no CDS found in $repid";
 }
 }
 unlink $outgbk;
}
warn "=> Done with all genomes. Total ",scalar @allout," of $total_num valid output files.\n\nSummarize the output now ...\n";

for(my $i=0;$i<=$#allout;$i+=10) #group to avoid long argument list
{
my $last=$i+9<$#allout?$i+9:$#allout;
my $input=sprintf join(' ',@allout[$i..$last]);
my $cmd="$summary $input ".($i?"| grep -v '#Genome' >":'')."> $output";
#!system($cmd) or die "error during run command: $cmd"; # the results might be null, can't check
system($cmd);
warn $last+1," ...\n" if ($last+1)%100==0 or $last==$#allout;
}
warn "All done.\n";
