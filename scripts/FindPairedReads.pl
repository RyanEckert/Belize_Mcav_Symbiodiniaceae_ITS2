#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Finds valid paired reads in a pair of files supplied by the user (FASTQ). 
Paired reads are written in the same order in a pair of files named
after the input files, with a prefix of 'paired.'
If the option "b" is chosen, orphans are also written to files named
after the input, prefixed with 'orphan.'
Usage: $scriptname -f F_reads -r R_reads <OPTIONS>
Required arguments:
	-f F_reads	input, forward reads (R1)(FASTQ)
	-r R_reads	input, reverse reads (R2)(FASTQ)
Options:
	-o option	choose 'p' to print only paired reads (default),
			or 'b' to print both paired and orphan reads.
	-g F_pair	a name for the paired forward (R1) reads output file (default: paired.input)	
	-i R_pair	a name for the paired reverse (R2) reads output file (default: paired.input)	
	-j F_orph	a name for the orphan forward (R1) reads output file (default: orphan.input)	
	-k R_orph	a name for the orphan reverse (R2) reads output file (default: orphan.input)	
USAGE
if ($#ARGV < 2 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('f:r:o:g:i:j:k:h');	
if (!$opt_f || !$opt_r || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
$fvar = $opt_f;		
$rvar = $opt_r;		
if ($opt_o)
	{
	if ($opt_o ne "p" && $opt_o ne "b") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
	$ovar = $opt_o;
	}
else {$ovar = "p";}

# rename variables and name output files
$if1 = $fvar;
$fvar =~ s/.+\///g;
$if2 = $rvar;
$rvar =~ s/.+\///g;
if ($opt_g) {$of1 = $opt_g;} else {$of1 = "paired.$fvar";}
if ($opt_i) {$of2 = $opt_i;} else {$of2 = "paired.$rvar";}
if ($opt_j) {$of3 = $opt_j;} else {$of3 = "orphan.$fvar";}
if ($opt_k) {$of4 = $opt_k;} else {$of4 = "orphan.$rvar";}

# loop through infile 1 and build a hash of sequence IDs
open(IN, $if1);
while(<IN>)
	{
	chomp;
	if ($_ !~ /\S/) {next;}
	$count++;
	$_ =~ s/\s.+//;
	if ($count==4) {$count=0;next;}
	if ($count==1)
	 	{
		$f1c++;
		$h1{$_}++;
		}
	}
close(IN);

# loop through infile 2 and build a hash of sequence IDs found in both files
open(IN, $if2);
$count=0;
while(<IN>)
	{
	chomp;
	if ($_ !~ /\S/) {next;}
	$count++;
	$_ =~ s/\s.+//;
	if ($count==4) {$count=0;next;}
	if ($count==1)
	 	{
		$f2c++;
		if (exists($h1{$_}))
			{
			$bh{$_}++;
			}
		}
	}
close(IN);

# loop through infile 1 and write out paired sequences to outfile 1 and unpaired to outfile 3
open(IN, $if1);
open(OUT, ">$of1");
if ($ovar eq "b") {open(UPO, ">$of3");}
$count=0; $switch=0;
while(<IN>)
	{
	chomp;
	if ($_ !~ /\S/) {next;}
	$count++;
	$sid = $_;
	$sid =~ s/\s.+//;
	if ($count==4) {$count=0;}
	if ($count==1)
	 	{
		$switch=0;
		if (exists($bh{$sid}))
			{
			$gc++;
			$switch++;
			}
		else	{$upn1++;}
		}
	if ($switch>0)
		{
		print OUT $_, "\n";
		}
	else
		{
		if ($ovar eq "b") {print UPO $_, "\n";}
		}
	}
close(IN);
close(OUT);
if ($ovar eq "b") {close(UPO);}

# loop through infile 2 and write out paired sequences to outfile 2 and unpaired to outfile 4
open(IN, $if2);
open(OUT, ">$of2");
if ($ovar eq "b") {open(UPO, ">$of4");}
$count=0; $switch=0;
while(<IN>)
	{
	chomp;
	if ($_ !~ /\S/) {next;}
	$count++;
	$sid = $_;
	$sid =~ s/\s.+//;
	if ($count==4) {$count=0;}
	if ($count==1)
	 	{
		$switch=0;
		if (exists($bh{$sid}))
			{
			$switch++;
			}
		else	{$upn2++;}
		}
	if ($switch>0)
		{
		print OUT $_, "\n";
		}
	else
		{
		if ($ovar eq "b") {print UPO $_, "\n";}
		}
	}
close(IN);
close(OUT);
if ($ovar eq "b") {close(UPO);}

print $f1c, " sequences in ", $if1, "\n";
print $f2c, " sequences in ", $if2, "\n";
print $gc, " valid pairs written to $of1 and $of2\n";
if ($ovar eq "b") 
	{
	print $upn1, " unpaired reads written to $of3\n";
	print $upn2, " unpaired reads written to $of4\n";
	}
