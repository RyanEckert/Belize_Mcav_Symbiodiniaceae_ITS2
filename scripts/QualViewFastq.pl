#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Shows descriptive statistics for quality scores by position in a file of 
DNA sequences provided by the user. The output is a box plot in text format,
following the same conventions as the html output from FASTQC. 
Usage:   $scriptname -i reads <OPTIONS> 
Required arguments:
	-i reads	the FASTQ file to be analyzed (the script expects reads of equal length)

Options:
	-w <NUMBER>	window	size of the region shown (basepairs); longer reads will be shown in 
			blocks of this size (default = 100)
USAGE
if (!$ARGV[0]|| $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;
$mod1="Bio::SeqIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SeqIO;
$mod1="Statistics::Descriptive";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Statistics::Descriptive;

use utf8;
binmode STDOUT, ":utf8";

# get variables from input
getopts('i:w:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
$ivar = $opt_i;		# define variable from required argument opt_a
if ($opt_w) {$wvar = $opt_w;} else {$wvar = 100;}

# rename inputs
my $sfile = $ivar;
my $winsize = $wvar;

# read in sequences and scores, store quality scores for each position in a hash of arrays
my $seqs = Bio::SeqIO->new(-file=>$sfile, -format=>'fastq-sanger');
while ($seq = $seqs->next_seq)
	{
	$readlen = $seq->length();
	$seqcount++;
	@qsa = @{$seq->qual};
	$loc = 0;
	foreach $b (@qsa)
		{
		$loc++;
		if (exists($qh{$loc}))
			{
			push @{$qh{$loc}}, $b;
			}
		else
			{
			$qh{$loc} = $b;
			}
		}
	}

# calculate descriptive statistics for each array of quality scores
foreach $b (sort{$a<=>$b}(keys(%qh)))
	{
#	print "Position ", $b, "\n";
#	print join(" ", @{$qh{$b}}), "\n";
	$stat = Statistics::Descriptive::Full->new();
	$stat->clear();
	$stat->add_data(@{$qh{$b}});
	$meani = $stat->mean;
	$medi = $stat->median;
	$p25 = $stat->quantile(1);
	$p75 = $stat->quantile(3);
	$p10 = $stat->percentile(10);
	$p90 = $stat->percentile(90);
	$stath{$b}{"median"} = $medi;
	$stath{$b}{"mean"} = $meani;
	$stath{$b}{"p10"} = $p10;
	$stath{$b}{"p25"} = $p25;
	$stath{$b}{"p75"} = $p75;
	$stath{$b}{"p90"} = $p90;
#	print join(", ", ($meani, $medi, $p10, $p25, $p75, $p90)), "\n\n";
	}

# print a text boxplot to visualize these statistics
$vbchar = '│';	# unicode character for side of plot
$hbchar = '─';	# unicode character for top and bottom of plot
$tcchar = '┌';	# unicode character for top corner
$bcchar = '└';	# unicode character for bottom corner
$boxchar = '║';	# unicode character for box in boxplot
if ($winsize>$readlen) {$winsize = $readlen;}
$readlen = @qsa;
$steps = int($readlen/$winsize+(($winsize-1)/$winsize));
if ($steps==0) {$steps = 1;}
@sa = (1..$steps);
print "\n";
# print descriptive statistics for each region of X basepairs (default = 100)
foreach $s (@sa)
	{
# print the top axis and labels
	@qa = reverse(15..41);
#	print "@qa\n";
	print "\t", $tcchar, $hbchar x ($winsize-1), "\n";
	foreach $q (@qa)
		{
# print the left border and axis labels		
		if ($q==29) {print " e\t", $vbchar;}
		elsif ($q==33) {print " S  \t", $vbchar;}
		elsif ($q==32) {print " c\t", $vbchar;}
		elsif ($q==31) {print " o\t", $vbchar;}
		elsif ($q==30) {print " r  30\t", $vbchar;}
		elsif ($q==20) {print "    20\t", $vbchar;}
		elsif ($q==40) {print "    40\t", $vbchar;}
		else {print "\t", $vbchar;}
# print the graphical representation of quality score statistics
		@locvec = (1+($s-1)*$winsize .. ($s-1)*$winsize+$winsize);
		foreach $p (@locvec)
			{
#			print $stath{$p}{"p10"}, " ", $stath{$p}{"p90"}, ", ";
			if ($q <= $stath{$p}{"p90"} & $q >= $stath{$p}{"p10"})
				{
				if (($q+1) > $stath{$p}{"median"} & $q >= $stath{$p}{"median"} & ($q-1) < $stath{$p}{"median"})
					{
					print "-";
					}
				elsif ($q <= $stath{$p}{"p75"} & $q >= $stath{$p}{"p25"})
					{
					print $boxchar;
					}
				else 
					{
					print "|";
					}
				}
			else
				{
				print " ";
				}
			}
		print "\n";
		}

# print the bottom axis and labels
	print "    \t", $vbchar, "_" x ($winsize-1), "\n";
	print "\t";
	for ($x=0; $x<=$winsize; $x+=25)
		{
		print $x+($s-1)*$winsize, " ";
		if ($x<$winsize) {print " " x 22;}
		} 
	print "\n";
	print " " x int($winsize/2), "Position (basepairs)\n\n"; 
	print "\t", "'-': median, '║': 25% & 75% percentiles, and '│': 10% & 90% percentiles.\n\n";
	}


