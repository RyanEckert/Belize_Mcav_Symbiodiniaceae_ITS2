#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Filters a set of short reads in FASTQ format, excluding any reads that lack the 
expected sequence at the 5' end of the read. This is typically used to identify
valid amplicons in amplicon sequence analysis. 
Usage: $scriptname -i input -s sequence -m mismatches <OPTIONS>
Required arguments:
	-i input	name of the input file, FASTQ
	-s sequence	expected sequence (primer sequence)
	-m mismatches	maximum number of mismatches allowed between expected and observed sequences	
Options:
	-p pass.name	a name for the output file of sequences passing this filter (default: target.input)
	-f fail.name	a name for the output file of sequences failing this filter (default: nontarget.input)
USAGE
if ($#ARGV < 3 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:s:m:p:f:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_s || !$opt_m || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_p) {$pvar = $opt_p;} else {$pvar = "target.$opt_i";}	
if ($opt_f) {$fvar = $opt_f;} else {$fvar = "nontarget.$opt_i";}	

$seqfile = $opt_i;	# input
$pstr = $opt_s;		# primer string
$mmm = $opt_m;		# max number of mismatches allowed to still count as a match
$outfile = $pvar;	# name for the output file containing reads matching the primer
$badfile = $fvar;	# name for the output file containing reads NOT matching the primer

# -- loop through fastq file
# -- record whether each read matches the primer sequence 
# -- print reads matching the primer to output file
open (IN, $seqfile);
open (OUT, ">$outfile");
open (BAD, ">$badfile");
$status = 0;
while(<IN>)
	{
	chomp;
	$count++; 
	if ($count eq 1) 
		{
		$seqid = $_;
		$readsin++;
		}
	elsif ($count eq 2)
		{
		# screen it and record
		@pa = split("", $pstr);
		$npa = @pa;
		$target = substr $_, 0, $npa;
		@rpa = split("", $target);
		$mm = 0; $bcount = 0;
#		print $target, "\n";
#		print "@rpa\n";
#		print "@pa\n";
		foreach $b (@rpa)
			{
			$bcount++;
			if ($pa[$bcount-1] ne $rpa[$bcount-1]) {$mm++;}
			}
		if ($mm <= $mmm) 
			{
			$status = 1;
			$match++;
			}
#		print $mm, " mismatches\n";
		$seqstr = $_;
		}
	elsif ($count eq 4)
		{
		if ($status == 1) 
			{
			print OUT $seqid, "\n", $seqstr, "\n", "+\n", $_, "\n";
			}
		else 
			{
			print BAD $seqid, "\n", $seqstr, "\n", "+\n", $_, "\n";
			}
		$status = $count = 0;
		}
	}
close(IN);

print $readsin, " reads in input.\n";
print $match, " of these matched the primer sequence $pstr with $mmm or fewer mismatches.\n";
