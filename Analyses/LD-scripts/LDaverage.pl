#!/usr/bin/perl -w
use warnings;
use strict;
#################################################################
###                                                           ###
###    Syntax: LDaverage.pl <input_filename> <bin-size>       ###
###                                                           ###
#################################################################
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nmissing parameter(s)\nUsage: $0 <input-filename> <bin-size(bp)>\n\n";
    exit;
}

use Data::Dumper qw(Dumper);
my $filename=$ARGV[0];
my $wSize=$ARGV[1];
my %Rbin;
my %Dbin;
my %LDvalue;
my $counter;
my $cSize;
my $maxDist = 100000; # max distance of markers to be considered
my $datestring = localtime();
open(FH, '<', $filename) or die $!;

my $firstLine = 1;
while(<FH>){
	if ($firstLine == 1) { $firstLine = 0;}
	 else {
	   my @fields = split /\s+/, $_;
	   $cSize = int(($fields[4] - $fields[1])/$wSize)*$wSize;
	   push ( @{$Rbin{ $cSize } }, $fields[6]);
	   push ( @{$Dbin{ $cSize } }, $fields[7]);
	}
}
close(FH);
my ($Rsum, $Dsum, $i, $Ravg,$Davg);
print "##############################################\n# inputfile: $filename\n# bin-size: ${wSize}bp\n";
print "# $datestring\n##############################################\n\n";
print "binStart\tbinEnd\tR-square\tD-value\tmrk-pairs\n";
while (my($index,@elem) = each %Rbin) {
	$Rsum = 0; $i = 0;
	foreach ( @{$Rbin{ $index } }) {
		$Rsum += $_;
		$i++ ;
	 }
	$Dsum = 0; $i = 0;
	foreach ( @{$Dbin{ $index } }) {
		$Dsum += $_;
		$i++ ;
	 }
	 $Ravg = $Rsum/$i;
	 $Davg = $Dsum/$i;
	 my $end = $index +$wSize -1;
	print "$index\t$end\t$Ravg\t$Davg\t$i\n";
}