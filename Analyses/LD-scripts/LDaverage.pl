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
my $cSize; ####### bin id = hash key
my $datestring = localtime();

open(FH, '<', $filename) or die $!;
my $firstLine = 1; ########## 1=header row in input file; 0=no header row
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
print "##############################################\n#  inputfile: $filename\n#  bin-size: ${wSize}bp\n";
print "#  date: $datestring\n##############################################\n\n";
print "binStart\tbinEnd\tR-square\tD-value\tmrk-pairs\n";

for my $key (sort { $Rbin{$a} <=> $Rbin{$b} } keys %Rbin) {
	$Rsum = 0;
	foreach ( @{$Rbin{ $key } }) {
		$Rsum += $_;
	 }
	$Dsum = 0;
	foreach ( @{$Dbin{ $key } }) {
		$Dsum += $_;
	 }
	 $i = scalar( @{ $Rbin{ $key } } );
	 $Ravg = $Rsum/$i;
	 $Davg = $Dsum/$i;
	 my $end = $key +$wSize -1;
	print "$key\t$end\t$Ravg\t$Davg\t$i\n";
}