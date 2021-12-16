#!/usr/bin/perl
#use warnings;
use strict;
######################################################################################
###                                                                                ###
###    Syntax: chr_ld.pl <input_filename> <window-size> [min-#-of-marker-pairs]    ###
###                                                                                ###
######################################################################################
my $filename = $ARGV[0]; # input filename
my $wSize = $ARGV[1]; # window size
my $mThreshold = 40; # Default minimum no. of marker-pairs to compute window
if ($ARGV[2]) {
	$mThreshold = $ARGV[2]; # Set minimum # of marker-pairs to compute window
}

#  ======= check command line parameters ===================
if (not defined $filename) {
	 die "\nMissing input filename ! \n\nsyntax: $0 <input_file_name> <window-size> [min-#-marker-pairs]\n";
	} else {
	die "\nThe file $filename does not exist, and I can't go on without it.\n\nsyntax: $0 <input_file_name> <window-size> [min-#-marker-pairs]\n" unless -e $filename;
}
if (not defined $wSize) {
	 die "Missing window-size parameter! \n";
	}
# ==========================================================
my @Dprime =();
my @R2 =();
my @binSize =();
my @field =();
my $pos;
my $chr = "";
my $maxDist = 100000; # max distance between markers to be considered
open(FH, '<', $filename) or die $!;

while(<FH>){
   my @fields = split /\s+/, $_;
   $pos = int($fields[1]/$wSize);
	$chr = $fields[0];
	if (int($fields[4]/$wSize) == $pos) {
		 $Dprime[$pos] +=  $fields[7];
		 $R2[$pos] +=   $fields[6];
		 $binSize[$pos]++;	 
	}
}
close(FH);
my $datestring = localtime();
print "##############################################\n# inputfile: $filename\n# window-size: ${wSize}bp\n";
print "# minimum no. of marker-pairs: $mThreshold\n";
print "# $datestring\n##############################################\n\n";
print "chr\tstart\tend\tR-squared\tD-prime\tmrk-pairs\n";
while (my($index,$elem) = each @Dprime) {
	my $start = $index * $wSize;
	my $end = $start + $wSize-1;
	if ($elem) {
		my $Dvalue = $elem/$binSize[$index];
		my $Rsquared = $R2[$index]/$binSize[$index];
		print "$chr\t$start\t$end\t$Rsquared\t$Dvalue\t$binSize[$index]\n" unless $binSize[$index] <= $mThreshold;
	}

}
