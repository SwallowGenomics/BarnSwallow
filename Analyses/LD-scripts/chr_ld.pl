#!/usr/bin/perl
#use warnings;
use strict;
my $filename = $ARGV[0]; # input filename
my $wSize = $ARGV[1]; # window size
my $mThreshold = 40; # Default minimum no. of marker-pairs for bin to be reported.
if ($ARGV[2]) {
	$mThreshold = $ARGV[2]; # Set minimum no. of marker-pairs for bin to be reported
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
		 #print "$pos\t$fields[1]\t$fields[4]\t$fields[6]\t$fields[7]\t$wSize\n";
		 $Dprime[$pos] +=  $fields[7];
		 $R2[$pos] +=   $fields[6];
		 $binSize[$pos]++;	 
	}
}
close(FH);
#print "Window size $wSize\n";
my $datestring = localtime();
print "##############################################\n# inputfile: $filename\n# bin-size: ${wSize}bp\n";
print "# minimum no. of marker-pairs: $mThreshold\n";
print "# $datestring\n##############################################\n\n";
print "chr\tstart\tend\tR-squared\tD-prime\tmrk-pairs\n";
while (my($index,$elem) = each @Dprime) {
	my $start = $index * $wSize;
	my $end = $start + $wSize-1;
	if ($elem) {
		my $Dvalue = $elem/$binSize[$index];
		my $Rsquared = $R2[$index]/$binSize[$index];
		#print "$index\t$elem\t$start\t$end\t$binSize[$index]\n";
		print "$chr\t$start\t$end\t$Rsquared\t$Dvalue\t$binSize[$index]\n" unless $binSize[$index] <= $mThreshold;
	}
	#print "$index\t$elem\t$binSize[$index]\n";
}
#print join ("\n", @LDvalue), "\n";
#print join ("\n", @binSize), "\n";
