#!/usr/bin/perl -w
use warnings;
use strict;
my $num_args = $#ARGV + 1;
if ($num_args != 3) {
    print "\nmissing parameter(s)\nUsage: $0 <input_filename> <chr_size_filename> <bin-size(bp)>\n\n";
    exit;
}

use Data::Dumper qw(Dumper);
my $filename=$ARGV[0]; # input file produced with plink containing LD estimates
	unless(-e $filename) {  die "$filename File Does not Exist!\n $! "; }
my $chrLenFile=$ARGV[1]; # input file containing names and size of chromosomes
	unless(-e $chrLenFile) {  die "$chrLenFile File Does not Exist!\n $!"; }
my $binSize=int($ARGV[2]); # Size of the bin (window) where average is calculated
	if ($binSize < 1 ) { print STDERR "binsize $binSize is not valid!\n $!"; exit;}

my %chrLen; # contain the size of all chrmosomes
my %Rbin; # stores the average R2 value per bin
my %Dbin; # stores the average D'-values value per bin
my $binName; # bin name calculated according to SNP position 
my $datestring = localtime;
my $Dist1; # distance from the distal telomere
my $snpAvgPos; #average position of SNP pairs
my @fields;
if (open(LOG, '>', $filename.".log")){
	print LOG "##############################################\n# scriptfile: $0\n# inputfile: $filename\n# chr-size file: $chrLenFile\n# bin-size: ${binSize}bp\n";
	print LOG "# $datestring\n##############################################\n\n";
	print LOG "Chromosome lenghts as read from: $chrLenFile\n";
	if (open(CH, '<', $chrLenFile)){
		while(<CH>){
		     my @fields = split /\s+/, $_;
		     $chrLen{ $fields[0] } = eval(int($fields[1])); 
		     print LOG "$fields[0]: ".int($fields[1])."\n";
		}
		print LOG "\n";
		close(CH);
	} else { die print STDERR "Error: Couldn't open the file - $chrLenFile\n";}
	# while ((my $key, my $value) = each (%chrLen))
	# {
	  # do whatever you want with $key and $value here ...
	#   $value = $chrLen{$key};
	#   print "  $key =====> $value\n";
	# }
	# exit;
	   
	if (open(FH, '<', $filename)) {
		my $firstLine = 1;
		while(<FH>){
			if ($firstLine == 1) { $firstLine = 0;}
			 else {
			   @fields = split /\s+/, $_;
			   $snpAvgPos = int((($fields[4] - $fields[1])/2)+$fields[1]); #  average SNP pair position ==> distance from chromosome start
			   if (defined $chrLen{$fields[0]}) {
			   	 $Dist1 = $chrLen{ $fields[0]} - $snpAvgPos;  # distance from chromosome end
			   } else {
			   	 print LOG "$fields[0] chromosome not found on line $. in file $filename!\n";
			   	 print STDERR "$fields[0] chromosome not found on line $. in file $filename!\n";
			   	 $Dist1 = "9999999999999"; # very large chromosome size in case chr name is not present in che_len input file
		        }
		       if ($snpAvgPos < $Dist1) {
		            eval ($binName = abs(int($snpAvgPos/$binSize))*$binSize);  # bin name when SNP are closer to chr start
		            if ($@) {print LOG "Error in the folling line:\n$_\n$@\n";}
		       } else {
		            eval($binName = abs(int($Dist1/$binSize))*$binSize); # bin name when SNP are closer to chr end
		            if ($@) {print LOG "Error in the folling line:\n$_\n$@\n";}
		       }
			   push ( @{$Rbin{ $binName } }, $fields[6]);
		       push ( @{$Dbin{ $binName } }, $fields[7]);
			}
		}
		close(FH);
	} else { die print STDERR "Error: Couldn't open the file - $filename\n";}
	# ============= PRINT OUTPUT ===================
	my ($sum,  $i, $avgR,$avgD);
	print "##############################################\n# scriptfile: $0\n# inputfile: $filename\n# chr-size file: $chrLenFile\n# bin-size: ${binSize}bp\n";
	print "# $datestring\n##############################################\n\n";
	print "binStart\tbinEnd\tR2\tD-value\tmrk-pairs\n";
	while (my($index,@elem) = each %Rbin) {
		$sum = 0; $i = 0;
		foreach ( @{$Rbin{ $index } }) {
			$sum += $_;
			$i++ ;
		 }
		 $avgR = $sum/$i;
	     $sum = 0; $i = 0;
	     foreach ( @{$Dbin{ $index } }) {
			$sum += $_;
			$i++ ;
		 }
	     my $avgD = $sum/$i;
		 my $end = $index +$binSize -1;
		print "$index\t$end\t$avgR\t$avgD\t$i\n";
	}
} else { die print STDERR "Error: Couldn't create log file\n";}
close (LOG);
exit;