#!/bin/bash/
#script to calculate the number of genomes covering each bHirRus1 chromosomes base
for file in coverage*.wig ; do 
    root=`basename $file .wig` 
	awk -v root="$root" 'BEGIN {print root}' 
    grep -v "fixedStep" $file | awk '{print $1}' | sort | uniq -c
done
